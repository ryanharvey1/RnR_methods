function df = ephys_tools_rnr_wrapper()
% get sessions to run
df = readtable('F:\Projects\PAE_PlaceCell\swr_data\post_processed\swr_df.csv');
data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
save_path = 'F:\Projects\PAE_PlaceCell\replay\';

% run through each session
WaitMessage = parfor_wait(length(unique(df.session)),'Waitbar',false,'ReportInterval',1);
sessions = unique(df.session);
parfor i = 1:length(sessions)
    if exist([save_path,sessions{i},'.mat'],'file')
        continue
    end
    data = load([data_path,sessions{i},'.mat'],'Spikes','events','ratemap',...
        'rat','sessionID','spikesID','linear_track','samplerate','frames');
    
    replayScores = main(data);
    
    save_data(save_path,sessions{i},replayScores)
    
    WaitMessage.Send;
end
WaitMessage.Destroy

df = load_all(save_path);

for i = find(abs(df.bayesLinearWeighted) > 0.6)'
    data = load([data_path,df.session{i},'.mat'],'Spikes','ratemap','frames','events');
    track_idx = data.frames(:,1) >= data.events(1,1) & data.frames(:,1) <= data.events(2,1);
    data.frames(track_idx,2) = rescale(data.frames(track_idx,2),1,40);
    idx = data.frames(:,1) >= df.start(i) & data.frames(:,1) <= df.stop(i);    

    figure;
    subplot(2,2,2)
    imagesc(df.Pr{i}');
    colormap magma
    hold on;
    [~,I] = max(df.Pr{i},[],2);
    plot(1:size(df.Pr{i},1),I,'w')
    
    plot(rescale(data.frames(idx,1),1,size(df.Pr{i},1)),data.frames(idx,2),'r')
    title([df.ep_type(i), num2str(df.bayesLinearWeighted(i))])
    
    
    subplot(2,2,4)
    % include place cells
    [include,~,template] = get_place_cells(data,1);  
    
    maps = template(include,:);
    [~,I] = max(maps,[],2);
    [~,I] = sort(I,'descend');
    maps = maps(I,:);
         
    data.Spikes = data.Spikes(include);
    data.Spikes = data.Spikes(I);

    spikes = format_spikes(data);
    
    idx = spikes.spindices(:,1) >= df.start(i) & spikes.spindices(:,1) <= df.stop(i);
    
    scatter(spikes.spindices(idx,1),spikes.spindices(idx,2),...
        20,spikes.spindices(idx,2),'filled')
    
    xlim([df.start(i),df.stop(i)])
    ylim([1,max(spikes.spindices(:,2))])
    
    subplot(2,2,3)

    for t = 1:size(maps,1)
        x = zscore(maps(t,:)) + ones(1,size(template,2))+t*2.25;
        plot(x,'linewidth',3)
        hold on;
    end
    darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
    pause(0.00001)
end

end

function save_data(save_path,session,replayScores)
save([save_path,session,'.mat'],'replayScores')
end

function df = load_all(save_path)
df = table();
file = dir([save_path,'*.mat']);
for i = 1:length(file)
    load([save_path,file(i).name],'replayScores');
    if ~isstruct(replayScores)
        continue
    end
   df=[df;construct_df(replayScores)];
end
end

function replayScores = main(data)
% format inputs from ephys_tools standards
% restrict to linear track and pedestals (pre/post)
data = restrict_track_pedestals(data);
% spikes
spikes = format_spikes(data);

%%
% read in and format ripples
[ripples,cur_df] = get_ripples(data);
if isempty(ripples.peaks)
    replayScores = NaN;
    return
end

% Estimates the Bayesian method for replay quantification
[replayScores] = replay_Bayesian(data,spikes,ripples);
% replayScores.peak_loc = peak_loc;
replayScores.df = cur_df;
% compareReplayMethods(spikes,ripples,template,include)
end

function data = restrict_track_pedestals(data)
if size(data.events,2) > 1
    stop = data.events(1,2);
else
    stop = data.events(2,1);
end
for i = 1:length(data.Spikes)
    data.Spikes{i} = data.Spikes{i}(data.Spikes{i} <= stop);
end
end

function spikes = format_spikes(data)
% spikes
spike_lengths = cellfun('length',data.Spikes);
uid = [];
for i = 1:length(data.Spikes)
    uid = [uid;repmat(i,spike_lengths(i),1)];
end
[spikes_times,idx] = sort(vertcat(data.Spikes{:}));
spikes.UID = uid(idx);
spikes.spindices = [spikes_times,spikes.UID];
spikes.times = data.Spikes';
end

function [ripples,cur_df] = get_ripples(data)
df = readtable('F:\Projects\PAE_PlaceCell\swr_data\post_processed\swr_df.csv');
idx = strcmp(df.ep_type,'track') | strcmp(df.ep_type,'pedestal_1') | strcmp(df.ep_type,'pedestal_2');
cur_df = df(strcmp(df.session,[data.rat,'_',data.sessionID]) & idx,:);
ripples.timestamps = [cur_df.start_time,cur_df.end_time];
ripples.peaks = cur_df.peak_time;
end

function ripples = get_candidate_events(data,spikes,include)
dt = 0.01;
% set up edges and bin centers
edges = min(spikes.spindices(:,1)):dt:max(spikes.spindices(:,1));
bin_centers = min(spikes.spindices(:,1))+dt/2:dt:max(spikes.spindices(:,1))-dt/2;
% spike density function of place cells
spf = histcounts(spikes.spindices(ismember(spikes.spindices(:,2),include),1),edges);
% smooth with a 15ms SD
spf = smoothdata(spf,'gaussian',3);
% find epochs where spike density function is > 2 SD
peaks = bin_centers(spf > mean(spf) + std(spf)*2);
% remove peaks in adjacent bins
peaks = peaks(find(diff(peaks) > dt+dt/2));

% restrict to epochs where speed is < 5cm/sec
speed = interp1(data.frames(:,1),data.frames(:,5),peaks);
peaks = peaks(speed < 5 | isnan(speed));
% expand window until event drops below mean (forwards and backwards)
for event = 1:length(peaks)
   idx = find(peaks(event) == bin_centers);
   for i = 1:length(bin_centers)
       if idx+i == length(bin_centers)
           timestamps(event,2) = idx+i;
           break
       end
       if spf(idx+i) < mean(spf)
           timestamps(event,2) = bin_centers(idx+i);
           break
       end
   end
   for i = 1:length(bin_centers)
       if idx-i < 1
           timestamps(event,1) = idx-i;
           break
       end
       if spf(idx-i) < mean(spf)
           timestamps(event,1) = bin_centers(idx-i);
           break
       end
   end
end
remove_idx = find(any(diff(timestamps)==[0,0],2));
peaks(remove_idx) = [];
timestamps(remove_idx,:) = [];

ripples.peaks = peaks;
ripples.timestamps = timestamps;
end

function [include,peak_loc,template] = get_place_cells(data,d)
% get_place_cells: pull rate map that has 1 place field and maximizes
% firing rate between the two running directions

% check for place field in both directions
for i = 1:size(data.ratemap,1)
        fields = place_cell_analysis.getPlaceFields('ratemap',data.ratemap{i,d}');
        % get peak rate
        peak_rate(i) = fields{1, 1}{1, 1}.peakFR;
        % check if field width is entire track
        place_field(i) = fields{1, 1}{1, 1}.width ~= size(data.ratemap{i,d},2);
        % check number of fields
        n_fields(i) = length(fields{1, 1});
        % peak location
        peak_loc(i) = fields{1, 1}{1, 1}.peakLoc;
end
% has one field 
include = find(place_field & (n_fields == 1));
template = vertcat(data.ratemap{:,d});
end

function df=construct_df(replayScores)
    df = replayScores.df;
    df.bayesLinearWeighted = replayScores.bayesLinearWeighted(:);
    df.bayesRadon = replayScores.bayesRadon(:);
    df.slope_hpc = replayScores.slope_hpc(:);
    df.pvalue_cellID_shuf = replayScores.pvalue_cellID_shuf(:);
    df.z_cellID_shuf = replayScores.z_cellID_shuf(:);
    df.pvalue_circular_shuf = replayScores.pvalue_circular_shuf(:);
    df.z_circular_shuf = replayScores.z_circular_shuf(:);

    df.nCells = replayScores.nCells(:);
    df.nSpks = replayScores.nSpks(:);
    df.Pr = replayScores.Pr(:);
    df.start = replayScores.start(:);
    df.stop = replayScores.stop(:);
    
    df.direction_used = replayScores.direction_used(:);

    idx = isnan(replayScores.bayesLinearWeighted);
    df(idx,:) = [];
end