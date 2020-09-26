function df = ephys_tools_rnr_wrapper()
% get sessions to run
df = readtable('F:\Projects\PAE_PlaceCell\swr_data\post_processed\swr_df.csv');
data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
save_path = 'F:\Projects\PAE_PlaceCell\replay\';

% run through each session
WaitMessage = parfor_wait(length(unique(df.session)),'Waitbar',false,'ReportInterval',1);
sessions = unique(df.session);
for i = 1:length(sessions)
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


for i = find(df.pvalue<0.05)'

    figure;
    subplot(2,1,1)
    imagesc(df.Pr{i}');
    colormap magma
    hold on;
    [~,I] = max(df.Pr{i},[],2);
    plot(1:size(df.Pr{i},1),I,'w')
    title([df.ep_type(i), num2str(df.bayesLinearWeighted(i))])
    pause(0.00001)
    
    subplot(2,1,2)
    data = load([data_path,df.session{i},'.mat'],'Spikes','ratemap');
    
    % template (mean both directions together)
    template = (vertcat(data.ratemap{:,1}) + vertcat(data.ratemap{:,2})) ./ 2;
    % include place cells
    [include,peak_loc] = get_place_cells(template);
    [~,idx] = sort(peak_loc,'descend');
    data.Spikes = data.Spikes(idx);
    
    spikes = format_spikes(data);
    idx = spikes.spindices(:,1) >= df.start(i) & spikes.spindices(:,1) <= df.stop(i);
    idx_include = ismember(spikes.spindices(:,2),include);
    scatter(spikes.spindices(idx & idx_include,1),spikes.spindices(idx & idx_include,2),...
        20,spikes.spindices(idx & idx_include,2),'filled')
    xlim([df.start(i),df.stop(i)])
    ylim([1,max(spikes.spindices(:,2))])
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
% template (mean both directions together)
template = (vertcat(data.ratemap{:,1}) + vertcat(data.ratemap{:,2})) ./ 2;

% include place cells
[include,peak_loc] = get_place_cells(data);           
if isnan(include) % if no place cells
    replayScores = NaN;
    return
elseif length(include) < 6 % at least 6 place cells
    replayScores = NaN;
    return
end

% % get candidate events based on Wu and Foster 2014 (spike density)
% ripples = get_candidate_events(data,spikes,include);
% if isempty(ripples.peaks)
%     replayScores = NaN;
%     return
% end

% Estimates the Bayesian method for replay quantification
[replayScores] = replay_Bayesian(spikes,ripples,template,include);
replayScores.peak_loc = peak_loc;
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
% expand window until event drops below mean (forwards and backwords)
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

function [include,peak_loc] = get_place_cells(data)
% % read in list of place cells
% df = readtable('D:\ryanh\github\harvey_et_al_2020\Rdata_pae_track_cylinder.csv');
% % get current session and track data
% df = df(strcmp(df.session,[data.rat,'_',data.sessionID,'.mat']) & df.mazetype == "track",:);
% % get index of place cells to include
% if ~isempty(df.tt)
%     include = find_cells(data,str2double(extractBetween(df.tt,'TT','.mat')),(df.cell));
% else
%     include = NaN;
% end


% straighten track
% [~,data.linear_track{1}.nonlinearFrames(:,2:3),~] =...
%     pca([data.linear_track{1}.nonlinearFrames(:,2),...
%      data.linear_track{1}.nonlinearFrames(:,3)]);
%  
% XY=interp1(data.linear_track{1}.nonlinearFrames(:,1),...
%     data.linear_track{1}.nonlinearFrames(:,2:3),spikes.spindices(:,1),'linear');
data.ratemap
for i = 1:size(data.ratemap,1)

    fields = place_cell_analysis.getPlaceFields('ratemap',template(i,:)');
    % check if field width is entire track
    place_field(i) = fields{1, 1}{1, 1}.width ~= size(template,2);
    % check number of fields
    n_fields(i) = length(fields{1, 1});
    % peak location
    peak_loc(i) = fields{1, 1}{1, 1}.peakLoc;
end
include = find(place_field & n_fields == 1);
end

function df=construct_df(replayScores)
    df = replayScores.df;
    df.bayesLinearWeighted = replayScores.bayesLinearWeighted(:);
    df.bayesRadon = replayScores.bayesRadon(:);
    df.slope_hpc = replayScores.slope_hpc(:);
    df.pvalue = replayScores.pvalue(:);
    df.z = replayScores.z(:);

%     df.bayesLinearWeighted_cellID_shuf = replayScores.bayesLinearWeighted_cellID_shuf(:);
%     df.bayesRadon_cellID_shuf = replayScores.bayesRadon_cellID_shuf(:);
%     df.bayesLinearWeighted_circular_shuf = replayScores.bayesLinearWeighted_circular_shuf(:);
%     df.bayesRadon_circular_shuf = replayScores.bayesRadon_circular_shuf(:);
    df.nCells = replayScores.nCells(:);
    df.nSpks = replayScores.nSpks(:);
    df.Pr = replayScores.Pr(:);
%     df.above_shuff = replayScores.above_shuff(:);
    df.start = replayScores.start(:);
    df.stop = replayScores.stop(:);

    idx = isnan(replayScores.bayesLinearWeighted);
    df(idx,:) = [];
end