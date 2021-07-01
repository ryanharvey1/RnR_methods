% decoding_accuracy


data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
df = readtable('F:\Projects\PAE_PlaceCell\replay\processed\replay_df.csv');
% Pr = load('F:\Projects\PAE_PlaceCell\replay\processed\replay_Pr.mat','Pr');
save_path = 'F:\Projects\PAE_PlaceCell\decoding_accuracy\';
% run through each session
WaitMessage = parfor_wait(length(unique(df.session)),'Waitbar',false,'ReportInterval',1);
sessions = unique(df.session);
parfor i = 1:length(sessions)
    if exist([save_path,sessions{i},'.mat'],'file')
        continue
    end
    data = load([data_path,sessions{i},'.mat'],'Spikes','events','ratemap',...
        'rat','sessionID','spikesID','linear_track','samplerate','frames');
    
    dat = decoded_position(data);
    
    save_data(save_path,sessions{i},dat)
    
    WaitMessage.Send;
end
WaitMessage.Destroy

df_dat = load_all(save_path);

% assign group
[~,b,~] = unique(df.session);
df_dat.group = df.group(b);

% save table
df_to_save = df_dat;
df_to_save.x_obs = [];
df_to_save.x_decoded = [];
df_to_save.mdl = [];

mkdir([save_path,'processed'])

% save df
writetable(df_to_save,[save_path,'processed\decoding_accuracy.csv'])

% quick compare
[p,h,stats] = ranksum(df_dat.r2(df_dat.group == "control"),df_dat.r2(df_dat.group == "pae"))
figure;
boxplot(df_dat.r2,df_dat.group)

function save_data(save_path,session,dat)
save([save_path,session,'.mat'],'dat')
end

function df = load_all(save_path)
df = table();
file = dir([save_path,'*.mat']);
for i = 1:length(file)
    load([save_path,file(i).name],'dat');
    if ~isstruct(dat)
        continue
    end
    df=[df;construct_df(dat,file(i).name)];
end
end

function df=construct_df(dat,session)
df = table();
df.session = {session};
df.r2 = dat.r2;
df.intercept = dat.intercept;
df.slope = dat.slope;
df.n_cells = dat.n_cells;
df.direction = dat.direction;
df.session_length = dat.session_length;

df.x_obs{1} = dat.x_obs;
df.x_decoded{1} = dat.x_decoded;
df.mdl{1} = dat.mdl;
end

% below are functions to make each panel in the plot
function dat = decoded_position(data)

[frames] = restrict_track(data);
frames(:,2) = rescale(frames(:,2),1,40);

% find direction with most place fields
for d = 1:2
    [include,~,~] = get_place_cells(data,d);
    n_cells(d) = sum(include);
end
[~,d] = max(n_cells);
[include,~,template] = get_place_cells(data,d);

% limit spikes to linear track
start = data.events(1,1);
stop = data.events(2,1);

for i = 1:length(data.Spikes)
    rip_idx = data.Spikes{i} >= start &...
        data.Spikes{i} <= stop;
    temp_spikes{i} = data.Spikes{i}(rip_idx);
end
binSize = .2;
overlap = 4;
spkmat = bz_SpktToSpkmat(temp_spikes(include),'overlap',overlap,'binSize',binSize * overlap);

idx = spkmat.timestamps >= start & spkmat.timestamps <= stop;
spkmat.data = spkmat.data(idx,:);
spkmat.timestamps = spkmat.timestamps(idx);
x_obs = interp1(frames(:,1),frames(:,2),spkmat.timestamps);

% generate posterior probability matrix using template and event FRs
[Pr, x_decoded] = placeBayes(spkmat.data, template(include,:), spkmat.dt); 

w = max(Pr,[],2);

mdl = fitlm(x_obs,x_decoded,'Weights',w);
r2 = mdl.Rsquared.Ordinary;
intercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);

dat.x_obs = x_obs;
dat.x_decoded = x_decoded;
dat.mdl = mdl;
dat.r2 = r2;
dat.intercept = intercept;
dat.slope = slope;
dat.n_cells = sum(include);
dat.direction = d;
dat.session_length = stop - start;

% figure;
% scatter(x_obs,x_decoded)
% hold on
% plot(x_obs,x_obs*slope + intercept)
% 
% 
% real(sqrt(x_obs.^2 - x_decoded.^2));
% 
% 
% 
% h = imagesc(Pr');
% axis xy
% colormap(gca,(magma(256)));
% hold on
% plot(x_obs,'w');
% xlim([1,length(x_obs)])
% box off
% set(gca,'xtick',[])
% ylabel('Position')
% 
% originalSize1 = get(gca, 'Position');
% cb = colorbar('box','off');
% set(gca, 'Position', originalSize1);
% cb.Title.String = 'Pr';
% cb.Position(3) = 0.5*cb.Position(3);

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
include = (place_field & (n_fields == 1));
template = vertcat(data.ratemap{:,d});
end
function [frames] = restrict_track(data)
% for i = 1:length(data.Spikes)
%     data.Spikes{i} = data.Spikes{i}(data.Spikes{i} >= data.events(1,1) &...
%         data.Spikes{i} <= data.events(2,1));
% end
frames = data.frames(data.frames(:,1) >= data.events(1,1) &...
    data.frames(:,1) <= data.events(2,1),:);
end
