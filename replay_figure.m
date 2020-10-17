data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
df = readtable('F:\Projects\PAE_PlaceCell\replay\processed\replay_df.csv');
Pr = load('F:\Projects\PAE_PlaceCell\replay\processed\replay_Pr.mat','Pr');

% df_temp = df(df.pvalue_cellID_shuf < 0.05,:);
% idx = [148,145,144,143,123,119,115,111,110,98,97,96,76,73,30,1];

df_temp = df(df.pvalue_circular_shuf < 0.05,:);
% idx = [441,440,429,426,425,424,423,422,414,412,390,385,384,381,376,375,373,...
%     367,366,361,360,358,357,355,354,351,348,345,335,334,333,329,304,278,273,...
%     270,259,258,241,226,216,211,168,38,35,33,6,4,2,1];


% ripple to show
ripples_to_use_idx = 23628;
idx = find(df_temp.ripple_number == ripples_to_use_idx);

% loop through potential replay events
for rip_n = idx
    data = load([data_path,df_temp.session{rip_n},'.mat']);
   
    plot_replay(data,df_temp,df,Pr,rip_n,df_temp.start(rip_n)-10,df_temp.stop(rip_n)+10)
end

% make plot
function plot_replay(data,df_temp,df,Pr_replay,rip_n,start,stop)
% set figure defaults
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','Serif','defaultTextFontName','Serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('thesis', 1, [4,3]);
set(fig,'Position',[835 270 fig_width_in fig_height_in])

% set replay window here (you can add a buffer if you want)
pr_start = df_temp.start(rip_n);
pr_stop = df_temp.stop(rip_n);

% create figure panels

subplot(3,4,[2.5,4]);
plot_decoded_position(data,df_temp,rip_n,start,stop)

subplot(3,4,[5.6,5.6])
plot_tuning_curves(data,df_temp,rip_n)

subplot(3,4,[6.5,8])
plot_raster(data,start,stop,df_temp,rip_n,pr_start,pr_stop)

subplot(3,4,[9,10])
plot_ripple_raster(data,start,stop,pr_start,pr_stop,df_temp,rip_n)
                       
subplot(3,4,[11,12]);
plot_posterior(Pr_replay,df,df_temp,rip_n)

% save pngs to folder (for final version, save to svg)
id = [df_temp.session{rip_n},'_',...
    df_temp.group{rip_n},'_',...
    df_temp.ep_type{rip_n},'_',...
    num2str(df_temp.ripple_number(rip_n))];
saveas(fig,['F:\Projects\PAE_PlaceCell\replay_temp_figs\',id,'.png'])
pause(0.00001)
close all
end

% below are functions to make each panel in the plot
function plot_decoded_position(data,df_temp,rip_n,start,stop)

[frames] = restrict_track(data);
frames(:,2) = rescale(frames(:,2),1,40);

[include,~,template] = get_place_cells(data,df_temp.direction_used(rip_n));

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

h = imagesc(Pr');
axis xy
colormap(gca,(magma(256)));
hold on
plot(x_obs,'w');
xlim([1,length(x_obs)])
box off
set(gca,'xtick',[])
ylabel('Position')

originalSize1 = get(gca, 'Position');
cb = colorbar('box','off');
set(gca, 'Position', originalSize1);
cb.Title.String = 'Pr';
cb.Position(3) = 0.5*cb.Position(3);

end

function plot_lfp(lfp,start,stop,df_temp,rip_n)
rip_idx = lfp.timestamps >= start &...
    lfp.timestamps <= stop;
x = lfp.timestamps(rip_idx);
y = lfp.data(rip_idx);
plot(x,y,'k')
hold on
% scatter(df_temp.peak_time(rip_n),max(y),'*r')
plot([df_temp.start_time(rip_n),df_temp.start_time(rip_n)],ylim,'r')
plot([df_temp.end_time(rip_n),df_temp.end_time(rip_n)],ylim,'r')
box off
axis tight
set(gca,'xtick',[])
xlim([start,stop])
end

function plot_raster(data,start,stop,df_temp,rip_n,pr_start,pr_stop)
for i = 1:length(data.Spikes)
    rip_idx = data.Spikes{i} >= start &...
        data.Spikes{i} <= stop;
    temp_spikes{i} = data.Spikes{i}(rip_idx)';
end

[include,~,template] = get_place_cells(data,df_temp.direction_used(rip_n));
maps = template(include,:);
[~,I] = max(maps,[],2);
[~,I] = sort(I);
maps = maps(I,:);
temp_spikes = temp_spikes(include);
temp_spikes = temp_spikes(I);

clear x y
for i=1:length(temp_spikes)
    if isempty(temp_spikes{i})
        x{i} = [];
        y{i} = [];
        continue
    end
    [x{i},y{i}]=plotSpikeRaster(temp_spikes(i),'PlotType','vertline');
end
cla(gca)
ax = gca;
ax.YDir = 'normal';
colors = colormap(gca,cool(length(temp_spikes)));

patch([pr_start-.10,pr_start-.10,pr_stop+.10,pr_stop+.10],...
    [.5,length(x)+.5,length(x)+.5,.5],'y','FaceAlpha',.5,'EdgeColor','none')

hold on
ylim([0.5,length(x)+.5])
for i=1:length(x)
    plot(x{i},y{i}+(i-1),'Color',colors(i,:),'linewidth',1);hold on
end
box off
xlim([start,stop])
xlabel('Time (s)')
ylabel('Unit #')
set(gca,'ytick',[])
end

function plot_ripple_raster(data,start,stop,pr_start,pr_stop,df_temp,rip_n)
for i = 1:length(data.Spikes)
    rip_idx = data.Spikes{i} >= start &...
        data.Spikes{i} <= stop;
    temp_spikes{i} = data.Spikes{i}(rip_idx)';
end

[include,~,template] = get_place_cells(data,df_temp.direction_used(rip_n));
maps = template(include,:);
[~,I] = max(maps,[],2);
[~,I] = sort(I);
maps = maps(I,:);
temp_spikes = temp_spikes(include);
temp_spikes = temp_spikes(I);

clear x y
for i=1:length(temp_spikes)
    if isempty(temp_spikes{i})
        x{i} = [];
        y{i} = [];
        continue
    end
    [x{i},y{i}]=plotSpikeRaster(temp_spikes(i),'PlotType','vertline');
end
colors = colormap(gca,cool(length(temp_spikes)));
xlim([pr_start,pr_stop])
for i=1:length(x)
    plot(x{i},y{i}+(i-1),'Color',colors(i,:),'linewidth',2);hold on
end
box off
xlabel('Time (s)')
ylabel('Unit #')
set(gca,'ytick',[])
xlim([pr_start,pr_stop])
end

function plot_tuning_curves(data,df_temp,rip_n)
[include,~,template] = get_place_cells(data,df_temp.direction_used(rip_n));
maps = template(include,:);
% maps = template;
[~,I] = max(maps,[],2);
[~,I] = sort(I);
maps = maps(I,:);
colors = colormap(gca,cool(size(maps,1)));

for t = 1:size(maps,1)
    x = rescale(maps(t,:),0.5,1.5) + t;
    plot(x,'Color',colors(t,:),'linewidth',1)
    hold on;
end
box off
set(gca,'ytick',[])
xlabel('Position')
ylabel('Unit #')
axis tight
end

function plot_posterior(Pr_replay,df,df_temp,rip_n)
map = Pr_replay.Pr{df.ripple_number == df_temp.ripple_number(rip_n)};
imagesc(map');
axis xy
colormap(gca,(magma(256)));
ylabel('Estimated position')
set(gca,'xtick',[])

originalSize1 = get(gca, 'Position');
cb = colorbar('box','off');
set(gca, 'Position', originalSize1);
cb.Title.String = 'Pr';
cb.Position(3) = 0.5*cb.Position(3);

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
%%




