data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
df = readtable('F:\Projects\PAE_PlaceCell\replay\processed\replay_df.csv');
Pr = load('F:\Projects\PAE_PlaceCell\replay\processed\replay_Pr.mat','Pr');
load('F:\Projects\PAE_PlaceCell\replay\processed\replay_inactive_bins.mat','inactive_bins');

df.ripple_number = df.Var1;
df_temp = df;


% ripple to show
ripples_to_use_idx = [23628,6547];
for i = 1:length(ripples_to_use_idx)
    idx(i) = find(df_temp.ripple_number == ripples_to_use_idx(i));
end

% loop through potential replay events
for rip_n = idx
    data = load([data_path,df_temp.session{rip_n},'.mat']);
   
    plot_replay(data,df_temp,df,Pr,rip_n,df_temp.start(rip_n)-10,df_temp.stop(rip_n)+10,inactive_bins)
end

% make plot
function plot_replay(data,df_temp,df,Pr_replay,rip_n,start,stop,inactive_bins)

pr_start = df_temp.start(rip_n);
pr_stop = df_temp.stop(rip_n);

% set figure defaults
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','Serif','defaultTextFontName','Serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('thesis', 0.5, [1.5,1]);
set(fig,'Position',[1 1 fig_width_in fig_height_in])

subplot(2,1,1)
plot_ripple_raster(data,start,stop,pr_start,pr_stop,df_temp,rip_n)

subplot(2,1,2)

map = Pr_replay.Pr{df.ripple_number == df_temp.ripple_number(rip_n)};

imagesc(map');
hold on
axis xy

colormap(gca,(magma(256)));
ylabel('Estimated position')
set(gca,'xtick',[])

originalSize1 = get(gca, 'Position');
cb = colorbar('box','off');
set(gca, 'Position', originalSize1);
cb.Title.String = 'Pr';
cb.Position(3) = 0.5*cb.Position(3);


x = 1:size(map,1);
[w,y] = max(map,[],2);

y(inactive_bins{df.ripple_number == df_temp.ripple_number(rip_n)}) = [];
x(inactive_bins{df.ripple_number == df_temp.ripple_number(rip_n)}) = [];

plot(x,y,'w','linewidth',2)
end

% below are functions to make each panel in the plot
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

%%




