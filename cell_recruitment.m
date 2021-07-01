% cell_recruitment
df_swr = readtable('F:\Projects\PAE_PlaceCell\swr_data\post_processed\swr_df.csv');

df_cell_class = readtable('F:\Projects\PAE_PlaceCell\cell_recruitment\processed\pyr_int_df.csv');
data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';

i = 1;
sessions = unique(df_swr.session);
for i = 1:length(sessions)
    data = load([data_path,df_swr.session{i},'.mat']);
    
    % ripple loop
    idx = contains(df_swr.session,sessions{i});
    
    start = df_swr.peak_time(idx) - .5;
    stop = df_swr.peak_time(idx) + .5;
    
    
    
    temp_spikes = data.Spikes;
    
    for r = 1:length(start)
        for s = 1:length(data.Spikes)
            rip_ts{s,r} = data.Spikes{s}(data.Spikes{s} >= start(r) & data.Spikes{s} <= stop(r)) - df_swr.peak_time(r);
        end
    end
    
    clear x y
    for i=1:size(rip_ts,1)
        if isempty(rip_ts{i})
            x{i} = [];
            y{i} = [];
            continue
        end
        [x{i},y{i}]=plotSpikeRaster(temp_spikes(i),'PlotType','vertline');
    end
end

figure
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
% xlim([pr_start,pr_stop])
for i=1:length(x)
    plot(x{i},y{i}+(i-1),'Color',colors(i,:),'linewidth',2);hold on
end
box off
xlabel('Time (s)')
ylabel('Unit #')
set(gca,'ytick',[])
xlim([pr_start,pr_stop])