data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
data = load([data_path,'LEM3206_S20190724151552.mat'],'Spikes','ratemap','events','frames','linear_track');
%%
% r_squared = [];
% for binSize = .01:.01:.6
binSize = .2;
overlap = 4;
direction = 2;

[data,frames] = restrict_track(data);
data = get_spike_times(data,direction);
spikes = format_spikes(data);


[include,~,template] = get_place_cells(data,direction);

template = template(include,:);

spkmat = bz_SpktToSpkmat(spikes.times(include),'overlap',overlap,'binSize',binSize * overlap);
% spkmat = bz_SpktToSpkmat(spikes.times(include),'binSize',binSize);

idx = sum(spkmat.data,2);
spkmat.data(idx == 0,:) = [];
spkmat.timestamps(idx == 0,:) = [];

[Pr, x_decoded] = placeBayes(spkmat.data, template, spkmat.dt); % generate posterior probability matrix using template and event FRs
frames(:,2) = rescale(frames(:,2),1,40);
x_obs = interp1(frames(:,1),frames(:,2),spkmat.timestamps);

figure;
subplot(2,2,1)
scatter(x_obs,x_decoded,'black','Filled')
subplot(2,2,[3,4])
plot(x_obs);hold on; plot(x_decoded);
mdl = fitlm(x_obs,x_decoded);

decoding_error = abs(x_obs-x_decoded);

time_window = [1,1000];
figure;imagesc(Pr(time_window(1):time_window(2),:)');
colormap magma
hold on;
plot(x_obs(time_window(1):time_window(2)),'w')

% r_squared = [r_squared;mdl.Rsquared.Ordinary];
% end

figure
X = cellfun(@transpose,data.Spikes(include),'un',0);
X = X(~cellfun('isempty',X));

% maps = template(include,:);
[~,I] = max(template,[],2);
[~,I] = sort(I,'descend');
% maps = maps(I,:);
%     
[~,I] = sort(cellfun(@length,X));
X = X(I);
clear x y
for i=1:length(X)
    [x{i},y{i}]=plotSpikeRaster(X(i),'PlotType','vertline');
end
for i=1:length(x)
    plot(x{i},y{i}+i,'Color',[rand(1),rand(1),rand(1)]);hold on
end
axis tight
grid on
xlabel('Time(sec)')
ylabel('Cells')
title('Raster')

function [data,frames] = restrict_track(data)
for i = 1:length(data.Spikes)
    data.Spikes{i} = data.Spikes{i}(data.Spikes{i} >= data.events(1,1) &...
        data.Spikes{i} <= data.events(2,1));
end
frames = data.frames(data.frames(:,1) >= data.events(1,1) &...
        data.frames(:,1) <= data.events(2,1),:);
end

function data = get_spike_times(data,direction)
if direction == 1
    direction = 'right';
else
    direction = 'left';
end
for i = 1:length(data.linear_track{1, 1}.(direction))
    idx = data.linear_track{1, 1}.(direction){1,i}.dataspks(:,6);
    data.Spikes{i} = data.linear_track{1, 1}.(direction){1, i}.dataspks(logical(idx),1);
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

function [include,peak_loc,template] = get_place_cells(data,d)
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
%%