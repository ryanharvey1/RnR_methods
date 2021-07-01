function [replayScores] = replay_Bayesian(df,spikes,ripples)
% Estimates the Bayesian method for replay quantification
%
% USAGE
%
% [] = compareReplayMethods(spikes,ripples,template,include)
%
% INPUTS
%
% spikes  - buzcode cellinfo file. A structure that only requires three fields:
%              -spikes.times, an Nx1 cell array of timestamps in seconds for each neuron
%              -spikes.UID, Nx1 vector of unique identifiers for each neuron in session
%              -spikes.spindices, a sorted vector of [spiketime UID], for all spikes.
%                               useful input to some functions and plotting rasters
% ripples - buzcode ripples events file. A structure that only requires two fields:
%             -ripples.timestamps, an Nx2 matrix of start/stop times
%             -ripples.peaks, an Nx1 vector of peak timestamps for each event
% template -NxD matrix of N cells and D positions, average firing ratemaps
%           for each cell
% include - indices (1:N) of cells (i.e.place cells) to keep
%
% OUTPUTS
%
% bayesRadon - integral under the line of best fit, using the Radon
%              transform (Davidson 2009)
% bayesLinearWeighted - linear weighted correlation of the posterior
%              probability matrix (Grosmark 2016)
%
% HELP
%
% See bz_GetSpikes from the buzcode repo for help with the spikes/ripples
%   data structures
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

nBinsThresh = 5; % the minimum number of bins, per event, to analyze event (binSize * nBinsThresh = total event duration)
binSize = .01;
overlap = 1;
% spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap); % converts spike timestamps into matrix format
% template = template(include,:);

% if ~isempty(include)
%     keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
% else
%     keep = find(sum(template')>0); % just remove cells w/ zeros spikes in template
% end

for event = 1:size(ripples.timestamps,1)
    % loop through each running direction
    for d = 1:2
        % find place cells
        [include,~,template] = get_place_cells(df,d);
        include = 1:size(template,1);
        % skip if no place cells
        if isempty(include)
            [bayesLinearWeighted(event,d),...
                slope_hpc(event,d),...
                bayesRadon(event,d),...
                Prs{event,d},...
                inactive_bins{event,d},...
                z_cellID_shuf(event,d),...
                pvalue_cellID_shuf(event,d),...
                z_circular_shuf(event,d),...
                pvalue_circular_shuf(event,d),...
                nCells(event,d),...
                nSpks(event,d)] = set_nan();
            start(event,d) = nan;
            stop(event,d) = nan;
            continue
        end
        
        % converts spike timestamps into matrix format
        spkmat = bz_SpktToSpkmat(spikes.times(include),'overlap',overlap,'binSize',binSize * overlap);
        
        % extract and preprocess ripple
        
        % gets start/stop timestamp of ripple event
        ts = ripples.timestamps(event,:);
        
        % processing func to get out the 'event' using common FR heuristics
        [data,counts,start(event,d),stop(event,d)] = process_replayData(spkmat, ts, include, binSize);
        
        % find time bins with no activity
        inactive_bins{event,d} = sum(data,2) == 0;
        
        % decode ripple
        [bayesLinearWeighted(event,d),...
            slope_hpc(event,d),...
            bayesRadon(event,d),...
            Prs{event,d},...
            z_cellID_shuf(event,d),...
            pvalue_cellID_shuf(event,d),...
            z_circular_shuf(event,d),...
            pvalue_circular_shuf(event,d),...
            z_column_cycle(event,d),...
            pvalue_column_cycle(event,d),...
            nCells(event,d),...
            nSpks(event,d)] = decode_ripple(data,counts,template(include,:),nBinsThresh,spkmat);
    end
end

% find direction that maximizes z_circular_shuf
[~,b]=min(pvalue_cellID_shuf, [], 2);

% create data struct to return
for i = 1:length(b)
    replayScores.bayesLinearWeighted(i,1) = bayesLinearWeighted(i,b(i));
    replayScores.bayesRadon(i,1) = bayesRadon(i,b(i));
    replayScores.slope_hpc(i,1) = slope_hpc(i,b(i));
    
    replayScores.pvalue_cellID_shuf(i,1) = pvalue_cellID_shuf(i,b(i));
    replayScores.z_cellID_shuf(i,1) = z_cellID_shuf(i,b(i));
    replayScores.pvalue_circular_shuf(i,1) = pvalue_circular_shuf(i,b(i));
    replayScores.z_circular_shuf(i,1) = z_circular_shuf(i,b(i));
    
    replayScores.pvalue_column_cycle(i,1) = pvalue_column_cycle(i,b(i));
    replayScores.z_column_cycle(i,1) = z_column_cycle(i,b(i));
    
    replayScores.nCells(i,1) = nCells(i,b(i));
    replayScores.nSpks(i,1) = nSpks(i,b(i));
    
    replayScores.Pr{i,1} = Prs{i,b(i)};
    replayScores.inactive_bins{i,1} = inactive_bins{i,b(i)};

    replayScores.start(i,1) = start(i,b(i));
    replayScores.stop(i,1) = stop(i,b(i));
    
    replayScores.direction_used(i,1) = b(i);
end
end

function z = get_rZ(obs,null)
z = (abs(obs) - abs(mean(null))) / abs(std(null));
end

function p = get_pvalue(obs,null)
p = (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1);
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

function  [bayesLinearWeighted,slope_hpc,bayesRadon,Prs,...
    z_cellID_shuf,pvalue_cellID_shuf,z_circular_shuf,pvalue_circular_shuf,z_column_cycle,pvalue_column_cycle,...
    nCells,nSpks] = decode_ripple(data,counts,template,nBinsThresh,spkmat)

if size(data,1) >= nBinsThresh && sum(any(data)) > 5
    [Pr, prMax] = placeBayes(data, template, spkmat.dt); % generate posterior probability matrix using template and event FRs
    Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
    Prs = Pr;
    [bayesLinearWeighted,outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
    [slope_hpc,bayesRadon,intercept] = Pr2Radon(Pr'); % Radon transform method of quantification (Davidson/Frank 2009)
    
    %% let's add some shuffling
    for i = 1:1000
        % cell ID
        [Pr, ~] = placeBayes(data, bz_shuffleCellID(template), spkmat.dt); % generate posterior probability matrix using template and event FRs
        Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
        [bayesLinearWeighted_cellID_shuf(i),outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
        
        % circular
        [Pr, ~] = placeBayes(data, bz_shuffleCircular(template), spkmat.dt); % generate posterior probability matrix using template and event FRs
        Pr(isnan(Pr)) = 0; % bad form... but some events have 0 spks in a particular time bin, doing this rather than filtering those events out
        [bayesLinearWeighted_circular_shuf(i),outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1)); % linear weight correlation method of quantification (Grosmark/Buzsaki 2016)
        
        % column cycle shuff
        r = floor(1 + (size(Prs,2)-1).*rand(size(Prs,1),1));
        for s = 1:size(Prs,1)
            temp_map(s,:) = circshift(Prs(s,:), r(s));
        end
        [bayesLinearWeighted_column_shuf(i),~] = makeBayesWeightedCorr1(temp_map,ones(size(temp_map,1),1));
    end

    
    
    % get z and pvalues based on shuffled null distribution
    z_cellID_shuf = get_rZ(bayesLinearWeighted,bayesLinearWeighted_cellID_shuf);
    pvalue_cellID_shuf = get_pvalue(bayesLinearWeighted,bayesLinearWeighted_cellID_shuf);
    
    z_circular_shuf = get_rZ(bayesLinearWeighted,bayesLinearWeighted_circular_shuf);
    pvalue_circular_shuf = get_pvalue(bayesLinearWeighted,bayesLinearWeighted_circular_shuf);
    
    z_column_cycle = get_rZ(bayesLinearWeighted,bayesLinearWeighted_column_shuf);
    pvalue_column_cycle = get_pvalue(bayesLinearWeighted,bayesLinearWeighted_column_shuf);
    
    nCells = sum(sum(counts)>0);
    nSpks = sum(sum(counts));
    
else
    [bayesLinearWeighted,slope_hpc,bayesRadon,Prs,~,...
        z_cellID_shuf,pvalue_cellID_shuf,...
        z_circular_shuf,pvalue_circular_shuf,...
        z_column_cycle,pvalue_column_cycle,...
        nCells,nSpks] = set_nan();
end
end

function  [bayesLinearWeighted,...
    slope_hpc,...
    bayesRadon,...
    Prs,...
    inactive_bins,...
    z_cellID_shuf,...
    pvalue_cellID_shuf,...
    z_circular_shuf,...
    pvalue_circular_shuf,...
    z_column_cycle,...
    pvalue_column_cycle,...
    nCells,...
    nSpks] = set_nan()

bayesLinearWeighted = nan;
slope_hpc = nan;
bayesRadon = nan;

Prs = nan;
inactive_bins = nan;

z_cellID_shuf = nan;
pvalue_cellID_shuf = nan;
z_circular_shuf = nan;
pvalue_circular_shuf = nan;
z_column_cycle= nan;
pvalue_column_cycle= nan;
nCells = nan;
nSpks = nan;
end