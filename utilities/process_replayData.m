function [data,counts,start,stop] = process_replayData(spkmat, ts, keep, binSize)
% function that 'pulls out' a start and stop time around each
% ripple/population event. This is an entirely arbitrary heuristic, I have
% attempted to match the behavior of this function to some of the most
% common heuristics used by the field. 
%
% INPUTS
%
% ts - timestamp of peak of event (in seconds)


start = round((round(ts(1) * 1000)-50) ./ (spkmat.dt*1000));
stop = round((round(ts(2) * 1000)+50) ./ (spkmat.dt*1000));


for spk = 1:size(spkmat.data,2)
    if stop > length(spkmat.data) | start < 1
        continue
    end
    data(:,spk) = (spkmat.data(start:stop,spk)')';  
    counts(:,spk) = (spkmat.data(start:stop,spk)')';  
end
if ~exist('data','var')
    data = NaN;
    counts = NaN;
    return
end
% cut 0 and 1 spk count bins from the start/end
while sum(counts(1,:)) < 2 && size(counts,1) > 1
    data = data(2:end,:);
    counts = counts(2:end,:);
    start = start + 1;
end
while sum(counts(end,:)) < 2 && size(counts,1) > 1
    data = data(1:end-1,:);
    counts = counts(1:end-1,:);
    stop = stop-1;
end
data = data ./ binSize;

start = (start / 1000) * (spkmat.dt*1000);
stop = (stop / 1000) * (spkmat.dt*1000);
end