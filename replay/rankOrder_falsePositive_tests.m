% this script attempts to calculate false positive rates for the rank order 
% correlation method given a number typical numbers of neurons and events
% recorded
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


seqRange = 3:10; % the typical number of unique neurons that spike in a given ripple or candidate replay event
numEvents = 100; % typical number of candidate events to examine, within a single session
seqLeng = nan(length(seqRange),1);
seqLengMax = nan(length(seqRange),1);
numIterations = 100;

for seqLen = seqRange
    co = nan(numIterations,numEvents);
    pval = nan(numIterations,numEvents);
    
    for iter = 1:numIterations
        seq1=1:seqLen;
        for rip = 1:numEvents
            seq2 = seq1(randperm(seqLen));
            [co(iter,rip),pval(iter,rip)] = corr(seq1',seq2');
        end
    end
    seqLeng(seqLen) = mean(sum(pval<.05,2)); % average FP rate
    seqLengMax(seqLen) = max(sum(pval<.05,2)); % the worst case FP rate
end

subplot(2,2,1)
histogram(sum(pval<.05,2))

subplot(2,2,2)
plot(seqLeng)
hold on
plot(seqLengMax)



