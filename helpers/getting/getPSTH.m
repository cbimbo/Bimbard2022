function [binnedArray, binCenters] = getPSTH(spikeTimes, eventTimes, window, binSize)
    %%% Adapted from psthAndBA from the 'spikes' toolbox (https://github.com/cortex-lab/spikes). 

    spikeTimes = spikeTimes(:);
    eventTimes = sort(eventTimes(:));

    spikeTimes = spikeTimes(spikeTimes>min(eventTimes+window(1)) & spikeTimes<max(eventTimes+window(2)));

    binBorders = window(1):binSize:window(2);
    numBins = length(binBorders)-1;

    if isempty(eventTimes)
        binnedArray = [];
        binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
        return;
    end

    binnedArray = zeros(length(eventTimes), numBins);

    if isempty(spikeTimes)
        binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
        return;
    end

    for r = 1:length(eventTimes)
        [n,binCenters] = histdiff(spikeTimes, eventTimes(r), binBorders);
        binnedArray(r,:) = n;
    end
end