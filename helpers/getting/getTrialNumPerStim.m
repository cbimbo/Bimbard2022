function [trialNumPerStim,labelsGroupsSnd,labelsGroupsVid,nGroupsSnd,nGroupsVid] = getTrialNumPerStim(audiovisuoCode)

% get labels
labelsGroupsSnd = unique(audiovisuoCode(2,:));
labelsGroupsVid = unique(audiovisuoCode(1,:));
nGroupsSnd = numel(labelsGroupsSnd);
nGroupsVid = numel(labelsGroupsVid);

% get number of trials pair stim
trialNumPerStim = nan(nGroupsVid,nGroupsSnd);
for v = 1:numel(labelsGroupsVid)
    for s = 1:numel(labelsGroupsSnd)
        trialNumPerStim(v,s) = sum((audiovisuoCode(1,:) == labelsGroupsVid(v)) & (audiovisuoCode(2,:) == labelsGroupsSnd(s)));
    end
end
    
end