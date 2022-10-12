function [dataA, dataM] = getOrgData(dat,code,P)
    
    % get organized version of the data
    % dataA contains all / dataM the mean over repeats
    % dataA is bins x important stim x unimportant stim x clusters x rep
    %%% supposes that the number of trials is the same for a and v
    
    if ~isempty(dat)
        trialNumPerStim = getTrialNumPerStim(code);
        
        s = size(dat);
        dataA = nan(s(1),P.nGroupsVid,P.nGroupsSnd,s(2),max(trialNumPerStim(:)));
        for gv = 1:P.nGroupsVid
            for ga = 1:P.nGroupsSnd
                tmp = dat(:,:,(code(1,:)==P.labelsGroupsVid(gv)) & (code(2,:)==P.labelsGroupsSnd(ga)));
                dataA(:,gv,ga,:,1:size(tmp,3)) = tmp;
            end
        end
        
        if isfield(P,'focusOn') && strcmp(P.focusOn,'sound')
            dataA = permute(dataA,[1 3 2 4 5]);
        end
        
        dataM = nanmean(dataA,5);
    else
        dataA = [];
        dataM = [];
    end
end