function U = projectData(d,M)

    %%% Will project each repetition of d on PC space through M
    %%% d must be time x clusters x repetitions
    %%% U is time x comps x repetitions
    
    % get a few parameters
    nTime = size(d,1);
    nClu = size(d,2);
    nTrials = size(d,3);
    nComp = size(M,2);
    
    d = reshape(permute(d, [1 3 2]), [nTime*nTrials,nClu]);
    U = permute(reshape(d*M,[nTime,nTrials,nComp]), [1 3 2]);
    
end