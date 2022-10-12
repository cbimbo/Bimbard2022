function [proj_train, proj_test, M, muTimeCourse, mu, sigma, spkData_test_res] = getProj(spkData_train,spkData_test,P)
    %%% Will compute the components on the train set, and project the test
    %%% set.

    %% Get components from the train set
    [proj_train, M, muTimeCourse, mu, sigma] = getComponents(spkData_train, [], P); % get spike data

    %% Process and then project the test set
    s = size(spkData_test);
    nBins = s(1);
    nStimImp = s(2);
    nStimUnImp = s(3);
    nClu = s(4);
    nRep = s(5);

    % Go through same processing as train data
    if ~isempty(muTimeCourse)
        spkData_test = spkData_test - muTimeCourse;
    end
    spkData_test = (spkData_test - mu)./sigma; 
    spkData_test_res = reshape(permute(spkData_test,[1 2 4 3 5]),[nBins*nStimImp,nClu,nStimUnImp*nRep]);
    
    % Project test data
    proj = projectData(spkData_test_res,M);
    proj_test = reshape(proj,[nBins,nStimImp,size(proj,2),nStimUnImp,nRep]);
    proj_test = permute(proj_test, [1 2 4 3 5]);
end