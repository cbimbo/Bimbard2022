function [UAll, M, muTimeCourse, mu, sigma, S, U] = getComponents(data, code, P)
    %%% Will compute the components of data
    %%% Supposes that data is timepoints x imp stim x unimp stim x clusters x repeats

    nBins = size(data,1); % can be different from P.nBins if take average response instead of timecourse
    P.nBins = nBins; % modify this internally
    if nBins == 1
        P.bins = 1;
    end
    
    % Get data
    if numel(size(data)) == 3
        trialNumPerStim = getTrialNumPerStim(code);
        dataA = getOrgData(data,code,P);
    else
        dataA = data;
    end
    nStimImp = size(dataA,2);
    nStimUnImp = size(dataA,3);
    nComps = size(dataA,4);
    nRep = size(dataA,5);
    
    if P.demeanTimeCourse
        % remove the timecourse in each neuron for each specific unimportant stimulus
        muTimeCourse = nanmean(dataA,[2,5]);
        dataA = dataA - muTimeCourse;
    else 
        muTimeCourse = [];
    end
    
    % zscore data
    if P.zscore
        mu = nanmean(dataA,[1 2 3 5]);
        sigma = nanstd(dataA,[],[1 2 3 5]); sigma(isnan(sigma)) = 1; sigma(sigma == 0) = 1;
    else
        mu = 0;
        sigma = 1;
    end
    dataA = (dataA - mu)./sigma;
    
    % reshape data
    dataA = reshape(permute(dataA,[1 2 4 3 5]), [nBins*nStimImp, nComps, nRep*nStimUnImp]);
    
    % get transition matrix for dimensionality reduction (DSS/SVD/none)
    [M,S,U] = getTransMat(dataA,P);
    
    % project everything
    if numel(size(data)) == 3
        % won't have the unimp stim subtraction
        data = (data - squeeze(mu)')./squeeze(sigma)';
        UAll = projectData(data,M);
    else
        UAll = projectData(dataA,M);
        UAll = reshape(UAll,[nBins,nStimImp,size(UAll,2),nStimUnImp,nRep]);
        UAll = permute(UAll,[1 2 4 3 5]);
    end
end