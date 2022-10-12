function Cov_comp = getCompCov_interactions(X,code,P)
    %%% A bit different from the classic audio/video component.
    %%% Would have needed a different architecture to homogeneize.
    
    P.focusOn = 'video'; % force it here
    dataA = getOrgData(X,code,P);

    if P.zscore
        sigma = nanstd(dataA,[],[1 2 3 5]); sigma(sigma == 0) = 1;
        dataA = (dataA - nanmean(dataA,[1 2 3 5]))./sigma;
    end

    s = size(dataA);
    nBins = s(1);
    nVid = s(2);
    nSnd = s(3);
    nClu = s(4);
    nRep = s(5);

    [nfold,folds] = getFolds(s,P.foldtype);

    for ff = 1:nfold
        idxfold = folds{ff};

        spkData_train = dataA;
        spkData_train(idxfold) = nan;

        % get components
        trialNumPerStim = getTrialNumPerStim(code);
        nBins = size(spkData_train,1); % can be different from P.nBins if take average response instead of timecourse

        % get data
        dataI = spkData_train;
        muV = nanmean(dataI,[2 5]);
        dataI = dataI - muV;
        muS = nanmean(dataI,[3 5]);
        dataI = dataI - muS;

        nNeu = size(spkData_train,4);

        % reshape data
        dataI = reshape(dataI, [nBins*nVid*nSnd, nNeu, nRep]);

        % get transition matrix for dimensionality reduction (DSS/SVD/none)
        [M,~,~] = getTransMat(dataI,P);
        nComp = size(M,2);

        % project everything

        UAll = projectData(dataI,M);
        UAll = reshape(UAll,[nBins,nVid,nSnd,size(UAll,2),nRep]);

        % get proj on train set
        % should be the same as going through the whole thing
        projpred = reshape(squeeze(nanmean(UAll(P.bins>0,:,:,:,:),5)),[sum(P.bins>0)*nVid*nSnd,nComp]);

        % going through the same processing pipeline as train set
        spkData_test = dataA;
        spkData_test(~idxfold) = nan;
        spkData_test = spkData_test - muV;
        spkData_test = spkData_test - muS;

        spkData_test = spkData_test(P.bins>0,:,:,:,:);
        spkData_test_res = reshape(spkData_test,[sum(P.bins>0)*nVid*nSnd, nClu,nRep]);
        spkData_test_res = spkData_test_res(:,:,~isnan(spkData_test_res(1,1,:)));
        spkData_test_res = nanmean(spkData_test_res,3);
        proj = projectData(spkData_test_res,M);

        for pc = P.dim.pcL(P.dim.pcL <= nComp)
            % test-retest covariance of the projection on each component
            Cov_comp(pc,ff) = getCov(proj(:,pc),projpred(:,pc));
        end
    end
end