function [a, dmatrix] = decodeDataXVal(data, audiovisuoCode, P)

    dataA = getOrgData(data,audiovisuoCode,P); % dataA is bins x imp stim x unimp stim x comp x rep

    s = size(dataA);
    nBins = s(1);
    nStimImp = s(2);
    nStimUnImp = s(3);
    nClu = s(4);
    nRep = s(5);

    [nfold,folds] = getFolds(s,P.foldtype);

    a = nan(numel(P.pc2Keep),nfold);
    dmatrix = nan(nStimImp,nStimImp,numel(P.pc2Keep),nfold);

    for ff = 1:nfold
        idxfold = folds{ff};

        % Get train set
        data_train = dataA;
        data_train(idxfold) = nan;

        % Get test set
        data_test = dataA;
        data_test(~idxfold) = nan;

        % Project both
        [proj_train, proj_test] = getProj(data_train,data_test,P);

        % Reshape it
        nComp = size(proj_train,4);
        UTrain_r = reshape(permute(proj_train,[1 2 4 3 5]),[nBins, nStimImp, nComp, nStimUnImp*nRep]);
        UTrain_r(:,:,:,isnan(squeeze(UTrain_r(1,1,1,:)))) = [];
        UTest_r = reshape(permute(proj_test,[1 2 4 3 5]),[nBins, nStimImp, nComp, nStimUnImp*nRep]);
        UTest_r(:,:,:,isnan(squeeze(UTest_r(1,1,1,:)))) = [];

        if ~iscell(P.pc2Keep)
            disp('P.pc2Keep must be cell')
        end
        for p = 1:numel(P.pc2Keep)
            if ~ischar(P.pc2Keep{p})
                [a(p,ff), dmatrix(:,:,p,ff)] = decodeTemplate(UTrain_r(P.bins>0,:,1:P.pc2Keep{p},:), UTest_r(P.bins>0,:,1:P.pc2Keep{p},:), 0);
            elseif strcmp(P.pc2Keep{p},'all')
                [a(p,ff), dmatrix(:,:,p,ff)] = decodeTemplate(UTrain_r(P.bins>0,:,1:end,:), UTest_r(P.bins>0,:,1:end,:), 0);
            else
                error('no idea what to do with that')
            end
        end
    end

    a = a*100;
end