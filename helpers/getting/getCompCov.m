function Cov_comp = getCompCov(X,code,P)

dataA = getOrgData(X,code,P);

s = size(dataA);
nStimImp = s(2);

[nfold,folds] = getFolds(s,P.foldtype);

Cov_comp = nan(max(P.dim.pcL),nfold);
for ff = 1:nfold
    idxfold = folds{ff};

    % Get train set
    spkData_train = dataA;
    spkData_train(idxfold) = nan;

    % Get test set
    spkData_test = dataA;
    spkData_test(~idxfold) = nan;

    % Project both
    [proj_train, proj_test] = getProj(spkData_train,spkData_test,P);

    % Reshape
    nComp = size(proj_train,4);
    proj_train = reshape(squeeze(nanmean(proj_train(P.bins>0,:,:,:,:),[3 5])),[sum(P.bins>0)*nStimImp,nComp]);
    proj_test = reshape(squeeze(nanmean(proj_test(P.bins>0,:,:,:,:),[3 5])),[sum(P.bins>0)*nStimImp,nComp]);

    for pc = P.dim.pcL(P.dim.pcL <= nComp)
        % Test-retest covariance of the projection on each component
        Cov_comp(pc,ff) = getCov(proj_test(:,pc),proj_train(:,pc));
    end
end