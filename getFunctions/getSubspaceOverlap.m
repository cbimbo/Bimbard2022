function [Cov_comp_vid,Cov_comp_snd,percExpl_vid,percExpl_snd] = getSubspaceOverlap(db,P,subOvP)

    %% Get the binned data
    [neurAll,~,motDat,baseIdx] = getDataBinnedToVideoFrames(db,P,subOvP);

    %% Get subspace overlap

    motionCompNum = subOvP.motionCompNum;
    lagbins = subOvP.lagbins;
    shuff = subOvP.shuff;
    meth = subOvP.meth;
    nModelComp2keep = subOvP.nModelComp2keep;

    Nc = nan(1,P.nMice);
    for k = 1:P.nMice
        Nc(k) = size(db(k).spikeData,2);
    end

    fprintf('*** Computing subspace overlap... ***\n')
    Cov_comp_vid.vid = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_vid.snd = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_vid.mov = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_vid.rand = nan(nModelComp2keep,2,P.nMice,shuff);
    Cov_comp_snd.vid = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_snd.snd = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_snd.mov = nan(nModelComp2keep,2,P.nMice);
    Cov_comp_snd.rand = nan(nModelComp2keep,2,P.nMice,shuff);
    for k = 1:P.nMice
        fprintf('Computing for mouse %s...', P.mouseRef{k})

        %% Get the movement-related subspace
        % Get spike and motion data
        spkD = neurAll{k};
        motD = getXpred(motDat{k}(:,1:motionCompNum),lagbins,1);

        % Get baseline period
        spkD = spkD(baseIdx{k},:);
        motD = motD(baseIdx{k},:);

        % Get movement-related subspace
        Mmov = movSubspace(spkD,motD,nModelComp2keep,meth);

        %% Get random subspaces
        Mrand = nan([size(Mmov),shuff]);
        for ss = 1:shuff
            Mrandtmp = randn([size(Mmov,1) size(Mmov,1)]);
            [~,~,Mrandtmp] = svd(Mrandtmp);
            Mrandtmp = Mrandtmp(:,1:size(Mmov,2));
            Mrand(:,:,ss) = Mrandtmp;
        end

        %% Get sound- and video-related responses
        Pvis = P;
        Pvis.focusOn = 'video';
        dataA_vid = getOrgData(db(k).spikeData,db(k).audiovisuoCode,Pvis);

        Psnd = P;
        Psnd.focusOn = 'sound';
        dataA_snd = getOrgData(db(k).spikeData,db(k).audiovisuoCode,Psnd);

        %% Loop through the folds

        [nfold,folds] = getFolds(size(dataA_vid),P.foldtype);
        for ff = 1:nfold
            %% Get train  and test sets for video- and sound-related activity
            idxfold = folds{ff};

            % Video-related
            spkData_train_vid = dataA_vid;
            spkData_train_vid(idxfold) = nan;

            spkData_test_vid = dataA_vid;
            spkData_test_vid(~idxfold) = nan;

            % Sound-related
            spkData_train_snd = dataA_snd;
            spkData_train_snd(idxfold) = nan;

            spkData_test_snd = dataA_snd;
            spkData_test_snd(~idxfold) = nan;

            %% Get video- and sound-related subspaces

            [~, Mvid, muTimeCourseVid, muVid, sigmaVid] = getComponents(spkData_train_vid, db(k).audiovisuoCode, Pvis); % get spike data
            Mvid = Mvid(:,1:nModelComp2keep); % Truncate it
            [~, Msnd, muTimeCourseSnd, muSnd, sigmaSnd] = getComponents(spkData_train_snd, db(k).audiovisuoCode, Psnd); % get spike data
            Msnd = Msnd(:,1:nModelComp2keep); % Truncate it

            %% Prepare video- and sound-related activity

            spkData_train_vid_res = prepareSpkData(spkData_train_vid,muTimeCourseVid,muVid,sigmaVid);
            spkData_test_vid_res = prepareSpkData(spkData_test_vid,muTimeCourseVid,muVid,sigmaVid);

            spkData_train_snd_res = prepareSpkData(spkData_train_snd,muTimeCourseSnd,muSnd,sigmaSnd);
            spkData_test_snd_res = prepareSpkData(spkData_test_snd,muTimeCourseSnd,muSnd,sigmaSnd);

            %% Get train-test covariance for both, on all subspaces
         
            % Video
            [Cov_comp_vid.vid(:,ff,k), Cov_comp_vid.snd(:,ff,k), Cov_comp_vid.mov(:,ff,k), Cov_comp_vid.rand(:,ff,k,:)] = ...
                getTrainTestCovariance(spkData_train_vid_res,spkData_test_vid_res,Mvid,Msnd,Mmov,Mrand);

            % Sound
            [Cov_comp_snd.vid(:,ff,k), Cov_comp_snd.snd(:,ff,k), Cov_comp_snd.mov(:,ff,k), Cov_comp_snd.rand(:,ff,k,:)] = ...
                getTrainTestCovariance(spkData_train_snd_res,spkData_test_snd_res,Mvid,Msnd,Mmov,Mrand);

        end
        fprintf('Done.\n')
    end

    %% Get percentage variance explained by each subspace
    % For video-related activity
    percExpl_vid = getPercExplVariance(Cov_comp_vid,'vid');

    % For sound-related activity
    percExpl_snd = getPercExplVariance(Cov_comp_snd,'snd');

    fprintf('*** Computing subspace overlap done. ***\n')
end


%% Helper functions

function Mmov = movSubspace(spkD,motD,nModelComp2keep,meth)
    switch meth
        case 'CCA'
            % Get canon corr
            A = canoncorr(spkD,motD);

            % Orthogonalize somehow
            [Q,~] = qr(A);
            A = Q; % should be orthonormal

            % Normalize
            A = A./sqrt(sum(A.^2,1));

            % Save in Mmov
            Mmov = A;

        case 'PLS'
            % Perform partial least square
            XL = plsregress(spkD,motD,nModelComp2keep);

            % Normalize
            XL = XL./sqrt(sum(XL.^2,1)); % not sure why I need this?

            % Save in Mmov
            Mmov = XL;

            % Orthogonalize somehow
            [Mmov,~] = qr(Mmov);

        case 'RRR'
            Xtrain = motD;
            ytrain = spkD;
            norm = std(ytrain);
            ytrain = ytrain./norm; % each neuron will have the same chance

            % Perform OLS prediction
            bols = pinv(Xtrain'*Xtrain)*Xtrain'*ytrain;

            % Get predicted y
            ytrain_pred = Xtrain*bols;

            % Get the svd decomposition
            [~,~,v] = svd(ytrain_pred,'econ');

            % Save it in Mmov
            Mmov = v(:,1:nModelComp2keep);
    end
end

function spkData_res = prepareSpkData(spkData,muTimeCourse,mu,sigma)
    s = size(spkData);
    nBins = s(1);
    nStimImp = s(2);
    nStimUnImp = s(3);
    nClu = s(4);
    nRep = s(5);

    if ~isempty(muTimeCourse)
        % remove the timecourse in each neuron for each specific unimportant stimulus
        spkData = spkData - muTimeCourse;
    end
    spkData = (spkData - mu)./sigma; 
    spkData_res = reshape(permute(spkData,[1 2 4 3 5]),[nBins*nStimImp, nClu,  nStimUnImp*nRep]);
    spkData_res = nanmean(spkData_res,3);
end

function [proj_vidPC, proj_sndPC, proj_movPC, proj_randPC] = projectOnAllSubspaces(spkData_res,Mvid,Msnd,Mmov,Mrand)
    proj_vidPC = projectData(spkData_res,Mvid);
    proj_sndPC = projectData(spkData_res,Msnd);
    proj_movPC = projectData(spkData_res,Mmov);
    proj_randPC = nan([size(proj_vidPC),size(Mrand,3)]);
    for ss = 1:size(Mrand,3)
        proj_randPC(:,:,ss) = projectData(spkData_res,Mrand(:,:,ss));
    end
end

function [Cov_comp_vid,Cov_comp_snd,Cov_comp_mov,Cov_comp_rand] = getTrainTestCovariance(spkData_train_res,spkData_test_res,Mvid,Msnd,Mmov,Mrand)

    % Project train set
    [proj_train_vidPC, proj_train_sndPC, proj_train_movPC, proj_train_randPC] = ...
        projectOnAllSubspaces(spkData_train_res,Mvid,Msnd,Mmov,Mrand);

    % Project test set
    [proj_test_vidPC, proj_test_sndPC, proj_test_movPC, proj_test_randPC] = ...
        projectOnAllSubspaces(spkData_test_res,Mvid,Msnd,Mmov,Mrand);

    % Compute covariance
    nComp = size(Mvid,2);
    shuff = size(Mrand,3);
    Cov_comp_vid = nan(nComp,1);
    Cov_comp_snd = nan(nComp,1);
    Cov_comp_mov = nan(nComp,1);
    Cov_comp_rand = nan(nComp,shuff);
    for pc = 1:nComp
        % test-retest covariance of the projection on each component
        Cov_comp_vid(pc) = getCov(proj_test_vidPC(:,pc),proj_train_vidPC(:,pc));
        Cov_comp_snd(pc) = getCov(proj_test_sndPC(:,pc),proj_train_sndPC(:,pc));
        Cov_comp_mov(pc) = getCov(proj_test_movPC(:,pc),proj_train_movPC(:,pc));
        for ss = 1:shuff
            Cov_comp_rand(pc,ss) = getCov(proj_test_randPC(:,pc,ss),proj_train_randPC(:,pc,ss));
        end
    end
end

function percExpl = getPercExplVariance(Cov_comp,ref)
    nMice = size(Cov_comp.vid,3);

    percExpl.vid = nan(nMice,1);
    percExpl.snd = nan(nMice,1);
    percExpl.mov = nan(nMice,1);
    percExpl.rand = nan(nMice,size(Cov_comp.rand,4));
    for k = 1:nMice
        refVarExpl = sum(Cov_comp.(ref)(:,:,k));

        percExpl.vid(k) = nanmean(sum(Cov_comp.vid(:,:,k))./refVarExpl)*100;
        percExpl.snd(k) = nanmean(sum(Cov_comp.snd(:,:,k))./refVarExpl)*100;
        percExpl.mov(k) = nanmean(sum(Cov_comp.mov(:,:,k))./refVarExpl)*100;
        for ss = 1:size(Cov_comp.rand,4)
            percExpl.rand(k,ss) = nanmean(sum(Cov_comp.rand(:,:,k,ss))./refVarExpl)*100;
        end
    end
end