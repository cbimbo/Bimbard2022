function [dataMPC,M,dataMPC_odd,dataMPC_even,dataA] = getMainFigData(db,P,mainFigP)
    
    %% get PCs for all mice
    P.focusOn = mainFigP.focusOn;
    
    dataMPC = nan(numel(P.bins),P.nGroupsSnd,numel(mainFigP.pc2Keep),P.nMice);
    dataMPC_odd = nan(numel(P.bins),P.nGroupsSnd,numel(mainFigP.pc2Keep),P.nMice);
    dataMPC_even = nan(numel(P.bins),P.nGroupsSnd,numel(mainFigP.pc2Keep),P.nMice);
    M = cell(1,P.nMice);
    for k = 1:P.nMice
        dataA = getOrgData(db(k).spikeData,db(k).audiovisuoCode,P);

        [~,folds] = getFolds(size(dataA),P.foldtype);
        
        ff = 2; % train is odd, test is even
        idxfold = folds{ff};
        
        % Get train set
        spkData_train = dataA;
        spkData_train(idxfold) = nan;

        % Get test set--all data for graphical purpose
        spkData_test = dataA;

        % Project everything
        [~, proj_test, M{k}] = getProj(spkData_train,spkData_test,P);
        
        dataAPC = proj_test(:,:,:,mainFigP.pc2Keep,:);
        
        % Little processing
        for pc = 1:numel(mainFigP.pc2Keep)
            % Normalize it
            sigma = nanstd(mat2vec(dataAPC(:,:,:,pc,:)));
            if sigma>0
                dataAPC(:,:,:,pc,:) = dataAPC(:,:,:,pc,:)/sigma;
            end
            
            % Resign it -- arbitrary
            si = sign(skewness(mat2vec(nanmean(dataAPC(:,:,:,pc,:),[3 5])))); % reorient it
            dataAPC(:,:,:,pc,:) = dataAPC(:,:,:,pc,:)*si;
            M{k}(:,pc) = M{k}(:,pc)*si;
        end

        % Average it
        dataMPC(:,:,:,k) = squeeze(nanmean(dataAPC,[3 5]));

        % get odd and even
        dataMPC_odd(:,:,:,k) = squeeze(nanmean(dataAPC(:,:,:,:,1:2:end),[3 5]));
        dataMPC_even(:,:,:,k) = squeeze(nanmean(dataAPC(:,:,:,:,2:2:end),[3 5]));
    end
    
    
    %% get raw data for main example mouse
    %%% ugly hack 
    switch mainFigP.focusOn
        case 'sound'
            P.focusOn = 'video';
        case 'video'
            P.focusOn = 'sound';
    end
    [dataA, ~] = getOrgData(db(mainFigP.k2plt(1)).spikeData,db(mainFigP.k2plt(1)).audiovisuoCode,P);
end