function [pc1_ipsi,pc1_contra,lags_xc] = getPCIpsiContra_1ms(db,P,ipsiContraP)

    window = ipsiContraP.window;
    psthBinSize = ipsiContraP.psthBinSize;

    mouseRefUni = unique(P.mouseRef);
    pc1_ipsi = cell(1,numel(mouseRefUni));
    pc1_contra = cell(1,numel(mouseRefUni));
    lags_xc = cell(1,numel(mouseRefUni));
    for k = 1:numel(mouseRefUni)
        ipsiRec = contains(P.mouseRef,mouseRefUni{k}) & [db.side] == 1;
        contraRec = contains(P.mouseRef,mouseRefUni{k}) & [db.side] == -1;
        pc1_ipsi{k} = get1msBinnedProj(db(ipsiRec),P,window,psthBinSize);
        pc1_contra{k} = get1msBinnedProj(db(contraRec),P,window,psthBinSize);

        lags_xc{k} = (-100:100)*psthBinSize;

    end
end

function pc1 = get1msBinnedProj(db,P,window,psthBinSize)
    % Do it on ipsi side
    spikeData_1ms = nan(numel(window(1):psthBinSize:window(2))-1, size(db.spikeData,2), size(db.spikeData,3));
    for c = 1:size(db.spikeData,2)
        clu = db.C.CluID(c);
        [~, ~, ~, ~, ~, ba] = psthAndBA(db.sp.st(db.sp.clu == clu),  ...
            db.eventTimes, window, psthBinSize);
        spikeData_1ms(:,c,:) = ba';
    end

    dataA = getOrgData(db.spikeData,db.audiovisuoCode,P);
    [~, M, ~, mu, sigma] = getComponents(dataA, db.audiovisuoCode, P); % get spike data

    dataA_1ms = getOrgData(spikeData_1ms,db.audiovisuoCode,P);
    muTimeCourse = nanmean(dataA_1ms,[2 3 5]);

    s = size(dataA_1ms);
    nBins = s(1);
    nStimImp = s(2);
    nStimUnImp = s(3);
    nClu = s(4);
    nRep = s(5);

    % reproject everything
    spkData_test = dataA_1ms;
    if P.demeanTimeCourse
        % remove the timecourse in each neuron for each specific unimportant stimulus
        spkData_test = spkData_test - muTimeCourse;
    end
    spkData_test = (spkData_test - mu)./sigma;
    spkData_test_res = reshape(permute(spkData_test,[1 2 4 3 5]),[nBins*nStimImp,nClu,nStimUnImp*nRep]);
    proj = projectData(spkData_test_res,M);
    proj_res = reshape(proj,[nBins,nStimImp,size(proj,2),nStimUnImp,nRep]);
    proj_res = permute(proj_res, [1 2 4 3 5]);
    pc1 = reshape(nanmean(proj_res(:,:,:,:,1),[3 4]),[nBins*nStimImp,1]);
end