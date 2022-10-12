function [xc,lags_xc,w,lags_w] = getLags(db,P,lagP)

    %% Get the binned data
    [~,neurProj,motDat,baseIdx] = getDataBinnedToVideoFrames(db,P,lagP);

    for k = 1:P.nMice
        % Take only baseline
        neurProj{k}(~baseIdx{k},:) = nan;

        % Subselect only 1st component
        neurProjPC1{k} = neurProj{k}(:,1);
        motDatPC1{k} = motDat{k}(:,1);
    end

    %% Get cross-correlogram
    for k = 1:P.nMice
        lags_xc{k} = (-10:10)*mean(diff(db(k).videoTimestampsEphysTime_face));
    end
    xc = getCrossCorrelogram(neurProjPC1,motDatPC1,lags_xc);

    %% Run model

    if lagP.runModel

        pc2predL = lagP.pc2predL; % number of neural PCs to predict
        motpc2keepfix = lagP.motpc2keepfix;
        model = lagP.model;

        lags_w = cell(1,P.nMice);
        w = cell(1,P.nMice);
        for k = 1:P.nMice
            lagbins = -ceil(0.2525/mean(diff(db(k).videoTimestampsEphysTime_face))):ceil(0.3787/mean(diff(db(k).videoTimestampsEphysTime_face)));
            lags_w{k} = lagbins*mean(diff(db(k).videoTimestampsEphysTime_face));

            y = nanzscore(neurProj{k}(:,pc2predL));
            c = corr(motDat{k}(~isnan(y(:,1)),1),y(~isnan(y(:,1))));
            X = getXpred(nanzscore(motDat{k}(:,1:motpc2keepfix)),lagbins);
            X = X*sign(c);
            w{k} = makePred_paper(y,X,[],model);
            w{k} = reshape(w{k}(1:end-1,:),[numel(lagbins),motpc2keepfix,numel(pc2predL)]); % last one is baseline?
        end
    else
        w = [];
        lags_w = [];
    end
end