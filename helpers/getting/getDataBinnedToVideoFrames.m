function [neurAll,neurProj,motDat,baseIdx] = getDataBinnedToVideoFrames(db,P,binP)

    if isfield(binP,'focusOn')
        P.focusOn = binP.focusOn;
    else
        P.focusOn = 'video';
    end

    %% clean time when flashes for alignments

    fprintf('*** Cleaning the flashes... ***\n')
    tic
    for k = 1:P.nMice
        avgIntensity = db(k).motSVD(:,1);
        [~,t] = max(abs(avgIntensity));
        avgIntensity = avgIntensity*sign(avgIntensity(t));
        intensMed = median(avgIntensity);
        intensMin = min(avgIntensity);

        vidIntensThresh = intensMin+(intensMed-intensMin)*[0.15 0.25];
        [~, ~, intensDown] = schmittTimes(1:numel(avgIntensity), avgIntensity, vidIntensThresh);
        success = 0;
        if numel(intensDown)==2 && (diff(intensDown)>1e4)
            success = 1;
        elseif (numel(intensDown)==4) && ((intensDown(2)-intensDown(1)<15) & (intensDown(4)-intensDown(3)<15))
            intensDown = intensDown([1 3]);
            success = 1;
        else
            fprintf('Didn''t work for mouse %s', P.mouseRef{k})
        end

        if success
            % replace with median values
            idx = [intensDown(1)-1:intensDown(1)+10 intensDown(2)-1:intensDown(2)+10];
            db(k).motSVD(idx,:) = repmat(median(db(k).motSVD),[numel(idx),1]);
        end
    end
    toc
    fprintf('*** Done. ***\n')
    % all the other won't be smoothed

    %% get time binned firing rate for each cell locked onto video frames

    fprintf('*** Binning... ***\n')
    neurAll = cell(1,P.nMice);
    neurProj = cell(1,P.nMice);
    motDat = cell(1,P.nMice);
    baseIdx = cell(1,P.nMice);
    for k = 1:P.nMice
        fprintf('Binning for mouse %s...', P.mouseRef{k})
        tic

        % Get spike histogram with video frame timings
        neurAll{k} = nan(numel(db(k).videoTimestampsEphysTime_face)-1,numel(db(k).C.CluID));
        for c = 1:numel(db(k).C.CluID)
            neurAll{k}(:,c) = histcounts(db(k).sp.st(db(k).sp.clu == db(k).C.CluID(c)),db(k).videoTimestampsEphysTime_face);
        end

        % Padd with 1st repeat
        neurAll{k} = cat(1,neurAll{k}(1,:),neurAll{k});

        % Find baseline period
        baseIdx{k} = ones(1,numel(db(k).videoTimestampsEphysTime_face));
        for e = 1:numel(db(k).eventTimes)
            idxTri = db(k).videoTimestampsEphysTime_face > db(k).eventTimes(e) & db(k).videoTimestampsEphysTime_face < db(k).eventTimes(e) + 4.5;
            baseIdx{k}(idxTri) = 0;
        end
        baseIdx{k} = logical(baseIdx{k});

        % Get proj of neural data on focusOn subspace
        [~, M, ~, ~, sigma] = getComponents(db(k).spikeData, db(k).audiovisuoCode, P);
        neurProj{k} = projectData((neurAll{k}-nanmean(neurAll{k}))./squeeze(sigma)',M);

        % Get mot data
        motDat{k} = db(k).motSVD(:,:,:,1);

        fprintf('Done in %ss.\n',toc)
    end
    fprintf('*** Binning done. ***\n')
end