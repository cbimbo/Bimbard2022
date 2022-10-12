%% Look at LFP

localLFPFolder = '\\zserver\Lab\Share\Celian\audioVisLFP';
clear LFP CSD depthChan timebinsAll
for k = 1:P.nMice
    rawFolder = fileparts(P.fileref.ksDir{k});
    [~,recRef] = fileparts(rawFolder);

    finalLFPFile = fullfile(localLFPFolder,recRef,[recRef, '_lfp_triggered.mat']);
    load(finalLFPFile);

    timebinsAll{k} = timebins;
    stimMap = stimMap - nanmean(stimMap(:, timebins<0, :), [2 3]);
    
    % Get channel maps
    xposUni = unique(chanPos(:,1));
    clear depthChanTmp LFPTmp
    for xposIdx = 1:numel(xposUni)
        xposChanIdx = find(chanPos(:,1) == xposUni(xposIdx));
        [depthChanTmp{xposIdx},sortChanIdx] = sort(chanPos(xposChanIdx,2));

        LFPTmp{xposIdx} = nanmean(stimMap(xposChanIdx(sortChanIdx),:,:),3);
    end

    gw = gausswin(5,1);
    smWin = gw./sum(gw);
    xposCeil = ceil((xposUni+1)/250)*250;
    xposCeilUni = unique(xposCeil);
    for xposIdx = 1:numel(xposCeilUni)
        depthChan{k}{xposIdx} = nanmean(cat(3,depthChanTmp{xposCeil == xposCeilUni(xposIdx)}),3);
        corresX = find(xposCeil == xposCeilUni(xposIdx));
        for xx = 1:numel(corresX)
            LFPTmp{corresX(xx)} = interp1(depthChanTmp{corresX(xx)},LFPTmp{corresX(xx)},depthChan{k}{xposIdx});
        end
        LFP{k}{xposIdx} = nanmean(cat(3,LFPTmp{xposCeil == xposCeilUni(xposIdx)}),3);
        ba = cat(1, repmat(LFP{k}{xposIdx}(1,:),[2 1]), LFP{k}{xposIdx}, repmat(LFP{k}{xposIdx}(end,:),[2 1]));
        bas = conv2(smWin,1,ba, 'same');
        LFP{k}{xposIdx} = bas(3:end-2,:);
    end

    for xposIdx = 1:numel(xposCeilUni)
        CSD{k}{xposIdx} = CSD_construction(LFP{k}{xposIdx},depthChan{k}{xposIdx},0);
    end
end

%% Plot

clear granBoundManual
for k = 1:P.nMice
    % Plot CSD
    plotLFPandCSD(LFP{k},CSD{k},depthChan{k},timebinsAll{k}, ...
        db_bu(k).C.XPos,db_bu(k).C.Depth,P.fileref.anatDir{k});

    for xx = 1:numel(LFP{k})
        tmp = ginput(2);
        granBoundManual{k}{xx} = tmp(:,2);
        % double click if unsure, will nan it
        if abs(diff(granBoundManual{k}{xx}))<10
            granBoundManual{k}{xx} = [nan; nan];
        end
    end
end


%% save

save(fullfile(preprocFolder,sprintf('granularLayerBounds_ephys%s',what2Load)),'granBoundManual')
