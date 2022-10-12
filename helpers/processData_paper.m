function [db,P] = processData_paper(db,P,recompute)

    if ~exist('recompute','var')
        recompute = 0;
    end

    % clean fields in case
    fnames = fieldnames(db);
    db = rmfield(db,fnames(contains(fnames,'Data')));

    P.bins = P.proc.window(1)+P.proc.binSize/2:P.proc.binSize:P.proc.window(2);
    P.nBins = numel(P.bins);

    fprintf('*** Preprocessing data... ***\n')
    for k = 1:P.nMice
        fprintf('Preprocessing mouse %s...',P.mouseRef{k})

        % minimal workspace should be saved somewhere
        if ~isempty(db(k).imecX)
            preproc_savingFileName = [P.pathSaving '\processedData\' P.exp '\' db(k).subject '_' db(k).date '_imec' num2str(db(k).imecX) '_preproc.mat'];
        else
            preproc_savingFileName = [P.pathSaving '\processedData\' P.exp '\' db(k).subject '_' db(k).date '_preproc.mat'];
        end
    
        if ~recompute && exist(preproc_savingFileName,'file')
            fprintf('\nMouse %s already preprocessed. Loading it.\n',P.mouseRef{k})
            load(preproc_savingFileName);
        else
            %% Bin spike and motion data
            
            nTrials = numel(db(k).eventTimes);
            nBins = numel(P.bins);

            % init
            spikeData = nan(nBins,numel(db(k).C.CluID),nTrials);
            motionData = nan(nBins,size(db(k).motSVD{2},2),nTrials); 
            pupilareaData = nan(nBins,1,nTrials);
            pupilcomData = nan(nBins,2,nTrials);
            pupilmotData = nan(nBins,2,nTrials);
            blinkData = nan(nBins,1,nTrials);

            % Get spikes psths
            for c = 1:numel(db(k).C.CluID)

                % spikes
                clu = db(k).C.CluID(c);
                st = db(k).sp.st(db(k).sp.clu == clu);

                % get psth
                [ba, bins] = getPSTH(st, db(k).eventTimes, P.proc.window, P.proc.binSize);

                spikeData(:,c,:) = ba'/P.proc.binSize;
            end

            if P.smooth
                % smooth each trial with a gaussian filter
                gw = gausswin(10,P.smooth);
                gw(1:numel(gw)/2) = 0; % causal filter
                smWin = gw./sum(gw);
                for tt = 1:size(spikeData,3)
                    ba = spikeData(:,:,tt);
                    ba = cat(1, repmat(ba(1,:),[2 1]), ba, repmat(ba(end,:),[2 1]));
                    bas = conv2(smWin,1,ba, 'same');
                    spikeData(:,:,tt) = bas(3:end-2,:);
                end
            end

            % Get movements psths
            for tr = 1:nTrials
                for b = 1:nBins
                    idxBin = (db(k).videoTimestampsEphysTime_face > db(k).eventTimes(tr)+bins(b)-P.proc.binSize/2) & ...
                        (db(k).videoTimestampsEphysTime_face < db(k).eventTimes(tr)+bins(b)+P.proc.binSize/2); % get time stamps for this trial
                    tmp = mean(db(k).motSVD{2}(idxBin,:),1);
                    motionData(b, 1:numel(tmp), tr) = tmp;
                end
            end

            % Get eye related variables psths
            for tr = 1:nTrials
                for b = 1:nBins
                    idxBin = (db(k).videoTimestampsEphysTime_eye > db(k).eventTimes(tr)+bins(b)-P.proc.binSize/2) & ...
                        (db(k).videoTimestampsEphysTime_eye < db(k).eventTimes(tr)+bins(b)+P.proc.binSize/2); % get time stamps for this trial
                    pupilareaData(b, 1, tr) = mean(db(k).pupilarea(idxBin));
                    pupilcomData(b,:,tr) = mean(db(k).pupilcom(idxBin,:));
                    pupilmotData(b,:,tr) = mean(db(k).pupilmot(idxBin,:));
                    blinkData(b, 1, tr) = mean(db(k).blink(idxBin));
                end
            end

            %% Little additional processing for the eye data

            if P.miceweye(k)

                % pupil related info:
                % X axis: >0 goes towards nose, <0 towards back
                % Y axis: >0 goes down, <0 goes up

                if strcmp(db(k).subject,'CB003')
                    blinkData(:) = detrend(blinkData);
                end

                % baseline correct pupil
                pupilareaData = pupilareaData - nanmean(pupilareaData(P.bins<0,:,:));
                pupilcomData = pupilcomData - nanmean(pupilcomData(P.bins<0,:,:));

                % reorient Y axis of pupil
                pupilcomData(:,2,:) = - pupilcomData(:,2,:);

                pupilmotData = diff(pupilcomData);
                pupilmotData = cat(1,pupilmotData(1,:,:),pupilmotData);
            end

            % save it
            save(preproc_savingFileName,'spikeData', ...
                'motionData', ...
                'pupilareaData','pupilcomData','pupilmotData','blinkData')
        end

        % save it in db for output
        db(k).spikeData = spikeData;
        db(k).motionData = motionData;
        db(k).pupilareaData = pupilareaData;
        db(k).pupilcomData = pupilcomData;
        db(k).pupilmotData = pupilmotData;
        db(k).blinkData = blinkData;
        fprintf('Done.\n')
    end
    fprintf('*** Preprocessing done. ***\n')
end