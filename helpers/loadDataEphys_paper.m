function [C, audiovisuoCode, P, ...
        eventTimes, sp, motSVD, videoTimestampsEphysTime_face, uMotMask, ...
        pupilarea, pupilcom, blink, videoTimestampsEphysTime_eye] = ...
        loadDataEphys_paper(subject, date, imecX, P, recompute)

    if ~exist('recompute','var')
        recompute = 0;
    end

    % minimal workspace should be saved somewhere
    dataPath = fullfile(P.pathSaving,'processedData',P.exp);
    if ~isempty(imecX)
        expRef = [subject '_' date '_imec' num2str(imecX)];
    else
        expRef = [subject '_' date];
    end
    
    if ~recompute && exist(fullfile(dataPath,expRef))
        [C, audiovisuoCode, P, ...
        eventTimes, sp, motSVD, videoTimestampsEphysTime_face, uMotMask, ...
        pupilarea, pupilcom, blink, videoTimestampsEphysTime_eye] = loadONE(dataPath, expRef);
    else
        %% Note
        % This part will only work on my computer. It loads the very raw
        % data and puts everything together.
        % Actual raw data (ephys file, spikesorting, etc.) are available
        % upon reasonable request.

        %% Define folders
        subjectsFolder = '\\znas.cortexlab.net\Subjects';
        if ~exist(fullfile(subjectsFolder, subject, date),'dir')
            subjectsFolder = '\\zinu\Subjects';
        end
        alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
        if ~exist(alignDir,'dir')
            mkdir(alignDir) 
        end

        root = getRootDir(subject, date);
        
        %% Get basic info
        
        [tags, hasEphys] = getEphysTags(subject, date);
        
        % focus on audioVis
        tags = tags(cellfun(@(x) ~isempty(x), strfind(tags,P.exp)));
        
        % determine what exp nums exist
        [expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
            whichExpNums(subject, date); % used to be dat.whichExpNums?
        
        useFlipper = true; % no idea why empty
        
        %% align times (timeline to ephys)
        
        % for any ephys, load the sync data
        if hasEphys
            for t = 1:length(tags)
                if isempty(tags{t})
                    [~, pdFlips, allET] = loadSyncChronic(subject, date);
                else
                    [~, pdFlips, allET] = loadSyncChronic(subject, date, tags{t});
                end
                if useFlipper
                    ephysFlips{t} = allET; % problem with big files, had to bypass spikeGLXdigitalParse
                    %             ephysFlips{t} = allET{7}{1};
                else
                    ephysFlips{t} = pdFlips;
                end
            end
        end
        
        %%% what if sequential recordings in the same day? -- should be
        %%% fixed with the tag selection
        % synchronize multiple ephys to each other
        if hasEphys
            if length(tags)>1
                for t2 = 2:length(tags)
                    fprintf(1, 'correct ephys %s to %s\n', tags{t2}, tags{1});
                    [~, b] = makeCorrection(ephysFlips{1}, ephysFlips{t2}, false);
                    writeNPY(b, fullfile(alignDir, sprintf('correct_ephys_%s_to_ephys_%s.npy', tags{t2}, tags{1})));
                end
            end
        end
        
        % detect sync events from timelines
        tlFlips = {};
        for e = 1:length(expNums)
            if hasTimeline(e)
                Timeline = tl{e};
                tt = Timeline.rawDAQTimestamps;
                if useFlipper
                    evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'flipper'));
                    evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
                else
                    evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                    evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
                    evT = evT([true; diff(evT)>0.2]);
                end
                tlFlips{e} = evT;
            end
        end
        
        % match up ephys and timeline events:
        % algorithm here is to go through each timeline available, figure out
        % whether the events in timeline align with any of those in the ephys. If
        % so, we have a conversion of events in that timeline into ephys
        %
        % Only align to the first ephys recording, since the other ones are aligned
        % to that
        if hasEphys
            ef = ephysFlips{1};
            if useFlipper && ef(1)<0.001
                % this happens when the flipper was in the high state to begin with
                % - a turning on event is registered at the first sample. But here
                % we need to drop it.
                ef = ef(2:end);
            end
            for e = 1:length(expNums)
                if hasTimeline(e)
                    fprintf('trying to correct timeline %d to ephys\n', expNums(e));
                    %Timeline = tl{e};
                    tlT = tlFlips{e};
                    
                    success=false;
                    if length(tlT)==length(ef)
                        % easy case: the two are exactly coextensive
                        [~,b] = makeCorrection(ef, tlT, true);
                        success = true;
                    elseif length(tlT)<length(ef) && ~isempty(tlT)
                        [~,b,success] = findCorrection(ef, tlT, false);
                    elseif length(tlT)>length(ef) && ~isempty(tlT)
                        [~,a,success] = findCorrection(tlT, ef, false);
                        if ~isempty(a)
                            b = [1/a(1); -a(2)/a(1)];
                        end
                    end
                    if success
                        writeNPY(b, fullfile(alignDir, ...
                            sprintf('correct_timeline_%d_to_ephys_%s.npy', ...
                            expNums(e), tags{1})));
                        fprintf('success\n');
                        eTimeline2keep = e;
                    else
                        fprintf('could not correct timeline to ephys\n');
                    end
                end
            end
        end
        
        TLexp = expNums(eTimeline2keep);
        Timeline = tl{eTimeline2keep};
        P.fileref.rootDir = fullfile(root, num2str(expNums(eTimeline2keep)));

        %% align times (timeline to blocks)
        
        %%% alignment doesn't work because block.stimWindowUpdateTimes is empty.
        %%% also video-0 is actually not displaying any flip in photodiod.
        
        %%% there seems to be a lot of jitter during the first presentation
        %%% round between the onsets of trials, as blocks.events and photodiode are
        %%% not well aligned.
        
        switch P.exp
            case 'audioVis'
                % get sound waveforms
                % get sound names first
                video_list = [25 22 11 18 17 26 15 0 27 29 9 7]; % hard coded this, should be able to retrieve it from somewhere... but hard if don't know the files yet
                labelsGroupsVid = unique(video_list);
                labelsGroupsSnd = labelsGroupsVid;
                nGroupsVid = numel(video_list);
                nGroupsSnd = nGroupsVid;
                sndPath = '\\zserver\Code\Rigging\ExpDefinitions\Tim\filmWorldwGaps\30-2019-11-22-ratio-edited\processed-merged-4-sec-w-corner-150';

                envSnd = zeros(Timeline.hw.daqSampleRate*4.022,nGroupsVid);
                for s = 1:nGroupsVid
                    [audio, Fs] = audioread(fullfile(sndPath, ['video-0-audio-' num2str(labelsGroupsVid(s)) '.mp4']));
                    %                 sndTime = (0:size(audio(:,1),1)-1)/Fs;
                    % take the envelope
                    e = envelope(audio,2000,'rms'); % need to adjust
                    % resample to match timeline daq sample
                    envSnd(:,s) = resample(e(:,1),Timeline.hw.daqSampleRate,Fs);
                    clear e
                end

                %% alignments
                % didn't get photodiode flips above, so get them now
                tt = Timeline.rawDAQTimestamps;
                evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));    
                sndOutput = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'audioMonitor')); 

                % find thresholds
                %%% there should be a way to simplify all
                %%% this...
                mini = min(evTrace((tt>tt(end)*0.5) & (tt<tt(end)*0.6))); % focus on 10% of time in the middle
                maxi = max(evTrace((tt>tt(end)*0.5) & (tt<tt(end)*0.6)));
                pdTOn = schmittTimes(tt, evTrace, [mini+1/4*(maxi-mini) maxi-0.05*(maxi-mini)]); % should contain onsets
                pdTOff = schmittTimes(tt, evTrace, [mini+0.05*(maxi-mini) mini+1/2*(maxi-mini)]); % should contain offsets / not sure it's needed

                pdTOn([find(diff(pdTOn)<1); find(diff(pdTOn)<1)+1]) = [];
                pdTOff([find(diff(pdTOff)<1); find(diff(pdTOff)<1)+1]) = [];
                if pdTOn(2) -pdTOn(1) > 10
                    % weird screen actualization?
                    pdTOn(1) = [];
                    pdTOff(1) = [];
                end

                if pdTOff(2) - pdTOff(1) > 10
                    % can still have something..
                    pdTOff(1) = [];
                end

                if pdTOn(1) - pdTOff(1) < 0.01
                    % weird screen actualization?
                    pdTOn(1) = [];
                    pdTOff(1) = [];
                end

                % other weird screen fluctuations before
                idx = find(pdTOn < pdTOff(1));
                pdTOn(idx(1:end-2)) = []; % first trial and first crossing

                if ~(pdTOn(2)-pdTOn(1)>1.1 || pdTOn(2)-pdTOn(1)<1.2)
                    pdTOn = pdTOn(2:end); % first crossing is screen starting?
                end

                pdTOn_bu = pdTOn(1:2:end);

                for e = 1:length(expNums)
                    %%% Should save pdT only when successfully aligned
                    if hasBlock(e) && isfield(blocks{e}.events, 'moviesStructValues')
                        fprintf('trying to correct block %d to timeline %d\n', expNums(e), expNums(eTimeline2keep));

                        nTrials = length(blocks{e}.events.playingMovieTimes);

                        trialsWVid = cell2mat(cellfun(@(x) ~contains(x, 'video-0'), {blocks{e}.events.moviesStructValues(1:nTrials).name}, 'uni', 0));
                        intermTrialNum = sum(trialsWVid);
                        if numel(pdTOn_bu) == intermTrialNum
                            disp('Interm step looks ok...')
                            pdTOn = pdTOn_bu;
                        elseif numel(pdTOn_bu) == intermTrialNum+1
                            disp(['Not the correct number of video only trials. Check if it''s the last trial? ' ...
                                'Remove it and see.'])
                            pdTOn = pdTOn_bu(1:end-1);
                        else
                            fprintf('Something weird')
                            continue
                        end
%                         figure; scatter(diff(blocks{e}.events.playingMovieTimes(trialsWVid)),diff(pdTOn))

                        pdT = pdTOn;

                        % add here the missing events
                        for trialNum = find(~trialsWVid)

                            % extract local envelope and fit the appropriate
                            % sound envelope to find beginning
                            lenMovie = 4; % works
                            idxEndPrevTrial = find(tt>pdT(trialNum-1)+lenMovie,1);
                            sndOutputSnippet = Timeline.rawDAQData(idxEndPrevTrial:idxEndPrevTrial + Timeline.hw.daqSampleRate*8, ...
                                strcmp({Timeline.hw.inputs.name}, 'audioMonitor')); % 8 is arbitrary here
                            env = envelope(sndOutputSnippet,100,'rms'); % need to adjust

                            % get sound id
                            tmp = regexp(blocks{e}.events.moviesStructValues(trialNum).name,'\d+\d*','match');

                            if ~strcmp(tmp{2},'0')
                                sndEnvTemplate = envSnd(:,labelsGroupsSnd == str2num(tmp{2})); % should have sampling rate 2500Hz

                                % find best correlation
                                c = xcorr2(env,sndEnvTemplate); % using this because fast?
                                [~,idx] = max(c); % will give the end point

                                sndStart = (idx-numel(sndEnvTemplate))/Timeline.hw.daqSampleRate - 0.135; % due to a fixed delay between the video and audio
                            else
                                sndStart = 2; % doesn't matter much, since no audio no video
                            end

                            pdT = [pdT(1:(trialNum-1)); pdT((trialNum-1))+lenMovie+sndStart; pdT(1+(trialNum-1):end)];
                        end

                        %%% supposes that it's played only once
                        if numel(pdT) == nTrials
                            % check whether the block is the correct one
                            eBlock2keep = e;
                            disp('success!')
                            % should be perfectly correlated
%                             figure; scatter(diff(blocks{e}.events.playingMovieTimes),diff(pdT))
                            pdTsuccess = pdT;
                        elseif numel(pdT) < nTrials
                            error('Missed one trial?')
                        elseif numel(pdT) > nTrials
                            error('Too many photodiode detected?')
                        end

                        % plot photodiode for checking
                        %%% last trial should be missing...
%                         figure; hold all
%                         plot(tt,zscore(evTrace));
%                         env = envelope(sndOutput,100,'rms'); % need to adjust
%                         plot(tt,zscore(env))
%                         scatter(pdT, ones(1, length(pdT)),80,'g')

                        % plot the photodiode and sndoutput locked for each category
                        % takes a bit of time but useful to compare
                        win = [-0.5 4.5];
                        photoDiode = nan(diff(win)*Timeline.hw.daqSampleRate,length(blocks{e}.events.playingMovieTimes));
                        soundOutput = nan(diff(win)*Timeline.hw.daqSampleRate,length(blocks{e}.events.playingMovieTimes));
                        for trialNum = 1:length(blocks{e}.events.playingMovieTimes)
                            winIdx = (tt>pdT(trialNum)+win(1)) & (tt<pdT(trialNum)+win(2));
                            photoDiode(1:sum(winIdx),trialNum) = evTrace(winIdx);
                            soundOutput(1:sum(winIdx),trialNum) = envelope(sndOutput(winIdx),100,'rms');
                        end

                        figure;
                        [~,sortidx] = sort({blocks{e}.events.moviesStructValues(1:nTrials).name});
                        ax(1) = subplot(221);
                        imagesc(photoDiode(:,sortidx)')
                        ylabel('Sorted by video')
                        title('Photodiode')
                        ax(2) = subplot(222);
                        imagesc(soundOutput(:,sortidx)')
                        title('Sound output')
                        linkaxes(ax(1:2),'y')
                        tmp = regexp({blocks{e}.events.moviesStructValues(1:nTrials).name},'\d-a','split');
                        [~,sortidx] = sort(cellfun(@(x) [x{2} x{1}], tmp, 'uni',0));
                        ax(3) = subplot(223);
                        imagesc(photoDiode(:,sortidx)')
                        ylabel('Sorted by audio')
                        ax(4) = subplot(224);
                        imagesc(soundOutput(:,sortidx)')
                        linkaxes(ax(3:4),'y')
                    end
                end
                if isempty(eBlock2keep)
                    error('Couldn''t find block file...')
                end   
                   
            case 'audioVisSimple'
                % maybe no need to rewrite the whole thing but easier to see
                % clearly
                for e = 1:length(expNums)
                    if hasBlock(e) && isfield(blocks{e}.events, 'soundFreqValues')
                        fprintf('trying to correct block %d to timeline %d\n', expNums(e), expNums(eTimeline2keep));
                        if useFlipper
                            % didn't get photodiode flips above, so get them now
                            tt = Timeline.rawDAQTimestamps;
                            evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                            
                            % find thresholds
                            tmp = medfilt1(evTrace,10);
                            tmp = diff(medfilt1(abs(diff(tmp)),10));
                            pdT = tt(find(diff(tmp>0.01)==1) + 3);
                            pdT = pdT([false (diff(pdT)>0.4)]);
                            pdT(1) = [];
                            pdT = pdT([true (diff(pdT)<1.4)]);
                            idxr = find(diff(pdT)<0.8,1)+1;
                            while ~isempty(idxr)
                                pdT(idxr) = [];
                                idxr = find(diff(pdT)<0.8,1)+1;
                            end
                        else
                            pdT = tlFlips{eTimeline2keep};
                        end
                        
                        nTrials = length(blocks{e}.events.newTrialTimes);
                        if nTrials == numel(pdT) 
                            eBlock2keep = e;
                            disp('success!')
                        elseif nTrials-1 == numel(pdT)
                            eBlock2keep = e;
                            disp('There''s one trial difference, but let''s call it a success!')
                        end
                    end
                end
                
            case 'sparseNoise'
                % not sure I have to do a new one
                for e = 1:length(expNums)
                    if hasBlock(e) && contains(blocks{e}.expDef,P.exp)
                        fprintf('trying to correct block %d to timeline %d\n', expNums(e), expNums(eTimeline2keep));
                        if useFlipper
                            % didn't get photodiode flips above, so get them now
                            tt = Timeline.rawDAQTimestamps;
                            evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                            
                            pdT = schmittTimes(tt, evTrace, [5 6]);
                        else
                            pdT = tlFlips{eTimeline2keep};
                        end
                        
                        nTrials = length(blocks{e}.stimWindowUpdateTimes);
                        if nTrials == numel(pdT) % don't know why -1???
                            % check whether the block is the correct one
                            eBlock2keep = e;
                            disp('success!')
                        elseif abs(nTrials-numel(pdT)) < 3
                            eBlock2keep = e;
                            disp('There''s a few trial difference, but let''s call it a success!')
                        end
                    end
                end
        end
        block = blocks{eBlock2keep};
        P.fileref.blockDir = fullfile(root, num2str(expNums(eBlock2keep))); 

        %% get spike data

        sp = loadAllKsDir(subject, date, P.exp);
        P.fileref.ksDir = sp.ksDir;

        % quick hack to select only one probe
        % suboptimal because has to recompute the same thing twice
        if numel(sp)>1
            if ~isempty(imecX)
                sp = sp(imecX+1);
            else
                error('choose which probe to look at; or modify the script')
            end
        end

        if all(sp.cgs == 3)
            tmp = tdfread(fullfile(sp.ksDir,'cluster_KSLabel.tsv'));
            sp.cgs = zeros(numel(tmp.KSLabel(ismember(tmp.cluster_id,sp.cids))),1);
            sp.cgs(strfind(tmp.KSLabel(ismember(tmp.cluster_id,sp.cids))','g')) = 3.2;
            sp.cgs(strfind(tmp.KSLabel(ismember(tmp.cluster_id,sp.cids))','m')) = 3.1;
        end

        %% get visual responses

        %%% Blocks aren't aligned so just keep the onset of each trial.
        bTLtoMaster = readNPY(fullfile(alignDir, sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{P.probeNumber})));

        eventTimes = pdTsuccess; % get all trial onsets
        eventTimes = applyCorrection(eventTimes, bTLtoMaster);

        % get stim code
        nTrials = length(block.events.newTrialTimes);
        if nTrials>length(eventTimes)
            % setting correct number of trials
            nTrials = length(eventTimes);
            disp(['Trial number in block and photodiode events are different. Keep the photodiode events (' num2str(nTrials) ').'])
        end

        switch P.exp
            case 'audioVis'
                audiovisuoCode = zeros(2,nTrials);
                for trialNum = 1:nTrials
                    tmp = regexp(block.events.moviesStructValues(trialNum).name,'\d+\d*','match');
                    audiovisuoCode(1,trialNum) = str2double(tmp{1});
                    audiovisuoCode(2,trialNum) = str2double(tmp{2});
                end
            case 'audioVisSimple'
                clear amp
                for trialNum = 1:nTrials
                    audiovisuoCode(1,trialNum) = block.events.gOriValues(trialNum);
                    audiovisuoCode(2,trialNum) = block.events.soundFreqValues(trialNum);
                    % find amplitude
                    amp(trialNum) = block.paramsValues(trialNum).audAmplitude;
                end
                % find discrete amplitude
                for f = unique(audiovisuoCode(2,:))
                    if ~ismember(f,[666])
                        idx = find(audiovisuoCode(2,:) == f);
                        b = unique(amp(idx));
                        for i = 1:4
                            audiovisuoCode(3,idx(amp(idx) == b(i))) = i*0.1;
                        end
                    end
                end
                audiovisuoCode(1,audiovisuoCode(1,:) == 1) = -1; % change blank to -1
                audiovisuoCode(2,audiovisuoCode(2,:) == 666) = -1; % change blank to -1
                audiovisuoCode(2,audiovisuoCode(2,:) == 123) = 100000; % change noise to a lot
                audiovisuoCode(2,:) = sum(audiovisuoCode(2:3,:),1);
                audiovisuoCode(3,:) = [];
                nGroupsSnd = numel(unique(audiovisuoCode(2,:)));
                nGroupsVid = numel(unique(audiovisuoCode(1,:)));
                labelsGroupsSnd = unique(audiovisuoCode(2,:));
                labelsGroupsVid = unique(audiovisuoCode(1,:));
            otherwise
                nTrials = length(pdT);
                audiovisuoCode = [];
                nGroupsVid = [];
                nGroupsSnd = [];
                labelsGroupsVid = [];
                labelsGroupsSnd = [];
        end

        %% get trials numbers

        switch P.exp
            case {'audioVis','audioVisSimple'}
                trialNumPerStim = getTrialNumPerStim(audiovisuoCode);
            otherwise
                trialNumPerStim = [];
        end

        %% get clusters depth

        [h,b] = hist(double(sp.clu),double(unique(sp.clu)));
        CluSpknum = h;

        clusterList = unique(sp.clu);
        nClusters = length(clusterList);

        ycoords = sp.ycoords;
        xcoords = sp.xcoords;
        pcFeat = sp.pcFeat;
        pcFeat = squeeze(pcFeat(:,1,:)); % take first PC only
        pcFeat(pcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there.
        pcFeatInd = sp.pcFeatInd;
        spikeTemps = sp.spikeTemplates;

        % which channels for each spike?
        spikeFeatInd = pcFeatInd(spikeTemps+1,:);
        % ycoords of those channels?
        spikeFeatYcoords = ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12
        % center of mass is sum(coords.*features)/sum(features)
        spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);

        spikeFeatXcoords = xcoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12
        spikeXPos = sum(spikeFeatXcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);

        for cl = clusterList'
            C.XPos(clusterList == cl) = nanmean(spikeXPos(sp.clu == cl)); % not sure why there can be nans here
            C.Depth(clusterList == cl) = nanmean(spikeDepths(sp.clu == cl));
        end
        C.CluID = clusterList';
        C.CluLab = sp.cgs;
        C.CluSpknum = CluSpknum;

        %% now get video

        expNum = expNums(eTimeline2keep);

        % get the real time at which video was played
        realtimeFileName_face = fullfile(root, num2str(expNum), 'face_timeStamps.mat');
        if ~exist(realtimeFileName_face, 'file')
            % check video file is here
            alignFace = 0;
            s = dir(fullfile(root, num2str(expNum), 'face.mj2'));
            if ~isempty(s) && (s.bytes>10^6)
                try
                    alignVideo(subject, date, expNum, 'face')
                    alignFace = 1;
                catch
                    print('Couldn''t align face')
                end
            else
                fprintf('Couldn''t find video file: %s\n', realtimeFileName_face)
            end
        else
            alignFace = 1;
        end

        procFileFace = fullfile(root, num2str(expNum), 'face_proc.mat');
        if exist(realtimeFileName_face, 'file') & exist(procFileFace,'file') & alignFace
            realtimeFile = load(realtimeFileName_face);
            realwindowFile = load(fullfile(root, num2str(expNum), 'face_avgIntensity.mat'));
            videoTimestampsEphysTime_face = applyCorrection(realtimeFile.tVid, bTLtoMaster);

            % load file processed by facemap
            load(procFileFace)

            motSVD = proc.motSVD;
            uMotMask = proc.uMotMask; clear proc

        else
            motSVD = [];
            videoTimestampsEphysTime_face = [];
            uMotMask = [];
        end

        %% get pupil

        % get the real time at which video was played
        realtimeFileName_eye = fullfile(root, num2str(expNum), 'eye_timeStamps.mat');
        if ~exist(realtimeFileName_eye, 'file')
            % check video file is here
            s = dir(fullfile(root, num2str(expNum), 'eye.mj2'));
            alignEye = 0;
            if ~isempty(s) && (s.bytes>10^6)
                try
                    alignVideo(subject, date, expNum, 'eye')
                    alignEye = 1;
                catch
                    print('Couldn''t align eye')
                end
            else
                fprintf('Couldn''t find video file: %s', realtimeFileName_eye)
            end
        else
            alignEye = 1;
        end

        procFileEye = fullfile(root, num2str(expNum), 'eye_proc.mat');
        if exist(realtimeFileName_eye, 'file') & exist(procFileEye,'file') & alignEye
            realtimeFile = load(realtimeFileName_eye);
            videoTimestampsEphysTime_eye = applyCorrection(realtimeFile.tVid, bTLtoMaster);

            % load file processed by facemap
            load(procFileEye);

            pupilarea = proc.pupil.area;
            pupilcom = proc.pupil.com;
            pupilcom = medfilt1(pupilcom,5,'truncate');
            pupilmot = abs(diff(pupilcom)); pupilmot = [pupilmot(1,:); pupilmot];
            blink = proc.blink.area;

            threshBlink = nanzscore(blink)<-10; % arbitrary here but seems to fit quite well when manually checking the videos
            X = ~(threshBlink)';
            Y = cumsum(X-diff([1,X])/2);
            pupilarea = interp1(1:nnz(X),pupilarea(X),Y);
            pupilcom(:,1) = interp1(1:nnz(X),pupilcom(X,1),Y);
            pupilcom(:,2) = interp1(1:nnz(X),pupilcom(X,2),Y);
            pupilmot(:,1) = interp1(1:nnz(X),pupilmot(X,1),Y);
            pupilmot(:,2) = interp1(1:nnz(X),pupilmot(X,2),Y);

            P.miceweye = 1;
        else
            pupilarea = [];
            pupilcom = [];
            pupilmot = [];
            blink = [];
            videoTimestampsEphysTime_eye = [];

            P.miceweye = 0;
        end

        P.fileref.blockDir = fullfile(root, num2str(expNums(eBlock2keep)));
        P.fileref.rootDir = fullfile(root, num2str(expNums(eTimeline2keep)));
        P.fileref.ksDir = sp.ksDir;
        %% get anat

        anatFolder = dir(fullfile(P.fileref.rootDir(1:strfind(P.fileref.rootDir,subject)+4),'histology'));
        if any(strcmpi({anatFolder.name},'slices')) % classic anatomy.
            %%% Old data, only one shank? Otherwise will have to include the x
            %%% coordinate...
            anatFolderClassic = fullfile(anatFolder(1).folder,'slices');

            % get aligned depth and corresponding areas
            load(fullfile(anatFolderClassic,'probe_ccf.mat'))
            C.area = interp1(3840-probe_ccf.probe_depths,probe_ccf.trajectory_areas,C.Depth,'nearest'); % works if bank0

            P.fileref.anatDir = anatFolderClassic;
        else % Anatomy hasn't been done through this pipeline, but manually aligned. Hardcoded threshold.
            warning('Can''t find processed anatomy for subject %s. Putting a hard-coded threshold.',subject)
            switch subject
                % everything here is done 'by hand' due to issues with
                % alignment programs
                case 'AL031'
                    depththreVIS = [4300,inf];
                    C.area(C.Depth>depththreVIS(1) & C.Depth<depththreVIS(2)) = 186; % VISp
                    C.area(C.Depth<depththreVIS(1)) = 455; % very likely subiculum, but 'HPF' for now
                case 'CB008'
                    depththreVIS = [2800,inf];
                    C.area(C.Depth>depththreVIS(1) & C.Depth<depththreVIS(2)) = 186; % VISp
                    C.area(C.Depth<depththreVIS(1)) = 455; % very likely subiculum, but 'HPF' for now
                case 'CR_transectomy7'
                    if imecX == 1
                        % on the cut side (ipsi, imec0), most medial shank is
                        % in posteromedial visual area VISpm
                        xposthreVIS = [100 inf]; % most medial is shank 0
                        C.area(C.XPos>xposthreVIS(1) & C.XPos<xposthreVIS(2)) = 186; % VISp
                        C.area(C.XPos<xposthreVIS(1)) = 200; % VISpm
                    else
                        % all the rest is in VISp
                        C.area(C.Depth>=0) = 186;
                    end
                case 'CR_transectomy8'
                    if imecX == 0
                        % all in lateral visual area (VISl)
                        C.area(C.Depth>=0) = 179; % VISl
                    else
                        % all in VISp
                        C.area(C.Depth>=0) = 186;
                    end
                otherwise
                    C.area(C.Depth>=0) = 186;
            end
            P.fileref.anatDir = [];
        end

        %% finish

        % put everything in parameters
        P.nGroupsVid = nGroupsVid;
        P.nGroupsSnd = nGroupsSnd;
        P.labelsGroupsVid = labelsGroupsVid;
        P.labelsGroupsSnd = labelsGroupsSnd;
        P.miceweye = logical(P.miceweye); % for some reason doesn't work when inputing directly true or false...

        % save data
        if ~exist(dataPath, 'dir')
            mkdir(dataPath)
        end

        dat.audiovisuoCode = audiovisuoCode;
        dat.eventTimes = eventTimes;
        dat.P = P;
        dat.sp = sp;
        dat.C = C;
        dat.blink = blink;
        dat.pupilarea = pupilarea;
        dat.pupilcom = pupilcom;
        dat.videoTimestampsEphysTime_eye = videoTimestampsEphysTime_eye;
        dat.motSVD = motSVD;
        dat.uMotMask = uMotMask;
        dat.videoTimestampsEphysTime_face = videoTimestampsEphysTime_face;
        save2ONE(dataPath,expRef,dat)
    end