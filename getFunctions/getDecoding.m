function acc = getDecoding(db,P,decP)
    
    % looping parameters
    focusOnL = decP.focusOnL;
    
    Pspk = P; % really suboptimal
    Pspk.reductionMethod = 'SVD';
    Pspk.pc2Keep = decP.pc2KeepSpk;

    Pmot = P;
    Pmot.reductionMethod = 'none';
    Pmot.pc2Keep = decP.pc2KeepMot;
    Pmot.zscore = 0; 
    
    Ppup = P;
    Ppup.reductionMethod = 'none';
    Ppup.pc2Keep = {'all'};

    % to save all
    aspk = cell(numel(focusOnL),P.nMice);
    amot = cell(numel(focusOnL),P.nMice);
    apupblink = cell(numel(focusOnL),P.nMice);
    apuparea = cell(numel(focusOnL),P.nMice);
    apupcom = cell(numel(focusOnL),P.nMice);
    apupmot = cell(numel(focusOnL),P.nMice);
    apupall = cell(numel(focusOnL),P.nMice);
    aspkshuf = cell(numel(focusOnL),P.nMice,decP.numshuf);
    
    % to save mean / won't work if not same number of video/sounds
    acc.spk = nan(numel(focusOnL),P.nMice,numel(Pspk.pc2Keep));
    acc.motonly = nan(numel(focusOnL),P.nMice,numel(Pmot.pc2Keep));
    acc.pupblink = nan(numel(focusOnL),P.nMice,1);
    acc.puparea = nan(numel(focusOnL),P.nMice,1);
    acc.pupcom = nan(numel(focusOnL),P.nMice,1);
    acc.pupmot = nan(numel(focusOnL),P.nMice,1);
    acc.pupalleye = nan(numel(focusOnL),P.nMice,1);
    acc.spkshuf = nan(numel(focusOnL),P.nMice,numel(Pspk.pc2Keep),decP.numshuf);
    
    fprintf('*** Running decoding analysis... ***\n') 
    for k = 1:P.nMice
        fprintf('Computing for mouse %s...', P.mouseRef{k})
        for f = 1:numel(focusOnL)
            Pspk.focusOn = focusOnL{f};
            Pmot.focusOn = focusOnL{f};
            Ppup.focusOn = focusOnL{f};
            aspk{f, k} = decodeDataXVal(db(k).spikeData, db(k).audiovisuoCode, Pspk);
            amot{f, k} = decodeDataXVal(db(k).motionData(:,1:200,:,1), db(k).audiovisuoCode, Pmot);
            
            %save means
            acc.spk(f, k, :) = mean(aspk{f, k},2);
            acc.motonly(f, k, :) = mean(amot{f, k},2);
            
            if P.miceweye(k)
                % for pupil / maybe change reduction method?
                apupblink{f, k} = decodeDataXVal(db(k).blinkData, db(k).audiovisuoCode, Ppup);
                apuparea{f, k} = decodeDataXVal(db(k).pupilareaData, db(k).audiovisuoCode, Ppup);
                apupcom{f, k} = decodeDataXVal(db(k).pupilcomData, db(k).audiovisuoCode, Ppup);
                apupmot{f, k} = decodeDataXVal(db(k).pupilmotData, db(k).audiovisuoCode, Ppup);
                apupall{f, k} = decodeDataXVal(cat(2,db(k).blinkData,db(k).pupilareaData,db(k).pupilcomData,db(k).pupilmotData), ...
                    db(k).audiovisuoCode, Ppup);
                
                acc.pupblink(f, k, :) = mean(apupblink{f, k},2);
                acc.puparea(f, k, :) = mean(apuparea{f, k},2);
                acc.pupcom(f, k, :) = mean(apupcom{f, k},2);
                acc.pupmot(f, k, :) = mean(apupmot{f, k},2);
                acc.pupalleye(f, k, :) = mean(apupall{f, k},2);
            else
                acc.pupblink(f, k, :) = nan;
                acc.puparea(f, k, :) = nan;
                acc.pupcom(f, k, :) = nan;
                acc.pupmot(f, k, :) = nan;
                acc.pupalleye(f, k, :) = nan;
            end
            
            if k == 1 % on first mouse only as a test..?
                for ss = 1:decP.numshuf
                    code = [];
                    for rep = 1:4
                        codetmp = dbart(k).audiovisuoCode(:,1:P.nGroupsSnd*P.nGroupsVid); % same over repeats
                        for g = 1:12 % shuffle imp
                            codetmp(strcmp(P.focusOn,{'video','sound'}),codetmp(~strcmp(P.focusOn,{'video','sound'}),:) == P.labelsGroupsSnd(g)) = ...
                                P.labelsGroupsSnd(randperm(12));
                        end
                        for g = 1:12 % shuffle unimp
                            codetmp(~strcmp(P.focusOn,{'video','sound'}),codetmp(strcmp(P.focusOn,{'video','sound'}),:) == P.labelsGroupsSnd(g)) = ...
                                P.labelsGroupsSnd(randperm(12));
                        end
                        code = [code, codetmp];
                    end
                    
                    code = code(:,1:size(db(k).audiovisuoCode,2));
                    
                    aspkshuf{f, k, ss} = decodeDataXVal(db(k).spikeData, code, P);
                    acc.spkshuf(f, k, :, ss) = mean(aspkshuf{f, k, ss},2);
                end
            end
        end
        fprintf('Done.\n')
    end
    fprintf('*** Decoding analysis done. ***\n') 
end