function [ypred_test, W_all] = predictComponent(db,P,predP,recompute)

     %%% This function needs the glmnet-matlab repo to work (https://github.com/junyangq/glmnet-matlab)

    if ~exist('recompute','var')
        recompute = 0;
    end

    P.focusOn = predP.focusOn;
    pc2predL = predP.pc2predL;
    motpc2keepL = predP.motpc2keepL;
    nPCmax = motpc2keepL(end);
    motpc2keepfix = predP.motpc2keepfix; 
    model = predP.model;
    lagbins = predP.lagbins;
    tmpSaveFolder = predP.tmpSaveFolder;
    mkdir_CB(tmpSaveFolder)

    for k = 1:numel(db)
        tmpSaveFile = fullfile(tmpSaveFolder,sprintf('tmpSave_%s_focusOn-%s_reg-%s_mouse-%s-%d.mat',P.focusOnArea,P.focusOn,model.type,P.mouseRef{k},db(k).imecX));
        if ~exist(tmpSaveFile,'file') | recompute
            clear dat
            
            %% Get data
            UAll = getComponents(db(k).spikeData, db(k).audiovisuoCode, P); % get spike data
            dat(1).dataA = getOrgData(UAll,db(k).audiovisuoCode,P); % dataA is bins x imp stim x unimp stim x comp x rep
            si = sign(skewness(mat2vec(nanmean(dat(1).dataA(:,:,:,1,:),[3 5])-nanmean(dat(1).dataA(:,:,:,1,:),[2 3 5])))); % reorient it
            dat(1).dataA = si*dat(1).dataA;
            dat(1).name = 'Neural';
            
            %% Get motion-related predictors
            dat(2).dataA = getOrgData(db(k).motionData,db(k).audiovisuoCode,P); % get motion data
            dat(2).name = 'Motion';
            
            %% Get eye-related predictors
            if P.miceweye(k)
                dattmp = getEyeData(db(k),P);
                for dd = 1:numel(dattmp)
                    dat(2+dd).dataA = dattmp(dd).dataA;
                    dat(2+dd).name = dattmp(dd).name;
                end
                clear dattmp
            end
            
            %% Recentering

            if P.demeanTimeCourse
                fprintf('Remove unimportant stim contribution.\n')
                for dd = 1:numel(dat)
                    dat(dd).dataA = dat(dd).dataA - nanmean(dat(dd).dataA,[2 5]);
                end
            end
            
            %% Shape data & predictors
            % get what to predict
            y = mat2pred(dat(1).dataA,numel(lagbins)); y = y(:,pc2predL);
            
            % motion related predictors
            xmot = mat2pred(dat(2).dataA,numel(lagbins)); xmot = nanzscore(xmot(:,1:nPCmax));
            Xmot = getXpred(xmot,lagbins); clear xmot
            
            % eye related predictors
            if P.miceweye(k)
                xeye = nan(size(y,1),numel(dat)-2);
                for c = 1:numel(dat)-2
                    xeye(:,c) = mat2pred(dat(2+c).dataA,numel(lagbins));
                end
                sig = nanstd(xeye(:,~all(xeye==0)));
                sig(sig==0) = 1;
                xeye(:,~all(xeye==0)) = (xeye(:,~all(xeye==0))-nanmean(xeye(:,~all(xeye==0))))./sig;
                Xeye = getXpred(xeye,lagbins); clear xeye
            else
                Xeye = [];
            end
            
            %% stimulus predictors
            if predP.stim
                xstim = zeros(size(y,1),size(dat(1).dataA,2));
                for st = 1:size(xstim,2)
                    xstimtmp = zeros(size(dat(1).dataA(:,:,:,1,:)));
                    xstimtmp(find(P.bins>0,1),st,:,:,:) = 1;
                    xstim(:,st) = mat2pred(xstimtmp,numel(lagbins));
                end
                xstim(isnan(xstim)) = 0; % here putting 0 should be fine as anyway these intertrials padded timepoints will be removed?
                Xstim = getXpred(xstim,predP.stimlagbins); clear xeye
            end
            
            %% initialize ypred
            s = size(dat(1).dataA);
            s(1) = s(1) + 2*numel(lagbins);
            for pc = 1:numel(pc2predL)
                ypred_test_eye_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                if predP.stim
                    ypred_test_stim_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                end
                ypred_test_mot_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                ypred_test_all_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                if predP.stim
                    ypred_test_allwstim_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                    ypred_test_allwstim_motonly_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                    ypred_test_allwstim_eyeonly_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                    ypred_test_allwstim_stimonly_tmp{pc} = nan([s(1),s(2),s(3),1,s(5)]);
                end
            end
            
            %% predict
            [nfold,folds] = getFolds(s,P.foldtype);
            
            W = nan(motpc2keepfix*numel(lagbins)+1, nfold, numel(pc2predL));
            
            for ff = 1:nfold
                disp(['Mouse ' P.mouseRef{k} ' / fold #' num2str(ff)])
                idxfold = folds{ff};
                idxfold = idxfold(:,:,:,1,:);
                
                idxfold = mat2pred(idxfold);
                idxfold(isnan(idxfold(:))) = 0;
                % predict
                for p2k = 1:numel(motpc2keepL)
                    % from motion pcs
                    motpc2keep = motpc2keepL(p2k);
                    [w,~,~,ytest_pred] = ...
                        makePred_paper(y,Xmot(:,1:motpc2keep*numel(lagbins)),idxfold,model);
                    
                    % save
                    if motpc2keep == motpc2keepfix
                        for pc = 1:numel(pc2predL)
                            ypred_test_mot_tmp{pc}(idxfold) = ytest_pred(:,pc); % for plt
                        end
                        W(:,ff,:,:) = w;
                    end
                end
                
                if P.miceweye(k)
                    
                    modeleye = model;
                    if strcmp(predP.model.type,'rrr') %%% not the best way to put it...
                        modeleye.type = 'ols';
                    end
                    % all eye pred together
                    [~,~,~,ytest_pred] = makePred_paper(y,Xeye,idxfold,modeleye);
                    for pc = 1:numel(pc2predL)
                        ypred_test_eye_tmp{pc}(idxfold) = ytest_pred(:,pc); % for plt
                    end
                else
                    for pc = 1:numel(pc2predL)
                        ypred_test_eye_tmp{pc}(:) = 1;
                    end
                end
                
                % altogether (with fixed number of motion PC...)
                X = [Xmot(:,1:motpc2keepfix*numel(lagbins)),Xeye];
                [~,~,~,ytest_pred] = ...
                    makePred_paper(y,X,idxfold,model);
                for pc = 1:numel(pc2predL)
                    ypred_test_all_tmp{pc}(idxfold) = ytest_pred(:,pc);
                end
                
                if predP.stim
                    X = Xstim;
                    modelstim = model;
                    [~,~,~,ytest_pred] = makePred_paper(y,X,idxfold,modelstim);
                    for pc = 1:numel(pc2predL)
                        ypred_test_stim_tmp{pc}(idxfold) = ytest_pred(:,pc); % for plt
                    end
                    
                    X = [Xmot(:,1:motpc2keepfix*numel(lagbins)),Xeye,Xstim];
                    [w,~,~,ytest_pred] = makePred_paper(y,X,idxfold,model);
                    for pc = 1:numel(pc2predL)
                        ypred_test_allwstim_tmp{pc}(idxfold) = ytest_pred(:,pc);
                        
                        % save individually for each predictor type
                        % motion
                        ypred_test_allwstim_motonly_tmp{pc}(idxfold) = Xmot(idxfold,1:motpc2keepfix*numel(lagbins))*w(1:motpc2keepfix*numel(lagbins),pc);
                        
                        % eye
                        if P.miceweye(k)
                            ypred_test_allwstim_eyeonly_tmp{pc}(idxfold) = Xeye(idxfold,:)*w(motpc2keepfix*numel(lagbins)+(1:size(Xeye,2)),pc);
                        else
                            ypred_test_allwstim_eyeonly_tmp{pc}(idxfold) = nan;
                        end
                        
                        % stimulus
                        ypred_test_allwstim_stimonly_tmp{pc}(idxfold) = Xstim(idxfold,:)*w(motpc2keepfix*numel(lagbins)+size(Xeye,2)+1:end-1,pc);
                    end
                end
            end
            
            W_all = reshape(nanmean(W(1:end-1,:,:,:),2),[numel(lagbins),motpc2keepfix,numel(pc2predL)]);
            s(1) = size(dat(1).dataA,1);

            clear ypred_test % if save it along the way crashes (memory issue?)
            ypred_test.real = dat(1).dataA(:,:,:,pc2predL,:);
            ypred_test.mot = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
            ypred_test.eye = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
            ypred_test.all = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
            if predP.stim
                ypred_test.stim = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
                ypred_test.allwstim = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
                ypred_test.allwstim_motonly = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
                ypred_test.allwstim_eyeonly = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
                ypred_test.allwstim_stimonly = nan([s(1),s(2),s(3),numel(pc2predL),s(5)]);
            end
            for pc = 1:numel(pc2predL)
                ypred_test.eye(:,:,:,pc,:) = ypred_test_eye_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                if predP.stim
                    ypred_test.stim(:,:,:,pc,:) = ypred_test_stim_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                end
                ypred_test.mot(:,:,:,pc,:) = ypred_test_mot_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                ypred_test.all(:,:,:,pc,:) = ypred_test_all_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                if predP.stim
                    ypred_test.allwstim(:,:,:,pc,:) = ypred_test_allwstim_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                    ypred_test.allwstim_motonly(:,:,:,pc,:) = ypred_test_allwstim_motonly_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                    ypred_test.allwstim_eyeonly(:,:,:,pc,:) = ypred_test_allwstim_eyeonly_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                    ypred_test.allwstim_stimonly(:,:,:,pc,:) = ypred_test_allwstim_stimonly_tmp{pc}(numel(lagbins)+1:end-numel(lagbins),:,:,:,:);
                end
            end
            
            save(tmpSaveFile,'ypred_test','W_all')
            % not ideal
            clear ypred_test W_all ypred_test_mot_tmp ypred_test_eye_tmp ypred_test_all_tmp
        end
    end
    
    % reload all
    clear W_all ypred_test
    for k = 1:numel(db)
        tmpSaveFile = fullfile(tmpSaveFolder,sprintf('tmpSave_%s_focusOn-%s_reg-%s_mouse-%s-%d.mat',P.focusOnArea,P.focusOn,model.type,P.mouseRef{k},db(k).imecX));
        
        tmp = load(tmpSaveFile);
        W_all{k} = tmp.W_all;
        fnames = fieldnames(tmp.ypred_test);
        for f = 1:numel(fnames)
            ypred_test(k).(fnames{f}) = tmp.ypred_test.(fnames{f});
        end
    end
end