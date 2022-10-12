function [quant,p,quant_sing] = plotPredResults(pred,sideTransec)
    
    if ~exist('sideTransec','var')
        sideTransec = -ones(1,numel(pred.ypred_test));
    end

    ypred_test = pred.ypred_test;
    predP = pred.predP;
    getQuant = predP.getQuant;
    
    %%
    
    fnames = fieldnames(ypred_test);
    for k = 1:numel(ypred_test)
        idxnan = isnan(ypred_test(k).allwstim(:,:,:,:,:,1));
        for f = 1:numel(fnames)
            ypred_test(k).(fnames{f})(repmat(idxnan,[1 1 1 1 1 size(ypred_test(k).(fnames{f}),6)])) = nan;
        end
    end
    
    clear quant_sing quant
    for k = 1:numel(ypred_test)
        s = size(ypred_test(k).real);
        
        [nfold,folds] = getFolds(s,predP.foldtype);
        
        for ff = 1:nfold
            testfolds = find(squeeze(folds{ff}(1,1,1,1,:)));
            trainfolds = find(~squeeze(folds{ff}(1,1,1,1,:)));
            
            testreal = reshape(nanmean(ypred_test(k).real(:,:,:,:,testfolds),[3 5]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2),size(ypred_test(k).real,4)]);
            %         testreal = reshape(permute(ypred_test(k).real(:,:,:,:,testfolds),[1 2 3 5 4]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2)*size(ypred_test(k).real,3)*numel(testfolds),size(ypred_test(k).real,4)]);
            
            for dd = 2:numel(fnames)
                testf = reshape(nanmean(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds,:),[3 5]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2),size(ypred_test(k).real,4),size(ypred_test(k).(fnames{dd}),6)]);
                %             testf = reshape(permute(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds,:),[1 2 3 5 4 6]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2)*size(ypred_test(k).real,3)*numel(testfolds),size(ypred_test(k).real,4),size(ypred_test(k).(fnames{dd}),6)]);
                
                for rrrpc = 1:size(ypred_test(k).(fnames{dd}),6)
                    % single PC variance explained
                    quant_sing.(fnames{dd})(k,:,rrrpc,ff) = getQuant(testreal,testf(:,:,rrrpc));
                    % whole subspace variance explained
                    quant.(fnames{dd})(k,rrrpc,ff) = getQuant(mat2vec(testreal-nanmean(testreal)), ...
                        mat2vec(testf(:,:,rrrpc)-nanmean(testf(:,:,rrrpc))));
                end
            end
            
            if all(ypred_test(k).eye(:) == 1)
                quant.all(k,:,:) = quant.mot(k,:,:);
                quant_sing.all(k,:,:,:) = quant_sing.mot(k,:,:,:);
            end
        end
    end
    
    % Main figure
    figure('Position', [ 696   726   390   257]);
    subplot(131); hold all
    scatter(nanmean(quant.stim(sideTransec == -1,:,:),3),nanmean(quant.all(sideTransec == -1,:,:),3),20,'k')
    scatter(nanmean(quant.stim(sideTransec == 1,:,:),3),nanmean(quant.all(sideTransec == 1,:,:),3),20,'k','x')
    m =  1;% 1.1*max([nanmean(ev_half.mot(:,end,:),3);nanmean(ev_half.allwstim(:,end,:),3)]);
    plot([0 m],[0 m],'k--')
    axis equal tight
    xlim([0 m])
    ylabel('Behavioral')
    xlabel('Sensory')
    p(1) = signrank(nanmean(quant.stim,3),nanmean(quant.all,3));
    title(sprintf('p = %d',p(1)));
    offsetAxes
    
    subplot(132); hold all
    scatter(nanmean(quant.stim(sideTransec == -1,:,:),3),nanmean(quant.allwstim(sideTransec == -1,:,:),3),20,'k')
    scatter(nanmean(quant.stim(sideTransec == 1,:,:),3),nanmean(quant.allwstim(sideTransec == 1,:,:),3),20,'k','x')
    m =  1;% 1.1*max([nanmean(ev_half.eye(:,end,:),3);nanmean(ev_half.allwstim(:,end,:),3)]);
    plot([0 m],[0 m],'k--')
    axis equal tight
    xlim([0 m])
    ylabel('Full (Beh. + Sens.)')
    xlabel('Sensory')
    p(2) = signrank(nanmean(quant.stim,3),nanmean(quant.allwstim,3));
    title(sprintf('p = %d',p(2)));
    offsetAxes

    subplot(133); hold all
    scatter(nanmean(quant.all(sideTransec == -1,:,:),3),nanmean(quant.allwstim(sideTransec == -1,:,:),3),20,'k')
    scatter(nanmean(quant.all(sideTransec == 1,:,:),3),nanmean(quant.allwstim(sideTransec == 1,:,:),3),20,'k','x')
    m =  1;% 1.1*max([nanmean(ev_half.mot(:,end,:),3);nanmean(ev_half.allwstim(:,end,:),3)]);
    plot([0 m],[0 m],'k--')
    axis equal tight
    xlim([0 m])
    ylabel('Full (Beh. + Sens.)')
    xlabel('Behavioral')
    p(3) = signrank(nanmean(quant.all,3),nanmean(quant.allwstim,3));
    title(sprintf('p = %d',p(3)));
    offsetAxes

    % to display the stuff that's around
    set(get(gcf,'children'),'clipping','off')

    % Eye specific figure
    miceweye = all(quant.eye>-inf,3);
    if any(miceweye)
        figure('Position', [696   752   266   231]);
        subplot(121); hold all
        scatter(nanmean(quant.all(miceweye,:,:),3),nanmean(quant.mot(miceweye,:,:),3),20,'k')
        m =  1;% 1.1*max([nanmean(ev_half.stim(:,end,:),3);nanmean(ev_half.allwstim(:,end,:),3)]);
        plot([0 m],[0 m],'k--')
        axis equal tight
        xlim([0 m])
        ylabel('Body only')
        xlabel('Eye and body')
        p(4) = signrank(nanmean(quant.all,3),nanmean(quant.mot,3));
        title(sprintf('p = %d',p(4)));
        offsetAxes

        subplot(122); hold all
        scatter(nanmean(quant.all(miceweye,:,:),3),nanmean(quant.eye(miceweye,:,:),3),20,'k')
        m =  1;% 1.1*max([nanmean(ev_half.stim(:,end,:),3);nanmean(ev_half.allwstim(:,end,:),3)]);
        plot([0 m],[0 m],'k--')
        axis equal tight
        xlim([0 m])
        ylabel('Eye only')
        xlabel('Eye and body')
        p(5) = signrank(nanmean(quant.all,3),nanmean(quant.eye,3));
        title(sprintf('p = %d',p(5)));
        offsetAxes
    end
end