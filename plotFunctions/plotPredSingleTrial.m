function [ev_half,z,p] = plotPredSingleTrial(predDat,mainFigP)

    %% Plot single trials

    % Parameters -- should be inputs!!
    k2plt = 1;
    trials2plt(1,:) = [4 5 6]; % videos
    trials2plt(2,:) = [1 1 1]; %  repeats

    % what stimuli
    idxstimsubset = [];
    for idx = 1:numel(mainFigP.stimsubset)
        idxstimsubset = [idxstimsubset, find(mainFigP.labelsGroupsSnd == mainFigP.stimsubset(idx))];
    end
    nstim = numel(idxstimsubset);

    ypred_test = predDat.ypred_test;

    %%
    figure('Position', [680   200   436   240]);

    % plot sound waveforms
    sndPath = '\\zserver\Code\Rigging\ExpDefinitions\Tim\filmWorldwGaps\30-2019-11-22-ratio-edited\processed-merged-4-sec-w-corner-150';
    for sidx = 1:nstim
        subplot(size(trials2plt,2)+1,nstim,sidx); hold all
        [audio, Fs] = audioread(fullfile(sndPath, ['video-0-audio-' num2str(mainFigP.labelsGroupsSnd(idxstimsubset(sidx))) '.mp4']));
        snd = sum(audio,2);
        plot(1/Fs:1/Fs:(size(audio,1)/Fs+0.135),cat(1,zeros(0.135*Fs,1),snd),'color',[.5 .5 .5],'LineWidth',1)
        xlim([mainFigP.bins(1),mainFigP.bins(end)])
        ylim([-2,2])
        ylabel('t') % to make sure all has the same size
        set(gca,'visible','off')
    end

    for tt = 1:size(trials2plt,2)

        for sidx = 1:numel(idxstimsubset)
            s = idxstimsubset(sidx);
            subplot(size(trials2plt,2)+1,nstim,nstim*(tt)+sidx); hold all

            clr = mainFigP.SndColors(s,:);
            plot(mainFigP.bins,ypred_test(k2plt).real(:,s,trials2plt(1,tt),1,trials2plt(2,tt)),'color',[128 0 128]/255,'LineWidth',1)

            plot(mainFigP.bins,ypred_test(k2plt).all(:,s,trials2plt(1,tt),1,trials2plt(2,tt),end),'color',[0 0 1],'LineWidth',1)
            plot(mainFigP.bins,ypred_test(k2plt).stim(:,s,trials2plt(1,tt),1,trials2plt(2,tt)),'color',[.5 .5 .5],'LineWidth',1)

            ylim([min(mat2vec(ypred_test(k2plt).real(:,s,trials2plt(1,:),1,trials2plt(2,:)))), ...
                max(mat2vec(ypred_test(k2plt).real(:,s,trials2plt(1,:),1,trials2plt(2,:))))])
            vline(0,'k')
            set(gca,'visible','off')
        end
    end
    set(get(gcf,'children'),'clipping','off')

    %% "noise" correlation
    % supposes that used oddeven fold type

    fnames = fieldnames(ypred_test(1));
    getQuant = @getCorr;
    for k = 1:mainFigP.nMice
        s = size(ypred_test(k).real);

        [nfold,folds] = getFolds(s,'oddeven');

        for ff = 1:nfold
            testfolds = find(squeeze(folds{ff}(1,1,1,1,:)));
            trainfolds = find(~squeeze(folds{ff}(1,1,1,1,:)));

            %         testreal = reshape(nanmean(ypred_test(k).real(:,:,:,:,testfolds),[3 5]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2),size(ypred_test(k).real,4)]);
            testreal = reshape(permute(ypred_test(k).real(:,:,:,:,testfolds)-nanmean(ypred_test(k).real(:,:,:,:,testfolds),[3 5]),[1 2 3 5 4]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2)*size(ypred_test(k).real,3)*numel(testfolds),size(ypred_test(k).real,4)]);

            for dd = 2:numel(fnames)
                %             testf = reshape(nanmean(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds,:),[3 5]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2),size(ypred_test(k).real,4),size(ypred_test(k).(fnames{dd}),6)]);
                testf = reshape(permute(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds,:)-nanmean(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds),[3 5]),[1 2 3 5 4 6]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2)*size(ypred_test(k).real,3)*numel(testfolds),size(ypred_test(k).real,4),size(ypred_test(k).(fnames{dd}),6)]);

                for rrrpc = 1:size(ypred_test(k).(fnames{dd}),6)
                    % single PC variance explained
                    ev_sing_half.(fnames{dd})(k,:,rrrpc,ff) = getQuant(testreal,testf(:,:,rrrpc));
                    % whole subspace variance explained
                    ev_half.(fnames{dd})(k,rrrpc,ff) = getQuant(testreal(:),mat2vec(testf(:,:,rrrpc)));
                end
            end

            if all(ypred_test(k).eye(:) == 1)
                ev_half.all(k,:,:) = ev_half.mot(k,:,:);
                ev_sing_half.all(k,:,:,:) = ev_sing_half.mot(k,:,:,:);
            end
        end
        c{k}(1,:) = testreal(:,1);
        testf = reshape(permute(ypred_test(k).all(:,:,:,:,testfolds,:)-nanmean(ypred_test(k).(fnames{dd})(:,:,:,:,testfolds),[3 5]),[1 2 3 5 4 6]),[size(ypred_test(k).real,1)*size(ypred_test(k).real,2)*size(ypred_test(k).real,3)*numel(testfolds),size(ypred_test(k).real,4),size(ypred_test(k).(fnames{dd}),6)]);
        c{k}(2,:) = testf(:,1);
    end

    figure('Position', [696   885   378    98]);
    subplot(121)
    Ptmp = mainFigP;
    Ptmp.nMice = 1;
    xedges = linspace(prctile(c{1}(1,:),1),prctile(c{1}(1,:),99),100);
    yedges = xedges;
    plotScatterHist(c(1),Ptmp,xedges,yedges,1)
    offsetAxes
    caxis([1 30])
    xlabel('Activity - noise')
    ylabel('Prediction - noise')
    % scatter(c{1}(1,:),c{1}(2,:),30,P.mouseColor(1,:))
    axis equal tight
    subplot(122)
    hold all
    z = ztransform(nanmean(ev_half.all,3),'mean');
    scatter(rand(1,size(ev_half.all,1))*0.5-0.2,nanmean(ev_half.all,3),20,'k')
    scatter(0,z,30,'k','filled')
    ylim([0 1])
    p = signrank(nanmean(ev_half.all,3));
    sigstar({[.9, 1.1]},p)