function plotCompDetails(spkVIS,spkHPF,spkPredVIS,spkPredHPF,focusOn)
    %%% Look at similarity of responses across PCs/regions, etc.

    %% Process a bit the data

    mainFigP = spkVIS.mainFigP;

    nMiceVIS = size(spkVIS.dataMPC_even,4);
    nMiceHPF = size(spkHPF.dataMPC_even,4);

    % reorient the pcs
    s = size(spkVIS.dataMPC);
    for pc = 1:s(3)
        % reorient VIS
        si = sign(corr(reshape(spkVIS.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceVIS]), ...
            nanmean(reshape(spkVIS.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceVIS]),2)));
        spkVIS.dataMPC(:,:,pc,:) = spkVIS.dataMPC(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);
        spkVIS.dataMPC_odd(:,:,pc,:) = spkVIS.dataMPC_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);
        spkVIS.dataMPC_even(:,:,pc,:) = spkVIS.dataMPC_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);

        % reorient HPF
        si = sign(corr(reshape(spkHPF.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceHPF]), ...
            nanmean(reshape(spkVIS.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceVIS]),2)));
        spkHPF.dataMPC(:,:,pc,:) = spkHPF.dataMPC(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
        spkHPF.dataMPC_odd(:,:,pc,:) = spkHPF.dataMPC_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
        spkHPF.dataMPC_even(:,:,pc,:) = spkHPF.dataMPC_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
    end

    P = mainFigP;

    %% Load pred

    model2keep = 'all';

    % VIS
    for k = 1:nMiceVIS
        spkPredVIS.dataMPC_real_odd(:,:,:,k) = squeeze(nanmean(spkPredVIS.ypred_test(k).real(:,:,:,mainFigP.pc2Keep,1:2:end),[3 5]));
        spkPredVIS.dataMPC_pred_odd(:,:,:,k) = squeeze(nanmean(spkPredVIS.ypred_test(k).(model2keep)(:,:,:,mainFigP.pc2Keep,1:2:end),[3 5]));

        spkPredVIS.dataMPC_real_even(:,:,:,k) = squeeze(nanmean(spkPredVIS.ypred_test(k).real(:,:,:,mainFigP.pc2Keep,2:2:end),[3 5]));
        spkPredVIS.dataMPC_pred_even(:,:,:,k) = squeeze(nanmean(spkPredVIS.ypred_test(k).(model2keep)(:,:,:,mainFigP.pc2Keep,2:2:end),[3 5]));
    end

    % HPF
    for k = 1:nMiceHPF
        spkPredHPF.dataMPC_real_odd(:,:,:,k) = squeeze(nanmean(spkPredHPF.ypred_test(k).real(:,:,:,mainFigP.pc2Keep,1:2:end),[3 5]));
        spkPredHPF.dataMPC_pred_odd(:,:,:,k) = squeeze(nanmean(spkPredHPF.ypred_test(k).(model2keep)(:,:,:,mainFigP.pc2Keep,1:2:end),[3 5]));

        spkPredHPF.dataMPC_real_even(:,:,:,k) = squeeze(nanmean(spkPredHPF.ypred_test(k).real(:,:,:,mainFigP.pc2Keep,2:2:end),[3 5]));
        spkPredHPF.dataMPC_pred_even(:,:,:,k) = squeeze(nanmean(spkPredHPF.ypred_test(k).(model2keep)(:,:,:,mainFigP.pc2Keep,2:2:end),[3 5]));
    end

    % reorient the pcs
    for pc = 1:size(spkVIS.dataMPC,3)
        % reorient VIS
        si = sign(corr(reshape(spkPredVIS.dataMPC_real_odd(:,:,pc,:),[s(1)*s(2),nMiceVIS]), ...
            nanmean(reshape(spkVIS.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceVIS]),2)));
        spkPredVIS.dataMPC_real_odd(:,:,pc,:) = spkPredVIS.dataMPC_real_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);
        spkPredVIS.dataMPC_pred_odd(:,:,pc,:) = spkPredVIS.dataMPC_pred_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);
        spkPredVIS.dataMPC_real_even(:,:,pc,:) = spkPredVIS.dataMPC_real_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);
        spkPredVIS.dataMPC_pred_even(:,:,pc,:) = spkPredVIS.dataMPC_pred_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceVIS]);

        % reorient HPF
        si = sign(corr(reshape(spkPredHPF.dataMPC_real_odd(:,:,pc,:),[s(1)*s(2),nMiceHPF]), ...
            nanmean(reshape(spkVIS.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceVIS]),2)));
        spkPredHPF.dataMPC_real_odd(:,:,pc,:) = spkPredHPF.dataMPC_real_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
        spkPredHPF.dataMPC_pred_odd(:,:,pc,:) = spkPredHPF.dataMPC_pred_odd(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
        spkPredHPF.dataMPC_real_even(:,:,pc,:) = spkPredHPF.dataMPC_real_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
        spkPredHPF.dataMPC_pred_even(:,:,pc,:) = spkPredHPF.dataMPC_pred_even(:,:,pc,:).*reshape(si,[1 1 1 nMiceHPF]);
    end

    stimsubset = P.labelsGroupsSnd;

    %% Get summary correlations

    % Compute it
    s = size(spkVIS.dataMPC_even);
    clear reliability_VIS similarity_VIS_acrossAnimals reliability_HPF similarity_HPF_acrossAnimals ...
        similarity_VISHPF_acrossRegions_withinAnimals similarity_VISHPF_acrossRegions_acrossAnimals ...
        predictability_VIS predictability_HPF similarity_VISHPF_acrossRegions

    for pc = 1:s(3)
        % Reliability + similarity across animals VIS
        c = corr(reshape(spkVIS.dataMPC_even(:,:,pc,:),[s(1)*s(2),nMiceVIS]), ...
            reshape(spkVIS.dataMPC_odd(:,:,pc,:),[s(1)*s(2),nMiceVIS]));
        c = (c+c')/2;
        reliability_VIS(:,pc) = diag(c);
        similarity_VIS_acrossAnimals(:,pc) = mat2vec(c(triu(true(size(c)),1)));

        % Reliability + similarity across animals HPF
        c = corr(reshape(spkHPF.dataMPC_even(:,:,pc,:),[s(1)*s(2),nMiceHPF]), ...
            reshape(spkHPF.dataMPC_odd(:,:,pc,:),[s(1)*s(2),nMiceHPF]));
        c = (c+c')/2;
        reliability_HPF(:,pc) = diag(c);
        similarity_HPF_acrossAnimals(:,pc) = mat2vec(c(triu(true(size(c)),1)));

        % Similarity across regions
        orderVIS = [find(ismember(spkVIS.mainFigP.mouseRef,spkHPF.mainFigP.mouseRef)) ...
            find(~ismember(spkVIS.mainFigP.mouseRef,spkHPF.mainFigP.mouseRef))];
        c = corr(reshape(spkVIS.dataMPC(:,:,pc,orderVIS),[s(1)*s(2),nMiceVIS]), ...
            reshape(spkHPF.dataMPC(:,:,pc,:),[s(1)*s(2),nMiceHPF]));
        ctrunc = c(1:nMiceHPF,1:nMiceHPF);
        ctrunc = (ctrunc+ctrunc')/2;
        similarity_VISHPF_acrossRegions(:,pc) = c(:);
        similarity_VISHPF_acrossRegions_withinAnimals(:,pc) = diag(ctrunc);
        similarity_VISHPF_acrossRegions_acrossAnimals(:,pc) = [mat2vec(ctrunc(triu(true(size(ctrunc)),1))); ...
            mat2vec(c(nMiceHPF+1:end,:))];

        % Predictability VIS
        x1 = reshape(spkPredVIS.dataMPC_real_odd(:,:,pc,:),[s(1)*s(2),nMiceVIS]);
        x2 = reshape(spkPredVIS.dataMPC_pred_odd(:,:,pc,:),[s(1)*s(2),nMiceVIS]);
        idxnan = any(isnan(x2'));
        c_odd = corr(x1(~idxnan,:),x2(~idxnan,:));
        x1 = reshape(spkPredVIS.dataMPC_real_even(:,:,pc,:),[s(1)*s(2),nMiceVIS]);
        x2 = reshape(spkPredVIS.dataMPC_pred_even(:,:,pc,:),[s(1)*s(2),nMiceVIS]);
        idxnan = any(isnan(x2'));
        c_even = corr(x1(~idxnan,:),x2(~idxnan,:));
        c = (c_odd+c_even)/2;
        predictability_VIS(:,pc) = diag(c);

        % Predictability HPF
        x1 = reshape(spkPredHPF.dataMPC_real_odd(:,:,pc,:),[s(1)*s(2),nMiceHPF]);
        x2 = reshape(spkPredHPF.dataMPC_pred_odd(:,:,pc,:),[s(1)*s(2),nMiceHPF]);
        idxnan = any(isnan(x2'));
        c_odd = corr(x1(~idxnan,:),x2(~idxnan,:));
        x1 = reshape(spkPredHPF.dataMPC_real_even(:,:,pc,:),[s(1)*s(2),nMiceHPF]);
        x2 = reshape(spkPredHPF.dataMPC_pred_even(:,:,pc,:),[s(1)*s(2),nMiceHPF]);
        idxnan = any(isnan(x2'));
        c_even = corr(x1(~idxnan,:),x2(~idxnan,:));
        c = (c_odd+c_even)/2;
        predictability_HPF(:,pc) = diag(c);
    end

    %% Plot actual predictions by movements

    pc2plt = 1;

    mainFigP.plt.alph = .4;
    wpatch = 0;
    mainFocusOn = mainFigP.focusOn;
    switch mainFocusOn
        case 'sound'
            mainCol = ones(size(mainFigP.SndColors)).*nanmean(mainFigP.SndColors);
            secCol =  ones(size(mainFigP.VidColors)).*nanmean(mainFigP.VidColors);
        case 'video'
            mainCol = ones(size(mainFigP.VidColors)).*nanmean(mainFigP.VidColors);
            secCol =  ones(size(mainFigP.SndColors)).*nanmean(mainFigP.SndColors);
    end

    idxstimsubset = [];
    for idx = 1:numel(stimsubset)
        idxstimsubset = [idxstimsubset, find(mainFigP.labelsGroupsSnd == stimsubset(idx))];
    end
    nstim = numel(idxstimsubset);

    % Plot it
    figure('Position',[230 591 762 219])
    numcol_splt = nstim;
    for c = 1:2*numel(pc2plt)+1
        for sidx = 1:numcol_splt
            subplt(c,sidx) = subplot(1+2*numel(pc2plt),numcol_splt,(c-1)*numcol_splt+sidx);
        end
    end

    % Plot stim
    sndPath = 'C:\Users\Hamish\OneDrive - University College London\Dossiers_UCL\Papers\2022-Bimbard_etal-Behavioural-origin-of-sound-evoked-activity-in-mouse-visual-cortex\Stimuli';
    switch focusOn
        case 'sound'
            for sidx = 1:nstim
                set(gcf,'CurrentAxes', subplt(1,sidx));
                [audio, Fs] = audioread(fullfile(sndPath, ['video-0-audio-' num2str(mainFigP.labelsGroupsSnd(idxstimsubset(sidx))) '.mp4']));
                snd = sum(audio,2);
                plot(1/Fs:1/Fs:(size(audio,1)/Fs+0.135),cat(1,zeros(0.135*Fs,1),snd),'color',[.5 .5 .5],'LineWidth',1)
                xlim([mainFigP.bins(1),mainFigP.bins(end)])
                ylim([-2,2])
                ylabel('t') % to make sure all has the same size
                set(gca,'visible','off')
            end
        case 'video'
            for sidx = 1:nstim
                set(gcf,'CurrentAxes', subplt(1,sidx));
                im = imread(fullfile(sndPath, sprintf('video%d.PNG',mainFigP.labelsGroupsSnd(idxstimsubset(sidx)))));
                imagesc(im)
                set(gca,'visible','off')
                axis equal tight
            end
    end

    % plot responses
    for pc = 1:numel(pc2plt)
        plotPSTH_fig1(squeeze(nanmean(spkPredVIS.dataMPC_real_odd(:,:,pc2plt(pc),:),4)),idxstimsubset,mainFigP,subplt(2+(pc-1)*2,1:nstim),wpatch,0,mainCol)
        plotPSTH_fig1(squeeze(nanmean(spkPredVIS.dataMPC_pred_odd(:,:,pc2plt(pc),:),4)),idxstimsubset,mainFigP,subplt(2+(pc-1)*2,1:nstim),wpatch,0,mainCol,'--')

        plotPSTH_fig1(squeeze(nanmean(spkPredHPF.dataMPC_real_odd(:,:,pc2plt(pc),:),4)),idxstimsubset,mainFigP,subplt(3+(pc-1)*2,1:nstim),wpatch,0,mainCol)
        plotPSTH_fig1(squeeze(nanmean(spkPredHPF.dataMPC_pred_odd(:,:,pc2plt(pc),:),4)),idxstimsubset,mainFigP,subplt(3+(pc-1)*2,1:nstim),wpatch,0,mainCol,'--')
    end


    %% Plot correlations

    nPC = size(reliability_VIS,2);

    %figure('Name',focusOn,'Position',[600   602   444   328]);
    figure('Name',focusOn,'Position',[600   602   557   328]);

    % Reliability
    % VIS
    subplot(241); hold all
    b = boxchart(reliability_VIS);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    h = findobj(gca,'Tag','Mean'); set(h,'Visible','off');
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(reliability_VIS,1))-0.2,reliability_VIS(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(reliability_VIS,'mean'),40,'k','filled');
    ylim([-0.2 1])
    ylabel({'Reliability'})
    hline(0)
    title('VIS')

    % HPF
    subplot(245); hold all
    b = boxchart(reliability_HPF);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(reliability_HPF,1))-0.2,reliability_HPF(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(reliability_HPF,'mean'),40,'k','filled');
    ylim([-0.2 1])
    ylabel({'Reliability'})
    hline(0)
    title('HPF')

    % Similarity across animals
    % VIS
    subplot(242); hold all
    b = boxchart(similarity_VIS_acrossAnimals);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(similarity_VIS_acrossAnimals,1))-0.2,similarity_VIS_acrossAnimals(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(similarity_VIS_acrossAnimals,'mean'),40,'k','filled');
    scatter(1:nPC,ztransform(reliability_VIS,'mean'),40,[.5 .5 .5],'filled','MarkerFaceAlpha', 0.5);
    ylim([-0.2 1])
    ylabel({'Similarity across animals'})
    hline(0)
    title('VIS')

    % HPF
    subplot(246); hold all
    b = boxchart(similarity_HPF_acrossAnimals);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(similarity_HPF_acrossAnimals,1))-0.2,similarity_HPF_acrossAnimals(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(similarity_HPF_acrossAnimals,'mean'),40,'k','filled');
    scatter(1:nPC,ztransform(reliability_HPF,'mean'),40,[.5 .5 .5],'filled','MarkerFaceAlpha', 0.5);
    ylim([-0.2 1])
    ylabel({'Similarity across animals'})
    hline(0)
    title('HPF')

    % Similarity across regions
    subplot(243); hold all
    b = boxchart(similarity_VISHPF_acrossRegions);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(similarity_VISHPF_acrossRegions,1))-0.2,similarity_VISHPF_acrossRegions(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(similarity_VISHPF_acrossRegions,'mean'),40,'r','filled');
    scatter(1:nPC+0.2,ztransform(similarity_VISHPF_acrossRegions,'mean'),40,'k','filled');
    ylim([-0.2 1])
    ylabel({'Similarity across regions'})
    hline(0)
    title('VIS x HPF')

    % Predictability
    % VIS
    subplot(244); hold all
    b = boxchart(predictability_VIS);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(predictability_VIS,1))-0.2,predictability_VIS(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(predictability_VIS,'mean'),40,'k','filled');
    scatter(1:nPC,ztransform(reliability_VIS,'mean'),40,[.5 .5 .5],'filled','MarkerFaceAlpha', 0.5);
    ylim([-0.2 1])
    ylabel({'Predictability'})
    hline(0)
    title('VIS')

    % HPF
    subplot(248); hold all
    b = boxchart(predictability_HPF);
    b.MarkerStyle = '.'; b.MarkerColor = 'k'; b.BoxFaceColor = 'k'; hold all
    for pc = 1:nPC
        scatter(pc+0.4*rand(1,size(predictability_HPF,1))-0.2,predictability_HPF(:,pc),10,[0.5 0.5 0.5]);
    end
    scatter(1:nPC,ztransform(predictability_HPF,'mean'),40,'k','filled');
    scatter(1:nPC,ztransform(reliability_HPF,'mean'),40,[.5 .5 .5],'filled','MarkerFaceAlpha', 0.5);
    ylim([-0.2 1])
    ylabel({'Predictability'})
    hline(0)
    title('HPF')


end