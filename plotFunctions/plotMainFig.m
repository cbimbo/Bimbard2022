function [corrval,pval] = plotMainFig(mainFigDat,eigDat)
    
    %% Rename everything

    mainFigP = mainFigDat.mainFigP;
    dataMPC = mainFigDat.dataMPC;
    M = mainFigDat.M; 
    dataA = mainFigDat.dataA;
    k2plt = mainFigP.k2plt;
    ii2plt = mainFigP.ii2plt;
    r = eigDat.r;
    compDimP = eigDat.compDimP;
    
    %% Parameters

    % plot parameters
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
    
    % what stimuli
    idxstimsubset = [];
    for idx = 1:numel(mainFigP.stimsubset)
        idxstimsubset = [idxstimsubset, find(mainFigP.labelsGroupsSnd == mainFigP.stimsubset(idx))];
    end
    nstim = numel(idxstimsubset);
    
    dataM = nanmean(dataA,5);
    
    %% 1st part of the figure: from single-cell activity to across mice responses
    %% Plot  example cells
    
    numrow_splt = 1 + numel(idxstimsubset) + numel(ii2plt) + numel(k2plt) + 1 + 1 + 1; % sounds + all stim + relative average + other neurons + allneurons imagesc + PC1 + PC1 other mice + summary + fig1 summary
    numcol_splt = nstim + 2; % all stim + delta video average/spectrum

    figure('Position',[549    32   800   964])
    % all subplots
    for c = 1:numrow_splt
        for sidx = 1:numcol_splt
            if (c == 1 && sidx < numcol_splt-1) || ...
                    (ismember(c,2:(1+numel(idxstimsubset)+numel(ii2plt))) && sidx < numcol_splt) || ...
                    (c == (2+numel(idxstimsubset)+numel(ii2plt)) && sidx < numcol_splt-1) || ...
                    ismember(c,(3+numel(idxstimsubset)+numel(ii2plt)):numrow_splt-1) || ...
                    (c == numrow_splt && sidx < numcol_splt) 
                subplt(c,sidx) = subplot(numrow_splt,numcol_splt,(c-1)*numcol_splt+sidx);
            end
        end
    end
    
    % plot sound waveforms
    sndPath = '\\zserver\Code\Rigging\ExpDefinitions\Tim\filmWorldwGaps\30-2019-11-22-ratio-edited\processed-merged-4-sec-w-corner-150';
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
    
    % responses to all stimuli
    % plot example cell
    mu = nanmean(dataA(mainFigP.bins<0,idxstimsubset,idxstimsubset,ii2plt(1),:),[1 2 3 5]);
    m = max(abs(mat2vec(dataA(:,idxstimsubset,idxstimsubset,ii2plt(1),:)-mu)));
    
    ma = max(mat2vec(dataM(:,idxstimsubset,idxstimsubset,ii2plt(1))));
    mi = min(mat2vec(dataM(:,idxstimsubset,idxstimsubset,ii2plt(1))));
    for vidx = 1:nstim
        v = idxstimsubset(vidx);
        if vidx == 1
            plotPSTH_fig1(squeeze(dataM(:,v,:,ii2plt(1))), ...
                idxstimsubset,mainFigP,subplt(1+vidx,1:nstim),-1,1,mainCol)
        else
            plotPSTH_fig1(squeeze(dataM(:,v,:,ii2plt(1))), ...
                idxstimsubset,mainFigP,subplt(1+vidx,1:nstim),-1,0,mainCol)
        end
        for sidx = 1:nstim
            set(gcf,'CurrentAxes', subplt(1+vidx,sidx));
            ylim([mi ma])
        end
    end
    
    % sounds for each
    clear yl
    for idx = 1:numel(ii2plt)
        plotPSTH_fig1(squeeze(nanmean(dataM(:,:,:,ii2plt(idx)),2)-nanmean(dataM(:,:,:,ii2plt(idx)),[2 3])), ...
            idxstimsubset,mainFigP,subplt(1+idx+nstim,1:nstim),wpatch,1,mainCol)
        set(gcf,'CurrentAxes', subplt(1+nstim+idx,1));
        yl{idx} = ylim; % save these limits for later
    end
    
    % videos for cell 1
    plotPSTH_fig1(squeeze(nanmean(dataM(:,:,:,ii2plt(1)),3)-nanmean(dataM(:,:,:,ii2plt(1)),[2 3])), ...
        idxstimsubset,mainFigP,subplt(1+(1:nstim),1+nstim),wpatch,1,secCol)
    ylvis = ylim;
    for s = 1:nstim % rescale
        set(gcf,'CurrentAxes', subplt(1+s,1+nstim)); hold all
        % ylim([ylvis(2)-(yl{1}(2)-yl{1}(1)), ylvis(2)])
        ylim(yl{1})
        xlim([mainFigP.bins(1),mainFigP.bins(end)])
        set(get(gca,'children'),'clipping','off')
    end
    
    % grand average
    for idx = 1:numel(ii2plt)
        plotPSTH_fig1(squeeze(nanmean(dataM(:,:,:,ii2plt(idx)),[2 3])), ...
            1,mainFigP,subplt(1+idx+nstim,nstim+1),0,1,[.0 .0 .0])
        % rescale
        ylim([min(nanmean(dataM(:,:,:,ii2plt(idx)),[2 3])),min(nanmean(dataM(:,:,:,ii2plt(idx)),[2 3]))+yl{idx}(2)-yl{idx}(1)])
        xlim([mainFigP.bins(1),mainFigP.bins(end)])
        set(get(gca,'children'),'clipping','off')
    end
    
    %% Plot raster for all cells

    % PSTH smoothing filter
    gw = gausswin(15,3);
    smWin = gw./sum(gw);
    dataM = (dataM - nanmean(dataM(mainFigP.bins<0,:,:,:),1))./nanstd(dataA(mainFigP.bins<0,:,:,:,:),[],[1 2 3 5]);
    dataM = squeeze(nanmean(dataM-nanmean(dataM,[2 3]),2));
    
    % get sorting indices
    ops.iPC = 1:min(20,size(dataM,3));
    [isort1, isort2, Sm] = mapTmap(reshape(dataM,[size(dataM,1)*size(dataM,2),size(dataM,3)])', ops);
    sortidx = isort1;
    
    for s = 1:nstim
        set(gcf,'CurrentAxes', subplt(1+idx+nstim+1,s)); hold all
        % get sound response
        m2plt = squeeze(dataM(:,idxstimsubset(s),sortidx))';
        % smooth it slightly for visu purpose
        paddedm2plt = cat(2,cat(2,repmat(m2plt(:,1),[1 numel(smWin)]),m2plt), repmat(m2plt(:,1),[1 numel(smWin)]));
        paddedm2pltSm = conv2(smWin,1,paddedm2plt', 'same')'./size(smWin,1);
        m2plt = paddedm2pltSm(:,numel(smWin)+1:end-numel(smWin));
        imagesc(mainFigP.bins,1:size(dataM,4),m2plt)
        colormap('RedBlue')
        %     caxis([-0.3,0.3])
        caxis([-prctile(mat2vec(dataM(:,idxstimsubset,:)),70) prctile(mat2vec(dataM(:,idxstimsubset,:)),70)])
        ylabel('t')
        set(gca,'visible','off')
    end
    
    %% Response for each example mice

    for kidx = 1:numel(k2plt)
        k = k2plt(kidx);
        plotPSTH_fig1(dataMPC(:,:,1,k),idxstimsubset,mainFigP,subplt(1+kidx+numel(ii2plt)+nstim+1,1:nstim),wpatch,0,mainCol)
    end
    
    %% Average over mice

    plotPSTH_fig1(squeeze(dataMPC(:,:,1,:)),idxstimsubset,mainFigP,subplt(1+1+numel(k2plt)+numel(ii2plt)+nstim+1,1:nstim),wpatch,0,mainCol)
    
    %% Sound test-retest covariance plots

    % for each example mouse
    fname = 'Cov_comp';
    focusOn = mainFocusOn;
    
    snd_cov = squeeze(nanmean(r.Cov_comp(:,:,:,strcmp(compDimP.focusOnList,focusOn)),2)./nansum(nanmean(r.Cov_comp,2), [1 4] ))*100;
    for kidx = 1:numel(k2plt)
        k = k2plt(kidx);
        set(gcf,'CurrentAxes', subplt(1+kidx+numel(ii2plt)+nstim+1,1+nstim)); hold all
        plot(compDimP.pcL(1:10),snd_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
        scatter(compDimP.pcL(1:10),snd_cov(1:10,k),30,'k','filled')
        
        hline(0,'k--')
        xlim([0 10])
        offsetAxes
        set(gca,'xticklabel',[])
    end
    
    % eigen spectrum for all mice
    set(gcf,'CurrentAxes', subplt(1+numel(k2plt)+1+numel(ii2plt)+nstim+1,1+nstim)); hold all
    for k = 1:mainFigP.nMice
        plot(compDimP.pcL(1:10),snd_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
    end
    scatter(compDimP.pcL(1:10),nanmean(snd_cov(1:10,:),2),30,'k','filled')
    hline(0,'k--')
    switch mainFocusOn
        case 'sound'
            xlabel('Auditory PCs')
        case 'video'
            xlabel('Visual PCs')
    end
    ylabel('Variance explained')
    xlim([0 10])
    xticks([1,5,10])
    offsetAxes

    % plot weights histograms
    wAll = [];
    clear f xi
    for k = 1:mainFigP.nMice
        wAll = cat(1,wAll,M{k}(:,1));
        [f{k},xi{k}] = ksdensity(M{k}(:,1));
    end
    
    % for each mouse
    for kidx = 1:numel(k2plt)
        k = k2plt(kidx);
        set(gcf,'CurrentAxes', subplt(1+kidx+numel(ii2plt)+nstim+1,1+nstim)); 
        pos = get(gca,'Position');
        axes('Parent', gcf, ...
            'Position',[pos(1) + .08, pos(2) + .05, .04, .02])
        hold all
        plot(xi{k},f{k},'color','k','Linewidth',1)
        vline(0,'k-')
        xlim([-0.2,0.4])
        set(gca,'xtick',0)
        set(gca,'yticklabel',[])
        set(gca,'ytick',[])
    end
    
    %% For all mice
    set(gcf,'CurrentAxes', subplt(1+numel(k2plt)+1+numel(ii2plt)+nstim+1,1+nstim));
    pos = get(gca,'Position');
    axes('Parent', gcf, ...
        'Position',[pos(1) + .08, pos(2) + .05, .04, .02])
    hold all
    [f_all,xi_all] = ksdensity(wAll);
    for k = 1:mainFigP.nMice
        plot(xi{k},f{k},'color',[0.5 0.5 0.5 .5],'Linewidth',1)
    end
    plot(xi_all,f_all,'color','k','Linewidth',2)
    vline(0,'k-')
    xlim([-0.2,0.4])
    set(gca,'xtick',0)
    set(gca,'yticklabel',[])
    set(gca,'ytick',[])

    %% Video test-retest covariance plots
    % for each example mouse
    fname = 'Cov_comp';
    focusOn = 'video';

    vid_cov = squeeze(nanmean(r.Cov_comp(:,:,:,strcmp(compDimP.focusOnList,focusOn)),2)./nansum(nanmean(r.Cov_comp,2), [1 4] ))*100;
    for kidx = 1:numel(k2plt)
        k = k2plt(kidx);
        set(gcf,'CurrentAxes', subplt(1+kidx+numel(ii2plt)+nstim+1,2+nstim)); hold all
        % plot real
        plot(compDimP.pcL(1:10),vid_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
        scatter(compDimP.pcL(1:10),vid_cov(1:10,k),30,'k','filled')

        hline(0,'k--')
        xlim([0 10])
        offsetAxes
        set(gca,'xticklabel',[])
    end

    % eigen spectrum for all mice
    set(gcf,'CurrentAxes', subplt(1+numel(k2plt)+1+numel(ii2plt)+nstim+1,2+nstim)); hold all
    for k = 1:mainFigP.nMice
        plot(compDimP.pcL(1:10),vid_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
    end
    scatter(compDimP.pcL(1:10),nanmean(vid_cov(1:10,:),2),30,'k','filled')
    hline(0,'k--')
    switch mainFocusOn
        case 'sound'
            xlabel('Auditory PCs')
        case 'video'
            xlabel('Visual PCs')
    end
    ylabel('Variance explained')
    % ylabel(regexprep(fname,'_',' '))
    xlim([0 10])
    xticks([1,5,10])
    offsetAxes

    % Replot the summary for main figure
    mainFigDat = load(fullfile('D:\ephys_results\processedData\audioVis\Mainfig\',sprintf('allMainfig_%s_ephysNormal_VIS.mat',mainFigP.focusOn)),'dataMPC');
    fig1Trace = squeeze(zscore(nanmean(mainFigDat.dataMPC(:,:,1,:)-nanmean(mainFigDat.dataMPC(mainFigP.bins<0,:,1,:)),4)));
    newTrace = squeeze(zscore(nanmean(dataMPC(:,:,1,:)-nanmean(dataMPC(mainFigP.bins<0,:,1,:)),4)));
    plotPSTH_fig1(fig1Trace,idxstimsubset,mainFigP,subplt(2+1+numel(k2plt)+numel(ii2plt)+nstim+1,1:nstim),wpatch,0,repmat([9 133 140]/256,[size(mainCol,1) 1]))
    plotPSTH_fig1(newTrace,idxstimsubset,mainFigP,subplt(2+1+numel(k2plt)+numel(ii2plt)+nstim+1,1:nstim),wpatch,0,mainCol)
    
    set(gcf,'CurrentAxes', subplt(2+1+numel(k2plt)+numel(ii2plt)+nstim+1,nstim+1)); hold all
    scatter(fig1Trace(:),newTrace(:),5,'k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    % colorbar
    % offsetAxes
    axis equal tight
    ylabel('New')
    xlabel('Fig1')
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    mdl = fitlm(fig1Trace(:),newTrace(:));
    corrval = corr(fig1Trace(:),newTrace(:));
    pval = mdl.anova.pValue(1);
    
    % to display the stuff that's around
    set(get(gcf,'children'),'clipping','off')
    
    %% To save as a proper svg
    
    set(gcf,'Renderer','painters');
    
    %% Plot second figure with population psths because of weird quality otherwise
    
    figure('Position',[526   313   753   379])
    % all subplots
    incr = 0;
    for c = 1
        for sidx = 1:numcol_splt
            subplt(c,sidx) = subplot(1,numcol_splt,(c-1+incr)*numcol_splt+sidx);
        end
    end
    
    prc = 65;
    for s = 1:nstim
        set(gcf,'CurrentAxes', subplt(1,s)); hold all
        % get sound response
        m2plt = squeeze(dataM(:,idxstimsubset(s),sortidx))';
        % smooth it slightly for visu purpose
        paddedm2plt = cat(2,cat(2,repmat(m2plt(:,1),[1 numel(smWin)]),m2plt), repmat(m2plt(:,1),[1 numel(smWin)]));
        paddedm2pltSm = conv2(smWin,1,paddedm2plt', 'same')'./size(smWin,1);
        m2plt = paddedm2pltSm(:,numel(smWin)+1:end-numel(smWin));
        % plot
        imagesc(mainFigP.bins,1:size(dataM,4),m2plt)
        colormap('RedBlue')
        %     caxis([-0.3,0.3])
        caxis([-prctile(mat2vec(dataM(:,idxstimsubset,:)),prc) prctile(mat2vec(dataM(:,idxstimsubset,:)),prc)])
        ylabel('t')
        set(gca,'visible','off')
        colorbar
    end