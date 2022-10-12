function [c_ContravsIpsi,p_ContravsIpsi,c_Contra,p_Contra,c_Ipsi,p_Ipsi] = plotTransFig(transFigDat,eigDat)
    
    %% Rename everything

    mainFigP = transFigDat.mainFigP;
    dataMPC = transFigDat.dataMPC;
    M = transFigDat.M; 
    dataA = transFigDat.dataA;
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
        case 'video'
            mainCol = ones(size(mainFigP.VidColors)).*nanmean(mainFigP.VidColors);
    end
    
    % what stimuli
    idxstimsubset = [];
    for idx = 1:numel(mainFigP.stimsubset)
        idxstimsubset = [idxstimsubset, find(mainFigP.labelsGroupsSnd == mainFigP.stimsubset(idx))];
    end
    nstim = numel(idxstimsubset);
    
    %%
    
    figure('Position',[549    32   800   300])
    % figure('Position',[549    32   630   297]);
    numcol_splt = nstim + 2;
    % all subplots
    for c = 1:4
        for sidx = 1:numcol_splt
            if (c == 1 && sidx < numcol_splt-1) || (ismember(c,[2 3]) && sidx < numcol_splt) || c == 4
                subplt(c,sidx) = subplot(4,numcol_splt,(c-1)*numcol_splt+sidx);
            end
        end
    end
    
    % plot stims
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

    % plot responses
    plotPSTH_fig1(squeeze(dataMPC(:,:,1,transFigDat.sideTransec==-1)),idxstimsubset,mainFigP,subplt(2,1:nstim),wpatch,0,mainCol)
    plotPSTH_fig1(squeeze(dataMPC(:,:,1,transFigDat.sideTransec==1)),idxstimsubset,mainFigP,subplt(3,1:nstim),wpatch,0,mainCol)

    %% Eigen spectrum
    fname = 'Cov_comp';
    focusOn = mainFocusOn;
    % Get weights
    clear f xi
    for k = 1:mainFigP.nMice
        [f{k},xi{k}] = ksdensity(M{k}(:,1));
    end
    
    snd_cov = squeeze(nanmean(r.Cov_comp(:,:,:,strcmp(compDimP.focusOnList,focusOn)),2)./nansum(nanmean(r.Cov_comp,2), [1 4] ))*100;

    % eigen spectrum ipsi
    set(gcf,'CurrentAxes', subplt(2,1+nstim)); hold all
    % % plot real
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == 1
            plot(compDimP.pcL(1:10),snd_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
        end
    end
    scatter(compDimP.pcL(1:10),nanmean(snd_cov(1:10,:),2),30,'k','filled')
    hline(0,'k--')
    xlim([0 10])
    xticks([1,5,10])
    offsetAxes
    
    % plot weights histograms
    wAll = [];
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == 1
            wAll = cat(1,wAll,M{k}(:,1));
        end
    end
    set(gcf,'CurrentAxes', subplt(2,1+nstim));
    pos = get(gca,'Position');
    axes('Parent', gcf, ...
        'Position',[pos(1) + .06, pos(2) + .1, pos(3)/2, pos(4)/2])
    hold all
    [f_all,xi_all] = ksdensity(wAll);
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == 1
            plot(xi{k},f{k},'color',[0.5 0.5 0.5 .5],'Linewidth',1)
        end
    end
    plot(xi_all,f_all,'color','k','Linewidth',2)
    vline(0,'k-')
    xlim([-0.2,0.4])
    set(gca,'xtick',0)
    set(gca,'yticklabel',[])
    set(gca,'ytick',[])
    
    % eigen spectrum contra
    set(gcf,'CurrentAxes', subplt(3,1+nstim)); hold all
    % % plot real
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == -1
            plot(compDimP.pcL(1:10),snd_cov(1:10,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
        end
    end
    scatter(compDimP.pcL(1:10),nanmean(snd_cov(1:10,:),2),30,'k','filled')
    hline(0,'k--')
    xlim([0 10])
    xticks([1,5,10])
    offsetAxes
    
    % plot weights histograms
    wAll = [];
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == -1
            wAll = cat(1,wAll,M{k}(:,1));
        end
    end
    set(gcf,'CurrentAxes', subplt(3,1+nstim));
    pos = get(gca,'Position');
    axes('Parent', gcf, ...
        'Position',[pos(1) + .06, pos(2) + .1, pos(3)/2, pos(4)/2])
    hold all
    [f_all,xi_all] = ksdensity(wAll);
    for k = 1:mainFigP.nMice
        if transFigDat.sideTransec(k) == 1
            plot(xi{k},f{k},'color',[0.5 0.5 0.5 .5],'Linewidth',1)
        end
    end
    plot(xi_all,f_all,'color','k','Linewidth',2)
    vline(0,'k-')
    xlim([-0.2,0.4])
    set(gca,'xtick',0)
    set(gca,'yticklabel',[])
    set(gca,'ytick',[])
    
    %%% replot the summary for main figure
    mainFigDat = load(fullfile('D:\ephys_results\processedData\audioVis\Mainfig\',sprintf('allMainfig_%s_ephysNormal_VIS.mat',mainFigP.focusOn)),'dataMPC');
    fig1Trace = squeeze(zscore(nanmean(mainFigDat.dataMPC(:,:,1,:)-nanmean(mainFigDat.dataMPC(mainFigP.bins<0,:,1,:)),4)));
    newTraceIpsi = squeeze(zscore(nanmean(dataMPC(:,:,1,transFigDat.sideTransec==1)-nanmean(dataMPC(mainFigP.bins<0,:,1,transFigDat.sideTransec==1)),4)));
    newTraceContra = squeeze(zscore(nanmean(dataMPC(:,:,1,transFigDat.sideTransec==-1)-nanmean(dataMPC(mainFigP.bins<0,:,1,transFigDat.sideTransec==-1)),4)));
    plotPSTH_fig1(fig1Trace,idxstimsubset,mainFigP,subplt(4,1:nstim),wpatch,0,repmat([9 133 140]/256,[size(mainCol,1) 1]))
    plotPSTH_fig1(newTraceIpsi,idxstimsubset,mainFigP,subplt(4,1:nstim),wpatch,0,mainCol,'--')
    plotPSTH_fig1(newTraceContra,idxstimsubset,mainFigP,subplt(4,1:nstim),wpatch,0,mainCol)
    
    % Plot correlation
    set(gcf,'CurrentAxes', subplt(4,1+nstim)); hold all
    scatter(newTraceContra(:),newTraceIpsi(:),5,'k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    axis equal tight
    ylabel('Ipsi')
    xlabel('Contra')
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    mdl_ContravsIpsi = fitlm(newTraceContra(:),newTraceIpsi(:));
    p_ContravsIpsi = mdl_ContravsIpsi.anova.pValue(1);
    c_ContravsIpsi = corr(newTraceContra(:),newTraceIpsi(:));

    % Plot correlation
    set(gcf,'CurrentAxes', subplt(4,2+nstim)); hold all
    scatter(fig1Trace(:),newTraceIpsi(:),5,'k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    axis equal tight
    ylabel('Ipsi')
    xlabel('Fig1')
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    mdl_Contra = fitlm(fig1Trace(:),newTraceContra(:));
    p_Contra = mdl_Contra.anova.pValue(1);
    c_Contra = corr(fig1Trace(:),newTraceContra(:));
    mdl_Ipsi = fitlm(fig1Trace(:),newTraceIpsi(:));
    p_Ipsi = mdl_Ipsi.anova.pValue(1);
    c_Ipsi = corr(fig1Trace(:),newTraceIpsi(:));
end