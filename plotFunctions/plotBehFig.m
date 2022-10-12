function [c,p] = plotBehFig(dat,behFigP)

    k2plt = behFigP.k2plt;
    
    %%
    
    behFigP.plt.alph = .4;
    wpatch = 0;
    mainFocusOn = 'sound';
    switch dat.name
        case 'Body motion'
            mainCol = ones(size(behFigP.SndColors)).*[0 0 1];
        otherwise
            mainCol = ones(size(behFigP.SndColors)).*[0.15  0.5 1.0];
    end
    % what stimuli
    idxstimsubset = [];
    for idx = 1:numel(behFigP.stimsubset)
        idxstimsubset = [idxstimsubset, find(behFigP.labelsGroupsSnd == behFigP.stimsubset(idx))];
    end
    nstim = numel(idxstimsubset);
    
    %% plot something like ~fig1 for blink, pupil area, saccades or motion
    
    numrow_splt = 4;
    numcol_splt = nstim+1;
    figure('Position',[549    32   630   350]);
    % all subplots
    for c = 1:numrow_splt
        for sidx = 1:numcol_splt
            if ~(sidx == numcol_splt) || c == numrow_splt
                subplt(c,sidx) = subplot(numrow_splt,numcol_splt,(c-1)*numcol_splt+sidx);
            end
        end
    end
    
    % plot sound waveforms
    sndPath = '\\zserver\Code\Rigging\ExpDefinitions\Tim\filmWorldwGaps\30-2019-11-22-ratio-edited\processed-merged-4-sec-w-corner-150';
    for sidx = 1:nstim
        set(gcf,'CurrentAxes', subplt(1,sidx));
        [audio, Fs] = audioread(fullfile(sndPath, ['video-0-audio-' num2str(behFigP.labelsGroupsSnd(idxstimsubset(sidx))) '.mp4']));
        snd = sum(audio,2);
        plot(1/Fs:1/Fs:(size(audio,1)/Fs+0.135),cat(1,zeros(0.135*Fs,1),snd),'color',[.5 .5 .5],'LineWidth',1)
        xlim([behFigP.bins(1),behFigP.bins(end)])
        ylim([-2,2])
        ylabel('t') % to make sure all has the same size
        set(gca,'visible','off')
    end
    
    % PLOT FACE MOTION
    dataMPC = dat.dataMPC;
    plotPSTH_fig1(squeeze(dataMPC(:,:,1,k2plt)),idxstimsubset,behFigP,subplt(2,1:nstim),wpatch,2,mainCol) % for example mouse
    plotPSTH_fig1(squeeze(dataMPC(:,:,1,:)),idxstimsubset,behFigP,subplt(3,1:nstim),wpatch,2,mainCol) % average over mice
    
    % PLOT FACE MOTION AND NEURAL ACTIVITY
    % replot the summary for main figure
    mainFigDat = load(fullfile('D:\ephys_results\processedData\audioVis\Mainfig\',sprintf('allMainfig_%s_ephysNormal_VIS.mat',behFigP.focusOn)),'dataMPC');
    fig1Trace = squeeze(zscore(nanmean(mainFigDat.dataMPC(:,:,1,:)-nanmean(mainFigDat.dataMPC(behFigP.bins<0,:,1,:)),4)));
    newTrace = squeeze(zscore(nanmean(dataMPC(:,:,1,:)-nanmean(dataMPC(behFigP.bins<0,:,1,:)),4)));
    plotPSTH_fig1(fig1Trace,idxstimsubset,behFigP,subplt(4,1:nstim),wpatch,0,repmat([9 133 140]/256,[size(mainCol,1) 1]))
    plotPSTH_fig1(newTrace,idxstimsubset,behFigP,subplt(4,1:nstim),wpatch,0,mainCol)
    
    % plot correlation
    set(gcf,'CurrentAxes', subplt(4,nstim+1)); hold all
    scatter(fig1Trace(:),newTrace(:),5,'k','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    % offsetAxes
    axis equal tight
    ylabel('Motion')
    xlabel('VIS')
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    mdl = fitlm(fig1Trace(:),newTrace(:));
    c = corr(fig1Trace(:),newTrace(:));
    p = mdl.anova.pValue(1);
    title(sprintf('Corr%d / p-value: %d',c,p))
    
    %%% to display the stuff that's around
    set(get(gcf,'children'),'clipping','off')