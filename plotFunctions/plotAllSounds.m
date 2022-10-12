function plotAllSounds(mainFigDat,behFigDat)
    neurTrace = squeeze(zscore(nanmean(mainFigDat.dataMPC(:,:,1,:)-nanmean(mainFigDat.dataMPC(mainFigDat.mainFigP.bins<0,:,1,:)),4)));
    behTrace = squeeze(zscore(nanmean(behFigDat.dat(find(strcmp({behFigDat.dat.name},'Body motion'))).dataMPC(:,:,1,:)- ...
        nanmean(behFigDat.dat(find(strcmp({behFigDat.dat.name},'Body motion'))).dataMPC(behFigDat.behFigP.bins<0,:,1,:)),4)));

    %% plot

    sndPath = '\\zserver\Code\Rigging\ExpDefinitions\Tim\filmWorldwGaps\30-2019-11-22-ratio-edited\processed-merged-4-sec-w-corner-150';
    sndOrder = [27 9 18 25 26 17 15 11 29 22 7 0];

    figure('Position',[680   226   560   752]);
    for s = 1:mainFigDat.mainFigP.nGroupsSnd
        subplot(floor(sqrt(mainFigDat.mainFigP.nGroupsSnd)),ceil(sqrt(mainFigDat.mainFigP.nGroupsSnd)),s); hold all

        sidx = find(mainFigDat.mainFigP.labelsGroupsSnd == sndOrder(s));
        [audio, Fs] = audioread(fullfile(sndPath, ['video-0-audio-' num2str(mainFigDat.mainFigP.labelsGroupsSnd(sidx)) '.mp4']));
        snd = sum(audio,2);
        [eupper,elower] = envelope(10*cat(1,zeros(0.135*Fs,1),snd),100,'rms');
        x = 1/Fs:1/Fs:(size(audio,1)/Fs+0.135);
        x = x(1:10:end); elower = elower(1:10:end); eupper = eupper(1:10:end);
        x2 = [x, fliplr(x)];
        inBetween = [elower', fliplr(eupper')]*0.5;
        h = fill(x2, inBetween,[.0 .0 .0]);
        set(h,'facealpha',.5,'edgealpha',.0)
        plot(mainFigDat.mainFigP.bins,zeros(1,numel(mainFigDat.mainFigP.bins)),'color',[0.5 0.5 0.5],'LineWidth',1)

        plot(mainFigDat.mainFigP.bins,neurTrace(:,sidx)-5,'color',nanmean(mainFigDat.mainFigP.SndColors),'LineWidth',1)
        plot(mainFigDat.mainFigP.bins,behTrace(:,sidx)-10,'color',[.0 .0 1.0],'LineWidth',1)
        ylim([-15,5])
        xlim([-1.0,3.8])
        set(gca,'Visible','off')

        if s == 1
            plot([-1 0], [-12 -12], 'k','LineWidth',2)
            text(-0.5,-13,'1s','horizontalalignment', 'center')
        end
    end

    set(get(gcf,'children'),'clipping','off')
    set(gcf,'Renderer','painters');
end