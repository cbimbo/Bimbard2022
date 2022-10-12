function [percExplSndMean,percExplVidMean,p] = plotSubspaceOverlap(subOvDat,sideTransec)

    Cov_compVid = subOvDat.Cov_comp_vid;
    percExplVid = subOvDat.percExpl_vid;
    Cov_compSnd = subOvDat.Cov_comp_snd;
    percExplSnd = subOvDat.percExpl_snd;

    vidColor = [0.0417    0.6250    0.0417];
    sndColor = [0.6250    0.0417    0.6250];
    behColor = [0 0 1];

    nMice = size(Cov_compVid.vid,3);
    if ~exist('sideTransec','var')
        sideTransec = -ones(1,nMice);
    end

    %%
    for k = 1:nMice
        figure('Position',[676   765   211   213]);

        %% VIDEO RELATED VARIANCE
        normFactor = sum(nanmean(Cov_compVid.vid(:,:,k),2))/100; % max video related variance
        xrand = [1:4 4:-1:1];
        yrand = [prctile(nanmean(Cov_compVid.rand(:,:,k,:),2)/normFactor,5,4)' prctile(nanmean(Cov_compVid.rand(:,:,k,:),2)/normFactor,95,4)'];
        
        % plot vid
        ax(1) = subplot(2,3,1); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compVid.vid,1),nanmean(Cov_compVid.vid(:,:,k),2)/normFactor,20,vidColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Vis. PCs')
        ylabel({'Video-related';'variance'})
        offsetAxes

        % plot snd
        ax(2) = subplot(2,3,2); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compVid.snd,1),nanmean(Cov_compVid.snd(:,:,k),2)/normFactor,20,sndColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Aud. PCs')
        offsetAxes

        % plot mov
        ax(3) = subplot(2,3,3); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compVid.mov,1),nanmean(Cov_compVid.mov(:,:,k),2)/normFactor,20,behColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Beh. PCs')
        linkaxes(ax,'y')
        offsetAxes

        %% SOUND RELATED VARIANCE
        normFactor = sum(nanmean(Cov_compSnd.snd(:,:,k),2))/100; % max sound related variance
        xrand = [1:4 4:-1:1];
        yrand = [prctile(nanmean(Cov_compSnd.rand(:,:,k,:),2)/normFactor,5,4)' prctile(nanmean(Cov_compSnd.rand(:,:,k,:),2)/normFactor,95,4)'];

        % plot vid
        ax(1) = subplot(2,3,4); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compSnd.vid,1),nanmean(Cov_compSnd.vid(:,:,k),2)/normFactor,20,vidColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Vis. PCs')
        ylabel({'Sound-related';'variance'})
        offsetAxes

        % plot snd
        ax(2) = subplot(2,3,5); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compSnd.snd,1),nanmean(Cov_compSnd.snd(:,:,k),2)/normFactor,20,sndColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Aud. PCs')
        offsetAxes

        % plot mov
        ax(3) = subplot(2,3,6); hold all
        patch(xrand,yrand,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.4,'EdgeAlpha',0.4)
        scatter(1:size(Cov_compSnd.mov,1),nanmean(Cov_compSnd.mov(:,:,k),2)/normFactor,20,behColor,'filled')
        hline(0,'k--')
        axis tight
        xlim([0 4])
        xlabel('Beh. PCs')
        linkaxes(ax,'y')
        offsetAxes
    end

    %% Individual mouse overlap
    figure('Position', [1120         814         120         164])
    subplot(121); hold all
    for k = 1:nMice
        if size(Cov_compSnd.rand,4)>0
            scatter(k+(0.3*rand(1,size(Cov_compSnd.rand,4))-0.15),percExplVid.rand(k,:),10,[0.5 0.5 0.5],'filled')
        end
        scatter(k,percExplVid.mov(k),20,'k')
    end
    ylim([0,100])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('mouse')
    ylabel('% overlap with video sub')

    subplot(122); hold all
    for k = 1:nMice
        if size(Cov_compSnd.rand,4)>0
            scatter(k+(0.3*rand(1,size(Cov_compSnd.rand,4))-0.15),percExplSnd.rand(k,:),10,[0.5 0.5 0.5],'filled')
        end
        scatter(k,percExplSnd.mov(k),20,'k')
    end
    ylim([0,100])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xlabel('mouse')
    ylabel('% overlap with sound sub')

    %% Summary

    figure('Position', [ 1120         827         120         151]); hold all
    % VIDEO
    scatter(rand(1,sum(sideTransec == -1))*0.5,percExplVid.mov(sideTransec == -1),20,vidColor)
    scatter(0.7,nanmean(percExplVid.mov(sideTransec == -1)),30,vidColor,'filled')
    scatter(rand(1,sum(sideTransec == 1))*0.5,percExplVid.mov(sideTransec == 1),20,vidColor,'x')
    scatter(0.7,nanmean(percExplVid.mov(sideTransec == 1)),30,vidColor,'x')
    if size(Cov_compVid.rand,4)>0
        for k = 1:nMice
            plot([0 0.5],ones(1,2)*prctile(percExplVid.rand(k,:),95),'color',vidColor)
        end
    end

    % SOUND
    scatter(1+rand(1,sum(sideTransec == -1))*0.5,percExplSnd.mov(sideTransec == -1),20,sndColor)
    scatter(1+0.7,nanmean(percExplSnd.mov(sideTransec == -1)),30,sndColor,'filled')
    scatter(1+rand(1,sum(sideTransec == 1))*0.5,percExplSnd.mov(sideTransec == 1),20,sndColor,'x')
    scatter(1+0.7,nanmean(percExplSnd.mov(sideTransec == 1)),30,sndColor,'x')
    if size(Cov_compVid.rand,4)>0
        for k = 1:nMice
            plot(1+[0 0.5],ones(1,2)*prctile(percExplSnd.rand(k,:),95),'color',sndColor)
        end
    end

    ylim([0,100])
    xlim([0 2.0])
    offsetAxes
    set(gca,'xtick',[0.3 1.3])
    set(gca,'xticklabel',{'vis.','aud.'})
    ylabel('% overlap')

    percExplSndMean = nanmean(percExplSnd.mov);
    percExplVidMean = nanmean(percExplVid.mov);
    p = signrank(percExplSnd.mov,percExplVid.mov);
    sigstar({[0.3 1.3]},p)