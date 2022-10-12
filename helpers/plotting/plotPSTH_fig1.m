function plotPSTH_fig1(dataM,idxstimsubset,P,subplt,wpatch,scale,col,line)

    %%% Expect a format dataM: time bins x stimulis x animals

    if ~exist('wpatch','var')
        wpatch = 1;
    end
    if ~exist('scale','var')
        scale = 0;
    end
    if ~exist('col','var')
        col = P.SndColors;
    end
    if ~exist('line','var')
        line = '-';
    end

    ma = max(mat2vec(nanmean(dataM,3)));
    mi = min(mat2vec(nanmean(dataM,3)));
    for sidx = 1:numel(idxstimsubset)
        s = idxstimsubset(sidx);
        set(gcf,'CurrentAxes', subplt(sidx));
        hold on
        switch wpatch
            case 1
                patch([0,P.bins(end),P.bins(end),0],[mi, mi, ma, ma],col(s,:),'FaceAlpha',P.plt.alph,'EdgeAlpha',.0);
                clr = [.0 .0 .0];
            case 0
                clr = col(s,:);
            case -1
                clr = [.0 .0 .0];
        end
        if size(dataM,3)>1
            for k = 1:size(dataM,3)
                plot(P.bins,dataM(:,s,k),line,'color',[clr,.1],'LineWidth',1)
            end
        end
        plot(P.bins,nanmean(dataM(:,s,:),3),line,'color',clr,'LineWidth',1)
        if ~(mi==0 && ma==0)
            ylim([mi,ma])
        end
        vline(0,'k')
        axis tight
        set(gca,'visible','off')
        if sidx==1 && scale
            if scale==1
                plot([-1.5;-1.5], [mi-0.5*abs(mi);mi-0.5*abs(mi)+20], '-k',  [-1.5;-0.5], [mi-0.5*abs(mi);mi-0.5*abs(mi)], '-k', 'LineWidth', 2)
                t = text(-1.7,mi-0.5*abs(mi), '20 sp./s', 'HorizontalAlignment','left');
            elseif scale==2
                plot([-1.5;-1.5], [mi-0.5*abs(mi);mi-0.5*abs(mi)+1], '-k',  [-1.5;-0.5], [mi-0.5*abs(mi);mi-0.5*abs(mi)], '-k', 'LineWidth', 2)
                t = text(-1.7,mi-0.5*abs(mi), '1 s.t.d', 'HorizontalAlignment','left');
            end
            set(t,'Rotation',90);
            text(-1.0,mi-0.5*abs(mi), '1 s', 'HorizontalAlignment','center')
        end
        xlim([P.bins(1),P.bins(end)])
        hold off
    end
end