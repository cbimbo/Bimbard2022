function [accVid,accSnd,p] = plotDecoding(acc,decP,sideTransec,chanceLevel,pltInset)
    
    if ~exist('pltInset','var')
        pltInset = 0;
    end
        
    idxviddim = size(acc,3); % usually the way it goes
    idxsnddim = 1;

    p(1) = signrank(acc(1,:,idxviddim),chanceLevel(1),'tail','right');
    p(2) = signrank(acc(2,:,idxsnddim),chanceLevel(2),'tail','right');
    p(3) = signrank(acc(1,:,idxviddim), acc(2,:,idxsnddim));

    figure('Position',[ 722   736   160 90]);
    hold all
    accVid = acc(strcmp(decP.focusOnL, 'video'),:,idxviddim);
    accSnd = acc(strcmp(decP.focusOnL, 'sound'),:,idxsnddim);
    scatter(accVid(sideTransec==1),accSnd(sideTransec==1),20,'k','x')
    scatter(accVid(sideTransec==-1),accSnd(sideTransec==-1),20,'k','o')
    ylabel('Sound')
    xlabel('Video')
    axis equal
    l = [0 max(ceil(max(acc(:))/5)*5,25)];
    xlim(l)
    ylim(l)
    hline(chanceLevel(2),'k--')
    vline(chanceLevel(1),'k--')
    offsetAxes
    
    if pltInset
        axes('Position',[.5 .7 .3 .3])
        hold all
        accVid = acc(strcmp(decP.focusOnL, 'video'),:,idxviddim);
        accSnd = acc(strcmp(decP.focusOnL, 'sound'),:,idxsnddim);
        scatter(accVid(sideTransec==1),accSnd(sideTransec==1),20,'k','x')
        scatter(accVid(sideTransec==-1),accSnd(sideTransec==-1),20,'k','o')
        offsetAxes
        axis equal
        xl = [floor(min(acc(1,:))/5)*5 ceil(max(acc(1,:))/5)*5];
        rangexl = diff(xl);
        yl = [floor(min(acc(2,:))/5)*5 ceil(max(acc(2,:))/5)*5];
        rangeyl = diff(xl);
        rangefinal = max(rangexl,rangeyl);
        if xl(2)>rangefinal
            xl = [xl(2)-rangefinal xl(2)];
        else
            xl = [0 rangefinal];
        end
        if yl(2)>rangefinal
            yl = [yl(2)-rangefinal yl(2)];
        else
            yl = [0 rangefinal];
        end
        xlim(xl)
        ylim(yl)
        hline(chanceLevel(2),'k--')
        vline(chanceLevel(1),'k--')
    end
    
    %%% to display the stuff that's around
    set(get(gcf,'children'),'clipping','off')
end