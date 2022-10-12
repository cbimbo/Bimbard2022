function [c,p] = plotVersus(acc1,acc2,sideTrans,chancelevel,pltMdl,l,inSbplt)
    
    if ~exist('pltMdl','var')
        pltMdl = 0;
    end

    if ~exist('l','var')
        l = [5 max(ceil(max([acc1(:); acc2(:)])/5)*5,30)];
    end

    if ~exist('inSbplt','var')
        inSbplt = 0;
    end

    if ~inSbplt
        figure('Position',[722   750   120   100]);
    end
    hold all
    scatter(acc1(sideTrans == -1),acc2(sideTrans == -1),20,'k')
    scatter(nanmean(acc1(sideTrans == -1)),nanmean(acc2(sideTrans == -1)),30,'k','filled')
    if sum(sideTrans == 1)>0
        scatter(acc1(sideTrans == 1),acc2(sideTrans == 1),20,'k','x')
        scatter(nanmean(acc1(sideTrans == 1)),nanmean(acc2(sideTrans == 1)),30,'k','x')
    end
    
    xlabel('accuracy 1')
    ylabel('accuracy 2')
    axis equal tight
    xlim(l)
    ylim(l)
    if exist('chancelevel','var') & ~isempty(chancelevel)
        hline(chancelevel(2),'k--')
        vline(chancelevel(1),'k--')
    end
    plot(l,l,'k--')
    if pltMdl
        mdl = fitlm(acc1',acc2');
        p = mdl.anova.pValue(1);
        c = corr(acc1',acc2');
        title(sprintf('p-value: %d',p))
        plot(acc1,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*acc1,'k')
    else 
        p = [];
    end
    offsetAxes
    
end