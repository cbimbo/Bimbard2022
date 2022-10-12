function cov = plotEigenspectra(r,compDimP,sideTransec)
    %%% This function will plot the eigenspectra 
   
    cov = squeeze(nanmean(r.Cov_comp ./ nansum(r.Cov_comp,[1 4]),2))*100;
    
    figure('Position',[680   380   500   330]);
    for f = 1:numel(compDimP.focusOnList)
        focusOn = compDimP.focusOnList{f};
        
        
        % Normal scale
        subplot(2,numel(compDimP.focusOnList),f)
        plotEigenspectrum(cov(:,:,strcmp(compDimP.focusOnList,focusOn)),compDimP,sideTransec)
        xlim([0 10])
        xticks([1,5,10,20])
        ylim([-2,ceil(max(cov(:))/5)*5])
        title(compDimP.focusOnList{f})
        ylabel('Test-retest covariance (%)')
        
        % Log scale
        subplot(2,numel(compDimP.focusOnList),numel(compDimP.focusOnList)+f)
        plotEigenspectrum(cov(:,:,strcmp(compDimP.focusOnList,focusOn)),compDimP,sideTransec)
        xlabel('PCs')
        set(gca,'yscale','log')
        xlim([0 10])
        xticks([1,5,10,20])
        ylim([0.01,100])
        yticks([0.1 1 10 50])
    end
end

function plotEigenspectrum(cov,compDimP,sideTrans)
    % plot real
    hold on
    for k = 1:size(cov,2)
        % plot each
        if sideTrans(k) == -1
            plot(compDimP.pcL,cov(:,k),'LineWidth',0.5, 'color', [.5 .5 .5 .5])
        else
            plot(compDimP.pcL,cov(:,k),'--','LineWidth',0.5, 'color', [.5 .5 .5 .5])
        end
    end
    scatter(compDimP.pcL,nanmean(cov(:,sideTrans == -1),2),30,'k','filled')
    scatter(compDimP.pcL,nanmean(cov(:,sideTrans == 1),2),30,'k','x')
    hline(0,'k--')
    hold off
end
    
    