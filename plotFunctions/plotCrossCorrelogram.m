function [delays,valAtDelay, valAt0] = plotCrossCorrelogram(xc,lags)

    nMice = numel(xc);

    % get the average across mice--interpolate in case not all the same time
    % bins?
    % take the first one as reference
    delays = nan(1,nMice);
    valAtDelay = nan(1,nMice);
    valAt0 = nan(1,nMice);
    lagsRef = lags{1};
    xcMean = nan(numel(lagsRef),nMice);
    for k = 1:nMice
        xcMean(:,k) = interp1(lags{k},xc{k},lagsRef);

        % get delay
        [m,idx] = max(xc{k});
        delays(k) = lags{k}(idx);

        % value at delay
        valAtDelay(k) = m;

        % value at 0
        valAt0(k) = xc{k}(lags{k}==0);
    end
    xcMean = tanh(mean(atanh(xcMean),2)); % take the Z-transform as a mean

    figure('Position',[843   870   240   130]);
    subplot(1,2,1); hold all
    for k = 1:nMice
        plot(lags{k},xc{k},'color',[.5 .5 .5 .5],'Linewidth',1)
    end
    plot(lagsRef,xcMean,'color',[.0 .0 .0],'Linewidth',2)
    vline(0)
    xlabel('Lag (s)')
    ylabel('correlation')
    xlim([-1 1])
    axis tight
    offsetAxes

    subplot(122); hold all
    histogram(delays,10,'EdgeColor','k','FaceColor','w')
    xlim([-0.05 0.15])
    vline(0)
    xlabel('Delay to peak (s)')
    ylabel('Number of mice')
end