function delays = plotLags(w,lags)
    
    nMice = numel(w);
    
    % get the average across mice--interpolate in case not all the same time
    % bins?
    % take the first one as reference
    delays = nan(1,nMice);
    lagsRef = lags{1};
    wMean = nan(numel(lagsRef),nMice);
    for k = 1:nMice
        wMean(:,k) = interp1(lags{k},w{k}(:,1,1),lagsRef);

        % get delay
        [~,idx] = max(w{k}(:,1,1));
        delays(k) = lags{k}(idx);
    end
    wMean = nanmean((wMean-nanmean(wMean,1))./nanstd(wMean),2);

    % plot first PC weights
    figure('Position',[680   872   347   124]);
    subplot(1,3,1:2); hold all
    for k = 1:numel(w)
        fir = w{k}(:,1,1);
        fir = fir * sign(fir(find(lags{k}>0,1))-fir(1));
        plot(lags{k},(fir-mean(fir))/std(fir),'color',[.5 .5 .5 .5],'Linewidth',1)
    end
    plot(lagsRef,wMean,'color',[.0 .0 .0],'Linewidth',2)
    vline(0)
    ylabel('Weights (z-score)')
    axis tight
    xlabel('Lag (s)')
    offsetAxes

    subplot(133); hold all
    histogram(delays,20,'EdgeColor','k','FaceColor','w')
    xlim([-0.05 0.15])
    vline(0)
    xlabel('Delay to peak (s)')
    ylabel('Number of mice')
end