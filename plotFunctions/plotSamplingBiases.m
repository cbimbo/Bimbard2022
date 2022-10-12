%% Check single cell responses

for  k = 1:P.nMice
    P.focusOn = 'video'; % so that order is unchanged
    dataA = getOrgData(db(k).spikeData,db(k).audiovisuoCode,P);

    % will recenter the data by removing the common response
    %%% note that this is not specific to each unimportant stim
    dataA = dataA - nanmean(dataA(:,:,:,:,1:2:end), [2 3 5]);

    nClu = size(dataA,4);

    c_acrossmod{k} = nan(2,nClu);
    
    % video reliability
    x1 = reshape(nanmean(dataA(:,:,:,:,1:2:end),[3 5]), [P.nBins*P.nGroupsVid,nClu]);
    x2 = reshape(nanmean(dataA(:,:,:,:,2:2:end),[3 5]), [P.nBins*P.nGroupsVid,nClu]);
    for c = 1:nClu
        c_acrossmod{k}(1,c) = corr(x1(:,c),x2(:,c));
    end

    % audio reliability
    x1 = reshape(nanmean(dataA(:,:,:,:,1:2:end),[2 5]), [P.nBins*P.nGroupsSnd,nClu]);
    x2 = reshape(nanmean(dataA(:,:,:,:,2:2:end),[2 5]), [P.nBins*P.nGroupsSnd,nClu]);
    for c = 1:nClu
        c_acrossmod{k}(2,c) = corr(x1(:,c),x2(:,c));
    end

    figure;
    prettyScatter(c_acrossmod{k});
    xlabel('corr. visual')
    ylabel('corr. auditory')

    % Compute other single cell variables
    Depths{k} = max(db(k).C.Depth) - db(k).C.Depth;
    DepthsL4{k} = db(k).C.DepthFromTopL4;
    FRs{k} = squeeze(nanmean(db(k).spikeData(P.bins<0,:,:),[1 3]));
    ISIviol{k} = db(k).C.ISIViol;
    ISIviol{k}(ISIviol{k} == 0) = 0.00001; % to avoid issue when plotting the log
end

%%

Q2plt = FRs;
% plot as a function of ISI violations
figure('Position',[ 893   781   347   197])
hold all
clear b slop corrCoeff
for k = 1:P.nMice
    x = Q2plt{k}';
    y = c_acrossmod{k}(2,:)';
    scatter(x,y,20,P.mouseColor(k,:))
    X = [x, ones(numel(x),1)];
    b{k} = X\y;
    plot(x,X*b{k},'color',P.mouseColor(k,:),'LineWidth',2)
    slope(k) = b{k}(1);
    corrCoeff(k) = corr(Q2plt{k}',c_acrossmod{k}(2,:)','Type','Spearman');
    % set(gca, 'XScale', 'log')
end
% xlabel('Depths')
xlabel('Baseline firing rate')
% xlabel('ISI violations')
% xlim([0 1]);
ylabel({'Test-retest correlation';'of auditory response'} )

% inset with mice histograms
pos = get(gca,'Position');
axes('Parent', gcf, ...
    'Position',[pos(1) + .5, pos(2) + .5, .2, .2])
hold all
histogram(slope,'BinEdges',linspace(min([slope 0]),max([slope 0]),10),'DisplayStyle','stairs','EdgeColor','k')
vline(0,'k--')
xlabel('Slope')
% scatter(cell2mat(cellfun(@(x) nanmean(x(2,:)), c_acrossmod, 'uni', 0)),slope,20,P.mouseColor)
% ylabel('Slope')
% xlabel('Average test-retest correlation')