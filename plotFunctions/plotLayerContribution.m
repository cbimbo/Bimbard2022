
load(fullfile(preprocFolder,sprintf('granularLayerBounds_ephys%s',what2Load)));
granBound = granBoundManual;

%% Save the cell's location

for k = 1:P.nMice
    shankPos = ceil(db(k).C.XPos/200);
    shanks = unique(shankPos);
    db(k).C.Layer = ones(1,numel(db(k).C.Depth))*2;
    for s = 1:numel(shanks)
        if ~any(isnan(granBound{k}{s}))
            db(k).C.Layer(shankPos == shanks(s) & db(k).C.Depth > granBound{k}{s}(1)) = 1;
            db(k).C.Layer(shankPos == shanks(s) & db(k).C.Depth < granBound{k}{s}(2)) = 3;
            db(k).C.DepthFromTopL4(shankPos == shanks(s)) = granBound{k}{s}(1)-db(k).C.Depth(shankPos == shanks(s));
        else
            db(k).C.Layer(shankPos == shanks(s)) = nan;
            db(k).C.DepthFromTopL4(shankPos == shanks(s)) = nan;
        end
    end
end

%% Plot example animals
% Example animal 1
k = 1;
plotLFPandCSD(LFP{k},CSD{k},depthChan{k},timebinsAll{k}, ...
    db_bu(k).C.XPos(db_bu(k).C.CluLab == 2),db_bu(k).C.Depth(db_bu(k).C.CluLab == 2), ...
    P.fileref.anatDir{k}, granBound{k});

% Example animal 2
k = 4;
plotLFPandCSD(LFP{k},CSD{k},depthChan{k},timebinsAll{k}, ...
    db_bu(k).C.XPos(db_bu(k).C.CluLab == 2),db_bu(k).C.Depth(db_bu(k).C.CluLab == 2), ...
    P.fileref.anatDir{k}, granBound{k});

%% Plot LFP, CSD and neuron count as a function of norm. depth

binLim = [-300 700];

% Get mean LFP and CSD
binGapLFP = 25;
binDepthsLFP = binLim(1):binGapLFP:binLim(2);
binPosLFP = binDepthsLFP(1:end-1)+binGapLFP/2;
timebins = timebinsAll{1}; % should all be the same?
LFPMean = [];
CSDMean = [];
d2Infra = [];
for k = 1:P.nMice
    for s = 1:numel(LFP{k})
        if ~any(isnan(granBound{k}{s}(1)))
            LFPMean = cat(3,LFPMean,interp1(granBound{k}{s}(1)-depthChan{k}{s},LFP{k}{s},binPosLFP));
            CSDMean = cat(3,CSDMean,interp1(granBound{k}{s}(1)-depthChan{k}{s}(3:end-2),CSD{k}{s},binPosLFP));
            d2Infra = cat(2,d2Infra,-diff(granBound{k}{s}));
        end
    end
end
LFPMean = nanmean(LFPMean,3);
CSDMean = nanmean(CSDMean,3);

figure('Position',[579   583   561   295]);
% Plot LFP
ax(1) = subplot(141);
imagesc(timebins,binLim,interp2(LFPMean,1))
set(gca,'YDir','normal')
colormap(ax(1),'RedBlue'); caxis([-max(abs(LFPMean(:))) max(abs(LFPMean(:)))])
hline(0,'k--')
hline(nanmean(d2Infra),'k--')
xlabel('Time (s)')
ylabel('Depths from top of L4 (μm)')
set(gca,'YDir','reverse')
set(gca,'box','off')
xlim([-0.01 0.08])

% Plot CSD
ax(2) = subplot(142);
imagesc(timebins,binLim,interp2(CSDMean,1))
set(gca,'YDir','normal')
colormap(ax(2),'RedBlue'); caxis([-max(abs(CSDMean(:))) max(abs(CSDMean(:)))])
hline(0,'k--')
hline(nanmean(d2Infra),'k--')
xlabel('Time (s)')
set(gca,'YDir','reverse')
set(gca,'box','off')
set(gca,'YTick',[])
xlim([-0.01 0.08])

% Plot neuron Count
ax(3) = subplot(143);
binGap = 50;
binDepths = binLim(1):binGap:binLim(2);
binPos = binDepths(1:end-1)+binGap/2;
neuronCount = zeros(1,numel(binPos));
hold all
for k = 1:P.nMice
    h = histcounts(db(k).C.DepthFromTopL4,binDepths);
    newNeuronCount = neuronCount + h;
    x = [binPos, fliplr(binPos)];
    inBetween = [neuronCount, fliplr(newNeuronCount)];
    fill(inBetween,x,P.mouseColor(k,:),'FaceAlpha',0.5)
    neuronCount = newNeuronCount;
end
set(gca,'YDir','reverse')
xlabel('Neuron count')
set(gca,'box','off')
set(gca,'YTick',[])

ax(4) = subplot(144);
% Plot sampling bias as a function of depth
hold all
clear b slop corrCoeff
for k = 1:P.nMice
    x = db(k).C.DepthFromTopL4';
    y = c_acrossmod{k}(2,:)';
    scatter(y,x,20,P.mouseColor(k,:),'MarkerEdgeAlpha',.2)
    X = [x, ones(numel(x),1)];
    b{k} = X\y;
    plot(X*b{k},x,'color',P.mouseColor(k,:),'LineWidth',2)
    slope(k) = b{k}(1);
    corrCoeff(k) = corr(x,y,'Type','Spearman');
end
vline(0,'k--')
set(gca, 'YDir', 'reverse')
xlabel({'Test-retest correlation';'of auditory response'} )
ylim(binLim)
set(gca,'box','off')
set(gca,'YTick',[])

% inset with mice histograms
pos = get(gca,'Position');
axes('Parent', gcf, ...
    'Position',[pos(1) + .1, pos(2) + .7, .1, .1])
hold all
histogram(slope,'BinEdges',linspace(min([slope 0]),max([slope 0]),10),'DisplayStyle','stairs','EdgeColor','k')
vline(0,'k--')
xlabel('Slope')

linkaxes(ax,'y')