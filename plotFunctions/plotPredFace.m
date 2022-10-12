function plotPredFace(predDat,weightsPixels,rootDir)

    if ~exist('rootDir','var')
        rootDir = [];
    end

    nMice = numel(predDat.ypred_test);
    W_all = predDat.W_all;

    % Weights for best predictions
    pc2plt = 1;
    figure;
    for k = 1:nMice
        subplot(ceil(sqrt(nMice)),ceil(sqrt(nMice)),k);
        plotFace(W_all{k}(predDat.predP.lagbins == 0,:,pc2plt), weightsPixels{k});
        title(['mouse ' num2str(k)])
    end

    % PC1
    figure;
    for k = 1:nMice
        subplot(ceil(sqrt(nMice)),ceil(sqrt(nMice)),k);
        w = zeros(1,size(weightsPixels{k},3));
        w(1) = 1;
        plotFace(w, weightsPixels{k});
        title(['mouse ' num2str(k)])
    end

    % Face
    if ~isempty(rootDir)
        figure;
        for k = 1:nMice
            subplot(ceil(sqrt(nMice)),ceil(sqrt(nMice)),k);
            vface = VideoReader(fullfile(rootDir{k},'face.mj2'));
            im = read(vface,1);
            imagesc(im)
            colormap(gray)
            colorbar
            axis equal tight
            set(gca,'XTick',[], 'YTick', [])
            title(['mouse ' num2str(k)])
        end
    end
end

function wmask = plotFace(w, mask)
    s = size(mask);
    mask = reshape(mask(:,:,1:numel(w)),[s(1)*s(2),numel(w)]);
    % mask should be already normalized
    wmask = mask*w';
    wmask = reshape(wmask,[s(1),s(2)]);
    imagesc(wmask)
    lmax = max(abs(wmask(:)));
    lmin = -lmax;
    caxis([lmin,lmax]);
    colormap('redblue')
    colorbar
    axis equal tight
    set(gca,'XTick',[], 'YTick', [])
end