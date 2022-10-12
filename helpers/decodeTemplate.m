function [accuracy, distanceMean] = decodeTemplate(UTrain_r, UTest_r, plt)

% will decode based on the template obtained on a train set
% NB: takes into account the strength of the projection on each pc (ie,
% the eigenvalues)
% UTrain_r and UTest_r should be nBins * nGroups * nPCs * nRep

if nargin < 3
    plt = 0;
end

nBins = size(UTrain_r,1);
nGroups = size(UTrain_r,2);
% nRepTrain = size(UTrain_r,4);
nRepTest = size(UTest_r,4);
nPCs = size(UTrain_r,3);

% compute template matching for each left out
Templates = reshape(permute(nanmean(UTrain_r,4), [1 3 2]), [nBins*nPCs, nGroups]);
what2Match = reshape(permute(UTest_r, [1 3 2 4]), [nBins*nPCs, nGroups, nRepTest]);

idxmin = nan(nGroups, nRepTest);
distances = nan(nGroups, nGroups, nRepTest);
cc = nan(nGroups, nGroups, nRepTest);
for r = 1:nRepTest % loop over nonimportant stimuli (=repeats)
    for rr = 1:nGroups % loop over important stimuli
        distances(:,rr,r) = sum((what2Match(:,rr,r) - Templates).^2 ,1);
%         cc(:,rr,r) = corr(what2Match(:,rr,r), Templates);
    end
    
    [~,idxmin(:,r)] = min(distances(:,:,r));
%     [~,idxmin(:,r)] = max(cc(:,:,r));
end
accuracy = sum(sum(idxmin' == 1:nGroups))/numel(idxmin);

distanceMean = nan(nGroups, nGroups);
for rr = 1:nGroups
    distanceMean(:,rr) = sum((mean(what2Match(:,rr,:),3) - Templates).^2);
end

if plt
    % plot for mean over all unimportant stim in last repeat
    figure('Position', [680   616   929   362]);
    colormap(magma)
    
    cmin = min([Templates(:); mat2vec(mean(what2Match,3))]);
    cmax = max([Templates(:); mat2vec(mean(what2Match,3))]);
    subplot(2,3,[1 4]); 
    imagesc(Templates); 
    caxis([cmin cmax])
    title('Template')
    
    subplot(2,3,[2 5]); 
    imagesc(mean(what2Match,3)); 
    caxis([cmin cmax])
    title('Test')
    
    subplot(2,3,3);
    [~,idxminMean] = sort(distanceMean,1);
    imagesc(distanceMean); hold on
    colormap(magma)
    c = colorbar;
    c.Label.String = 'Error';
    scatter(1:nGroups, idxminMean(1,:), 80, 'o', 'w')
    scatter(1:nGroups, idxminMean(2,:), 40, 'o', 'w')
    axis equal tight
    
    subplot(2,3,6);
    imagesc(idxmin)
    axis equal tight
    ylabel('stim')
    xlabel('presentation #')
    c = colorbar;
    c.Label.String = 'Decoded stim';
    title(sprintf('Decoding accuracy: %s%', num2str(accuracy)))
end

end