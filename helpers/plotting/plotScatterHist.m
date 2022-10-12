function plotScatterHist(c,P,xedges,yedges,nonewfig)

tmpAll = [];
for  k = 1:P.nMice
    tmpAll = [tmpAll, c{k}];
end

if ~exist('xedges','var')
    l = prctile(tmpAll(1,:),[1,99]);
    xedges = linspace(l(1),l(2),30);
end
if ~exist('yedges','var')
    l = prctile(tmpAll(2,:),[1,99]);
    yedges = linspace(l(1),l(2),30);
end

if ~exist('nonewfig','var')
    figure;
end
hold all
H = hist3(tmpAll','edges',{xedges, yedges});
imagesc(xedges,yedges,H')
set(gca,'YDir','Normal')
colormap(flipud(bone(256)));
axis tight 
% for  k = 1:P.nMice
%     scatter(mean(c{k}(1,:)),mean(c{k}(2,:)),80,P.mouseColor(k,:),'filled')
%     eh = prctile(c{k}(1,:),[10 90]);
%     ev = prctile(c{k}(2,:),[10 90]);
%     errorbar(mean(c{k}(1,:)),mean(c{k}(2,:)), ...
%         abs(ev(1)-mean(c{k}(2,:))), abs(ev(2)-mean(c{k}(2,:))),abs(eh(1)-mean(c{k}(1,:))),abs(eh(2)-mean(c{k}(1,:))),'color','k')
% end
hline(0,'k--')
vline(0,'k--')
plot([yedges(1),yedges(end)],[yedges(1),yedges(end)],'k--')
% set(gca,'visible','off')
colorbar
axis equal tight
end