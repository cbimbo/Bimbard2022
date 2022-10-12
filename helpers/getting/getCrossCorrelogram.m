function xc = getCrossCorrelogram(dat1,dat2,lags,noResign)

if ~exist('noResign','var')
    noResign = 0;
end

% get cross-correlogram
nMice = numel(dat1);
xc = cell([1,nMice]);
for k = 1:nMice
    lagbins = int64(lags{k}/(lags{k}(2)-lags{k}(1)));
    xc_tmp = nan(numel(lagbins),1);
    sTot = size(dat1{k},1);
    for l = 1:numel(lagbins)
        if lagbins(l)<0
            idx1 = 1:sTot+lagbins(l);
            idx2 = -lagbins(l)+1:sTot;
        else
            idx1 = lagbins(l)+1:sTot;
            idx2 = 1:sTot-lagbins(l);
        end
        idxnan = isnan(dat1{k}(idx1)) | isnan(dat2{k}(idx2));
        xc_tmp(l) = corr(dat1{k}(idx1(~idxnan)),dat2{k}(idx2(~idxnan)));
    end
    xc{k} = xc_tmp;
end

% resign it--sign is initially arbitrary
if ~noResign
    for k = 1:nMice
        xc{k} = xc{k}*sign(xc{k}(find(lags{k}>0,1)));
    end
end
