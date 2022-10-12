function X = getXpred(x,lagbins,pad)

% will shape X predictor matrix with certain lags
% pad = 1 if need to pad with last/bext value, nan if want to pad w/ nan

if nargin<3 
    pad = nan;
end

nT = size(x,1);
nPC = size(x,2);

X = nan(nT,numel(lagbins)*nPC);
for l = 1:numel(lagbins)
    if lagbins(l) < 0
        pred = cat(1, x(-lagbins(l)+1:end,:), pad*ones(-lagbins(l),nPC).*x(-lagbins(l)+1,:));
    elseif lagbins(l) > 0
        pred = cat(1, pad*ones(lagbins(l),nPC).*x(1,:), x(1:end-lagbins(l),:));
    else
        pred = x;
    end
    X(:,l:numel(lagbins):(numel(lagbins))*(nPC-1)+l) = pred;
end

end