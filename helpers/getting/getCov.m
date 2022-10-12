function c = getCov(y,ypred)

idxnan_y = any(isnan(y),2);
idxnan_ypred = any(isnan(ypred),2);
idxnan = idxnan_y | idxnan_ypred;
y(idxnan,:) = [];
ypred(idxnan,:) = [];

c = sum((y-mean(y)).*(ypred-mean(y)))/(numel(y)-1);

end
