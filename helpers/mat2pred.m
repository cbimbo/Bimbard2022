function pred = mat2pred(mat,lag,paddingval)

% lag will padd on each side
if ~exist('lag','var')
    lag = 0;
end
    
pred = permute(mat,[1 2 3 5 4]);
s = size(pred); s(end+1) = 1;

if lag>0
   pred = cat(1,nan(lag,s(2),s(3),s(4),s(5)),pred);
   pred = cat(1,pred,nan(lag,s(2),s(3),s(4),s(5)));
   s(1) = size(pred,1);
end

pred = reshape(pred,[s(1)*s(2)*s(3)*s(4),s(5)]);