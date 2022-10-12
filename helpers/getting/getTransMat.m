function [M,S,U] = getTransMat(data,P)

% gets transition matrix M
% data should be timecourse x clusters x repetitions
% bins is binning of the stimuli
% P contains diverse general parameters

if ~isfield(P,'reducingMethod')
    P.reducingMethod = 'SVD';
end

if ~isfield(P,'recenterTimePoint')
    P.recenterTimePoint = false;
end

% get a few parameters
s = size(data); s(end+1) = 1;
nTime = s(1);
nBins = P.nBins;
nGroups = nTime/nBins;
nComps = s(2);
nRep = s(3);

switch P.reductionMethod
    case 'DSS'
        D_white = nan(nTime, nComps, nRep);
        for r = 1:nRep
            Z = data(:,:,r);
            
            missing_tps = any(isnan(Z),2);
            
            % reshape and whiten data
            Z = reshape(Z(~missing_tps,:), [sum(~missing_tps), nComps]);
            [U,~,V] = svd(Z-mean(Z), 'econ');
            Z_white = U*V';
            % Z_white = U(:,1:min(1,nComps))*V(:,1:min(1,nComps))'; % reduce dim
            
            D_white(~missing_tps,:,r) = Z_white;
        end
        
        Xm = nanmean(D_white,3); % average over train set
        Xmr = reshape(Xm, [size(Xm,1)/nGroups, nGroups, nComps]);
        Xbias = Xmr((P.bins>0) & (P.bins<3.85),:,:);
        Xbias = reshape(Xbias, [sum((P.bins>0) & (P.bins<3.85))*nGroups, nComps]);
        [U,S,V] = svd(Xbias-mean(Xbias), 'econ');

        % orthogonalize somehow
        [Q,~] = qr(V);

        % find transition matrix
        M = Q;
        
        %%% check if orthonormal?
        
    case 'SVD'
        X = nanmean(data,3);
        missing_tps = isnan(X(:,1));
        X = X(~missing_tps,:);
        [U,S,V] = svd(X-mean(X),'econ');
        M = V;
        
    case 'none'
        M = eye(nComps,nComps);
        S = [];
        U = [];
        
    otherwise
        disp('Sorry, can''t deal with that.')
end