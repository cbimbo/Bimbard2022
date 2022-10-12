function [nfold, folds] = getFolds(s,foldtype)

nBins = s(1);
nImp = s(2);
nUnImp = s(3);
nClu = s(4);
nRep = s(5);
    
switch foldtype
    case 'oddeven'
        
        nfold = 2;
        for ff = 1:nfold
            fold = zeros(nBins,nImp,nUnImp,1,nRep);
            fold(:,:,:,:,ff:2:end) = 1;
            %     idxfold(:,:,ff,:,:) = 1;
            fold = logical(repmat(fold,[1 1 1 nClu 1]));
            folds{ff} = fold;
        end
        
    case 'oneunimp'
        
        nfold = nUnImp;
        for ff = 1:nfold
            fold = zeros(nBins,nImp,nUnImp,1,nRep);
            fold(:,:,ff,:,:) = 1;
            fold = logical(repmat(fold,[1 1 1 nClu 1]));
            folds{ff} = fold;
        end
        
    case 'onerep'
    
        nfold = nRep;
        for ff = 1:nfold
            fold = zeros(nBins,nImp,nUnImp,1,nRep);
            fold(:,:,:,:,ff) = 1;
            fold = logical(repmat(fold,[1 1 1 nClu 1]));
            folds{ff} = fold;
        end
end