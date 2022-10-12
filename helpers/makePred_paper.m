function [w,ytest,ytrain,ytest_pred,ytrain_pred,RRRdim] = makePred_paper(y,X,idxfold,model)

% y is the data you want to predict (timepoint x components)
% X is the set of predictors you want to predict with (timepoint x
% predictors)
% idxfold are the timepoints you want to leave out in the test set
% model is a structure that contains what prediction model to use, etc.

if ~exist('model','var')
    model.type = 'ols';
end

if exist('idxfold','var') && ~isempty(idxfold)
    ytest = y(idxfold,:);
    ytrain = y(~idxfold,:);
    
    Xtest = X(idxfold,:);
    Xtrain = X(~idxfold,:);
else
    ytest = y;
    ytrain = y;
    
    Xtest = X;
    Xtrain = X;
end

idxnanX = any(isnan(Xtrain),2);
idxnany = any(isnan(ytrain),2);
idxnan = idxnanX | idxnany;
Xtrain(idxnan,:) = [];
ytrain(idxnan,:) = [];

switch model.type
    case {'ridge','lasso'}
        % cross-validated model
        switch model.type
            case 'ridge'
                options.alpha = 0;
            case 'lasso'
                options.alpha = 1;
        end
        options.nlambda = 30;
        
        w = nan(size(Xtrain,2)+1,size(ytrain,2));
        ytrain_pred = nan(size(ytrain));
        ytest_pred = nan(size(ytest));
        Xtrain = double(Xtrain); % otherwise crashes pretty badly..?
        ytrain = double(ytrain);
        tic
        parfor pc = 1:size(ytrain,2)
            CVerr = cvglmnet(Xtrain,ytrain(:,pc),'gaussian',options,'mse',model.nfold);
            Beta = CVerr.glmnet_fit.beta(:,find(CVerr.glmnet_fit.lambda == CVerr.lambda_min));
            Bias = CVerr.glmnet_fit.a0(find(CVerr.glmnet_fit.lambda == CVerr.lambda_min));

            w(:,pc) = [Beta;Bias];
            
            % plot prediction
            ytrain_pred(:,pc) = Xtrain*Beta + Bias;
            ytest_pred(:,pc) = Xtest*Beta + Bias;
        end
        toc
        
    case 'ols'
        Xtrain = cat(2,Xtrain,ones(size(Xtrain,1),1));
        Xtest = cat(2,Xtest,ones(size(Xtest,1),1));
        
        w = pinv(Xtrain'*Xtrain)*Xtrain'*ytrain;
        
        ytrain_pred = Xtrain*w;
        ytest_pred = Xtest*w;
        
    case {'rrr','RRR'}
        Xtrain = cat(2,Xtrain,ones(size(Xtrain,1),1));
        Xtest = cat(2,Xtest,ones(size(Xtest,1),1));

        % loop over folds to find rank
        nfold = model.nfold; % should be between 1 and numel(unimportant stim)?
        foldsize = floor(size(Xtrain,1)/nfold);
        ev_test_tot = nan(numel(model.RRRdimL),1);
        ev_train_tot = nan(numel(model.RRRdimL),1);

        for nf = 1:nfold
            idxfoldrrr = zeros(1,size(Xtrain,1));
            idxfoldrrr((nf-1)*foldsize+(1:foldsize)) = 1; % do it pretty abruptly
            idxfoldrrr = logical(idxfoldrrr);

            % get train and test within train
            Xtraintrain = Xtrain(~idxfoldrrr,:);
            ytraintrain = ytrain(~idxfoldrrr,:);
            Xtraintest = Xtrain(idxfoldrrr,:);
            ytraintest = ytrain(idxfoldrrr,:);

            bols = pinv(Xtraintrain'*Xtraintrain)*Xtraintrain'*ytraintrain;
            
            ytraintrain_pred = Xtraintrain*bols;
            
            [~,~,v] = svd(ytraintrain_pred,'econ');

            for pc = 1:numel(model.RRRdimL)
                pc2keep = model.RRRdimL(pc);
                brrr = bols*v(:,1:pc2keep)*v(:,1:pc2keep)';
                
                ytraintest_pred = Xtraintest*brrr;
                ytraintrain_pred = Xtraintrain*brrr;
                
                % get goodness of fit
                ev_test_tot(pc,nf) = getEV(ytraintest(:),ytraintest_pred(:));
                ev_train_tot(pc,nf) = getEV(ytraintrain(:),ytraintrain_pred(:));

            end
        end

        % choose best hyperparameter 
        [~,RRRdim] = max(mean(ev_test_tot,2));
    
        % recompute with all folds
        bols = pinv(Xtrain'*Xtrain)*Xtrain'*ytrain;
        ytrain_pred = Xtrain*bols;
        [~,~,v] = svd(ytrain_pred,'econ');
        pc2keep = model.RRRdimL(RRRdim);
        brrr = bols*v(:,1:pc2keep)*v(:,1:pc2keep)';
        % save 
        ytest_pred = Xtest*brrr;
        ytrain_pred = Xtrain*brrr;
        w = brrr;

    otherwise
        disp(['Method ''' model.type ''' unknown'])
end

end