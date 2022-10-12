function dat = getBehFigData(db,P,behP)
    
    P.focusOn = behP.focusOn;
    
    %% get pupil part
    
    % for each mouse
    clear dat
    for k = 1:P.nMice
        if P.miceweye(k)
            dattmp = getEyeData(db(k),P);
            
            for dd = 1:numel(dattmp)
                dataA = dattmp(dd).dataA ;
                
                % NB: all 1D
                sigma = nanstd(mat2vec(dataA(:,:,:,1,:)));
                if sigma>0
                    dataA(:,:,:,1,:) = dataA(:,:,:,1,:)/sigma;
                end
                
                dat(dd).dataA = dataA - nanmean(dataA,[2 5]);
                dat(dd).dataMPC(:,:,1,k) = nanmean(dat(dd).dataA,[3 5]);
                
                dat(dd).dataMPC_odd(:,:,1,k) = squeeze(nanmean(dataA(:,:,:,:,1:2:end) ...
                    -nanmean(dataA(:,:,:,:,1:2:end),[2 5]),[3 5]));
                dat(dd).dataMPC_even(:,:,1,k) = squeeze(nanmean(dataA(:,:,:,:,2:2:end) ...
                    -nanmean(dataA(:,:,:,:,1:2:end),[2 5]),[3 5])); % remove the odd average--more similar to what we do to the spikes

                dat(dd).name = dattmp(dd).name;
            end
        end
    end
    
    %% get motion
    
    % for each mouse
    dd = numel(dat)+1;
    for k = 1:P.nMice
        dataA = getOrgData(db(k).motionData(:,behP.pc2Keep,:),db(k).audiovisuoCode,P); % dataA is bins x imp stim x unimp stim x comp x rep

        for pc = 1:numel(behP.pc2Keep)
            si = sign(skewness(mat2vec(dataA(:,:,:,pc,:)))); % reorient it
            sigma = nanstd(mat2vec(dataA(:,:,:,pc,:)));
            if sigma>0
                dataA(:,:,:,pc,:) = dataA(:,:,:,pc,:)/sigma*si;
            end
        end

        dat(dd).dataA = (dataA - nanmean(dataA,[2 5]));
        
        dat(dd).dataMPC(:,:,:,k) = nanmean(dat(dd).dataA,[3 5]);

        dat(dd).dataMPC_odd(:,:,:,k) = squeeze(nanmean(dataA(:,:,:,:,1:2:end) ...
            -nanmean(dataA(:,:,:,:,1:2:end),[2 5]),[3 5]));
        dat(dd).dataMPC_even(:,:,:,k) = squeeze(nanmean(dataA(:,:,:,:,2:2:end) ...
            -nanmean(dataA(:,:,:,:,1:2:end),[2 5]),[3 5])); % remove the odd average--more similar to what we do to the spikes

    end
    dat(dd).name = 'Body motion';