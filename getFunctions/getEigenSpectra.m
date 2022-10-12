function [r,r_shuff] = getEigenSpectra(db,P,compDimP)
    
    focusOnList = compDimP.focusOnList;
    doitOn = compDimP.doitOn;
    pcL = compDimP.pcL;
    shuff = compDimP.shuff;
    
    r.Cov_comp = nan(numel(pcL),2,P.nMice,numel(focusOnList));
    r_shuff.Cov_comp = nan(numel(pcL),2,P.nMice,numel(focusOnList),shuff);
    
    for f = 1:numel(focusOnList)
        
        switch focusOnList{f}
            case {'sound','video'}
                whicDimMeasure = @getCompCov;
                P.focusOn = focusOnList{f};
            case {'interaction'}
                whicDimMeasure = @getCompCov_interactions;
                P.focusOn = 'video';
        end

        fprintf('*** Focusing on %s... ***\n', focusOnList{f})
        for  k = 1:P.nMice
            fprintf('Computing for mouse %s...', P.mouseRef{k})
            
            P.dim.pcL = pcL(pcL<size(db(k).(doitOn),2));
            
            % real covariance
            tmp = whicDimMeasure(db(k).(doitOn)(:,:,:,1),db(k).audiovisuoCode,P);
            r.Cov_comp(1:size(tmp,1),1:size(tmp,2),k,f) = tmp;
            
            % shuffled
            if shuff > 0
                fprintf('Shuffling...\n')
                for sh = 1:shuff
                    code = [];
                    for rep = 1:4
                        codetmp = db(k).audiovisuoCode(:,1:P.nGroupsSnd*P.nGroupsVid); % same over repeats
                        for g = 1:12 % shuffle imp
                            codetmp(strcmp(P.focusOn,{'video','sound'}),codetmp(~strcmp(P.focusOn,{'video','sound'}),:) == P.labelsGroupsSnd(g)) = ...
                                P.labelsGroupsSnd(randperm(12));
                        end
                        for g = 1:12 % shuffle unimp
                            codetmp(~strcmp(P.focusOn,{'video','sound'}),codetmp(strcmp(P.focusOn,{'video','sound'}),:) == P.labelsGroupsSnd(g)) = ...
                                P.labelsGroupsSnd(randperm(12));
                        end
                        code = [code, codetmp];
                    end
                    
                    code = code(:,1:size(db(k).audiovisuoCode,2));
                    
                    tmp = whicDimMeasure(db(k).(doitOn),code,P);
                    r_shuff.Cov_comp(1:size(tmp,1),1:size(tmp,2),k,f,sh) = tmp;
                end
                fprintf('Done shuffling.\n')
            end
            fprintf('Done.\n')
        end
        fprintf('*** Focusing on %s done. ***\n', focusOnList{f})
    end

end