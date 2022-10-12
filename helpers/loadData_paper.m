function [db,P] = loadData_paper(whichData,recompute)
    
    fprintf('*** Loading data... ***\n')
    if contains(whichData,'ephys')
        n = 1;        
        switch whichData
            case {'ephys', 'ephysNormal'}
                % V1/HPF
                normalColor = hot(16);
                db(n).subject = 'CB007'; db(n).color = normalColor(1,:); db(n).date = '2021-02-19'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2020-09-15'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'AL030'; db(n).color = normalColor(2,:); db(n).date = '2020-02-19'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2019-06-26'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'AL031'; db(n).color = normalColor(3,:); db(n).date = '2020-02-19'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2019-06-03'; db(n).sex = 'F'; n = n+1;
                db(n).subject = 'CB005'; db(n).color = normalColor(4,:); db(n).date = '2020-12-07'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2020-08-11'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CB008'; db(n).color = normalColor(5,:); db(n).date = '2021-03-10'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2020-09-15'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CB003'; db(n).color = normalColor(6,:); db(n).date = '2020-11-11'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2020-06-02'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'AL036'; db(n).color = normalColor(7,:); db(n).date = '2020-02-25'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2019-10-17'; db(n).sex = 'F'; n = n+1;
                db(n).subject = 'AL032'; db(n).color = normalColor(8,:); db(n).date = '2020-03-05'; db(n).imecX = []; db(n).side = -1; db(n).birthdate = '2019-07-24'; db(n).sex = 'M'; n = n+1;
                
            case 'ephysTransec'
                tmp1 = winter(4); tmp2 = autumn(3); tmp3 = spring(6);
                transColor = [tmp1(1:2,:); tmp2(2,:); tmp3(1:3,:)];
                db(n).subject = 'CR_transectomy3'; db(n).color = transColor(1,:); db(n).date = '2021-10-26'; db(n).transecdate = '2021-10-12'; db(n).imecX = 1; db(n).side = -1; db(n).birthdate = '2021-08-18'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy3'; db(n).color = transColor(1,:); db(n).date = '2021-10-26'; db(n).transecdate = '2021-10-12'; db(n).imecX = 0; db(n).side = 1; db(n).birthdate = '2021-08-18'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy3'; db(n).color = transColor(2,:); db(n).date = '2021-10-27'; db(n).transecdate = '2021-10-12'; db(n).imecX = 0; db(n).side = -1; db(n).birthdate = '2021-08-18'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy3'; db(n).color = transColor(2,:); db(n).date = '2021-10-27'; db(n).transecdate = '2021-10-12'; db(n).imecX = 1; db(n).side = 1; db(n).birthdate = '2021-08-18'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy7'; db(n).color = transColor(3,:); db(n).date = '2021-11-25'; db(n).transecdate = '2021-11-18'; db(n).imecX = 0; db(n).side = -1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy7'; db(n).color = transColor(3,:); db(n).date = '2021-11-25'; db(n).transecdate = '2021-11-18'; db(n).imecX = 1; db(n).side = 1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(4,:); db(n).date = '2021-12-02'; db(n).transecdate = '2021-10-12'; db(n).imecX = 0; db(n).side = -1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(4,:); db(n).date = '2021-12-02'; db(n).transecdate = '2021-10-12'; db(n).imecX = 1; db(n).side = 1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(5,:); db(n).date = '2021-12-03'; db(n).transecdate = '2021-10-12'; db(n).imecX = 0; db(n).side = -1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(5,:); db(n).date = '2021-12-03'; db(n).transecdate = '2021-10-12'; db(n).imecX = 1; db(n).side = 1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(6,:); db(n).date = '2021-12-04'; db(n).transecdate = '2021-10-12'; db(n).imecX = 0; db(n).side = -1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
                db(n).subject = 'CR_transectomy8'; db(n).color = transColor(6,:); db(n).date = '2021-12-04'; db(n).transecdate = '2021-10-12'; db(n).imecX = 1; db(n).side = 1; db(n).birthdate = '2021-06-28'; db(n).sex = 'M'; n = n+1;
        end
        
        % put it in subparts of P?
        P.probeNumber = 1; % fixed?
        P.pathSaving = 'D:\ephys_results'; % fixed?
        P.nMice = numel(db);
        P.exp = 'audioVis';
        P.whichData = whichData;

        for k = 1:P.nMice
            fprintf('Loading mouse %s...',[db(k).subject '_' db(k).date])
            [db(k).C, db(k).audiovisuoCode, db(k).trialNumPerStim, Ptmp, ...
                db(k).eventTimes, db(k).sp, ...
                db(k).motSVD, db(k).videoTimestampsEphysTime_face, db(k).uMotMask, ...
                db(k).pupilarea, db(k).pupilcom, db(k).pupilmot, db(k).blink, db(k).videoTimestampsEphysTime_eye] = ...
                loadDataEphys_paper(db(k).subject, db(k).date, db(k).imecX, P, recompute);
            
            % save individual mouse info / not super elegant but easier
            % to have P for each mouse in general
            if k == 1
                P = Ptmp;
                P = rmfield(P,'fileref');
            else
                %%% Should add a check that it's the same? Very unlikely
                %%% that it isn't.
            end
            P.nMice = numel(db);
            P.mouseRef{k} = [db(k).subject '_' db(k).date];
            P.fileref.blockDir{k} = regexprep(Ptmp.fileref.blockDir,'128.40.224.65','zinu'); % in case temp server hasn't been changed
            P.fileref.rootDir{k} = regexprep(Ptmp.fileref.rootDir,'128.40.224.65','zinu');
            P.fileref.ksDir{k} = regexprep(Ptmp.fileref.ksDir,'128.40.224.65','zinu');
            if ~isempty(Ptmp.fileref.anatDir)
                P.fileref.anatDir{k} = regexprep(Ptmp.fileref.anatDir,'128.40.224.65','zinu');
            else
                P.fileref.anatDir{k} = [];
            end
            P.miceweye(k) = Ptmp.miceweye;
            fprintf('Done.\n')
        end
        
        % get mouse colors
        % P.mouseColor = hot(numel(db)*2); P.mouseColor = P.mouseColor(1:numel(db),:);
        for k = 1:P.nMice
            P.mouseColor(k,:) = db(k).color;
        end
        
    else
        % no other cases for now
    end
    
    
    %% get colors for the stims used... should maybe be inside the loop?
    P.SndColors = zeros(P.nGroupsSnd,3);
    P.SndColors(:,1) = linspace(0.2, 1, P.nGroupsSnd);
    P.SndColors(:,3) = linspace(0.2, 1, P.nGroupsSnd);
    P.SndColors(1,:) = [0.5,0.5,0.5];
    P.labelsGroupsSnd = [0 27 9 18 26 17 15 11 29 22 7 25];
    P.VidColors = zeros(P.nGroupsVid,3);
    P.VidColors(:,2) = linspace(0.2, 1, P.nGroupsVid);
    P.VidColors(P.labelsGroupsVid == 0,:) = [0.5,0.5,0.5];
    P.labelsGroupsVid = [0 27 9 18 26 17 15 11 29 22 7 25];
    fprintf('*** Loading done. ***\n')
end