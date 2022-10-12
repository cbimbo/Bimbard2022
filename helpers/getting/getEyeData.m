function dat = getEyeData(db,P)
    
    n = 1;
    Blink = db.blinkData;
    Blink(:) = nanzscore(Blink(:));
    threshBlink = double(Blink<-10); % arbitrary here but seems to fit quite well when manually checking the videos
    if all(isnan(threshBlink(:)))
        threshBlink(:) = 0;
    end
    threshBlink(isnan(Blink)) = nan;
    [dat(n).dataA, dat(n).dataM] = getOrgData(threshBlink,db.audiovisuoCode,P); % get blink data
    dat(n).name = 'Blink';
    n = n+1;
    
    [dat(n).dataA, dat(n).dataM] = getOrgData(db.pupilareaData,db.audiovisuoCode,P); % get pupil area data
    dat(n).name = 'Pupil area';
    n = n+1;
    
    [dataApc, dataMpc] = getOrgData(db.pupilcomData,db.audiovisuoCode,P); % get pupil com data
    dat(n).dataA = dataApc(:,:,:,1,:);
    dat(n).dataM = dataMpc(:,:,:,1);
    dat(n).name = 'Pupil position (X)';
    n = n+1;
    dat(n).dataA = dataApc(:,:,:,2,:);
    dat(n).dataM = dataMpc(:,:,:,2);
    dat(n).name = 'Pupil position (Y)';
    n = n+1;
    
    [dataApm, dataMpm] = getOrgData(db.pupilmotData,db.audiovisuoCode,P); % get motion of pupil data
    dat(n).dataA = dataApm(:,:,:,1,:);
    dat(n).dataM = dataMpm(:,:,:,1);
    dat(n).name = 'Pupil motion (X, sign.)';
    n = n+1;
    dat(n).dataA = dataApm(:,:,:,2,:);
    dat(n).dataM = dataMpm(:,:,:,2);
    dat(n).name = 'Pupil motion (Y, sign.)';
    n = n+1;
    
    dat(n).dataA = abs(dataApm(:,:,:,1,:));
    dat(n).dataM = nanmean(dat(n).dataA,5);;
    dat(n).name = 'Pupil motion (X)';
    n = n+1;
    dat(n).dataA = abs(dataApm(:,:,:,2,:));
    dat(n).dataM = nanmean(dat(n).dataA,5);
    dat(n).name = 'Pupil motion (Y)';
    n = n+1;
    
    tmp = sqrt(db.pupilmotData(:,1,:).^2 + db.pupilmotData(:,2,:).^2);
    [dat(n).dataA, dat(n).dataM] = getOrgData(tmp,db.audiovisuoCode,P); % get distance of pupil data to center
    dat(n).name = 'Pupil motion';
    n = n+1;
end