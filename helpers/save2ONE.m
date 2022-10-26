function save2ONE(dataPath,expRef,dat)
    %%% This function will save the data in the ONE format.

    savePath = fullfile(dataPath,expRef);
    mkdir(savePath)

    %%

    % trials
    saveONEFormat(dat.audiovisuoCode(1,:)',savePath,'_cb_trials','_cb_videoID','npy',expRef);
    saveONEFormat(dat.audiovisuoCode(2,:)',savePath,'_cb_trials','_cb_audioID','npy',expRef);
    saveONEFormat(dat.eventTimes,savePath,'_cb_trials','_cb_onsetTimes','npy',expRef);

    % spikes
    saveONEFormat(dat.sp.st,savePath,'spikes','times','npy',expRef);
    saveONEFormat(dat.sp.clu,savePath,'spikes','clusters','npy',expRef);
    saveONEFormat(dat.sp.spikeAmps,savePath,'spikes','amps','npy',expRef);
    saveONEFormat(dat.sp.spikeDepths,savePath,'spikes','depths','npy',expRef);

    % clusters
    saveONEFormat(dat.C.CluID',savePath,'clusters','_cb_IDs','npy',expRef);
    saveONEFormat(dat.C.CluLab',savePath,'clusters','_cb_labels','npy',expRef);
    saveONEFormat(dat.C.Depth',savePath,'clusters','depths','npy',expRef);
    saveONEFormat(dat.C.XPos',savePath,'clusters','_cb_xpos','npy',expRef);
    saveONEFormat(dat.C.CluSpknum',savePath,'clusters','_cb_spkNums','npy',expRef);
    saveONEFormat(dat.C.area',savePath,'clusters','_cb_areas','npy',expRef);

    % cameras
    % eye
    saveONEFormat(dat.blink,savePath,'eye','blink','npy',expRef);
    saveONEFormat(dat.pupilarea',savePath,'eye','_cb_area','npy',expRef);
    saveONEFormat(dat.pupilcom,savePath,'eye','centerPos','npy',expRef);
    saveONEFormat(dat.videoTimestampsEphysTime_eye,savePath,'eye','times','npy',expRef);
    % face
    saveONEFormat(dat.motSVD{2},savePath,'camera','_cb_motionPCs','npy',expRef);
    saveONEFormat(dat.videoTimestampsEphysTime_face,savePath,'camera','times','npy',expRef);
    saveONEFormat(permute(dat.uMotMask{2},[3 1 2]),savePath,'_cb_motionPCs','weights','npy',expRef);

    % parameters
    P = dat.P;
    save(fullfile(savePath,'parameters.m'),'P')
end