function [C, audiovisuoCode, P, ...
        eventTimes, sp, motSVD, videoTimestampsEphysTime_face, uMotMask, ...
        pupilarea, pupilcom, blink, videoTimestampsEphysTime_eye] = loadONE(dataPath, expRef)

    savePath = fullfile(dataPath,expRef);

    %%

    % trials
    audiovisuoCode(1,:) = readNPY(fullfile(savePath,sprintf('_cb_trials._cb_videoID.%s.npy',expRef)));
    audiovisuoCode(2,:) = readNPY(fullfile(savePath,sprintf('_cb_trials._cb_audioID.%s.npy',expRef)));
    eventTimes = readNPY(fullfile(savePath,sprintf('_cb_trials._cb_onsetTimes.%s.npy',expRef)));

    % spikes
    sp.st = readNPY(fullfile(savePath,sprintf('spikes.times.%s.npy',expRef)));
    sp.clu = readNPY(fullfile(savePath,sprintf('spikes.clusters.%s.npy',expRef)));
    sp.spikeAmps = readNPY(fullfile(savePath,sprintf('spikes.amps.%s.npy',expRef)));
    sp.spikeDepths = readNPY(fullfile(savePath,sprintf('spikes.depths.%s.npy',expRef)));

    % clusters
    C.CluID = readNPY(fullfile(savePath,sprintf('clusters._cb_IDs.%s.npy',expRef)));
    C.CluLab = readNPY(fullfile(savePath,sprintf('clusters._cb_labels.%s.npy',expRef)));
    C.Depth = readNPY(fullfile(savePath,sprintf('clusters.depths.%s.npy',expRef)));
    C.XPos = readNPY(fullfile(savePath,sprintf('clusters._cb_xpos.%s.npy',expRef)));
    C.CluSpknum = readNPY(fullfile(savePath,sprintf('clusters._cb_spkNums.%s.npy',expRef)));
    C.area = readNPY(fullfile(savePath,sprintf('clusters._cb_areas.%s.npy',expRef)));

    % cameras
    % eye
    blink = readNPY(fullfile(savePath,sprintf('eye.blink.%s.npy',expRef)));
    pupilarea = readNPY(fullfile(savePath,sprintf('eye._cb_area.%s.npy',expRef)));
    pupilcom = readNPY(fullfile(savePath,sprintf('eye.centerPos.%s.npy',expRef)));
    videoTimestampsEphysTime_eye = readNPY(fullfile(savePath,sprintf('eye.times.%s.npy',expRef)));
    % face
    motSVD = readNPY(fullfile(savePath,sprintf('camera._cb_motionPCs.%s.npy',expRef)));
    videoTimestampsEphysTime_face = readNPY(fullfile(savePath,sprintf('camera.times.%s.npy',expRef)));
    uMotMask = readNPY(fullfile(savePath,sprintf('_cb_motionPCs.weights.%s.npy',expRef)));
    uMotMask = permute(uMotMask,[2 3 1]);

    % parameters
    load(fullfile(savePath,'parameters.mat'),'P')