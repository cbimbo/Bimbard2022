close all

%% Load preprocessed data

% Choose what dataset to use
% what2Load = 'Normal';
what2Load = 'Transec';

preprocFolder = 'D:\ephys_results\processedData\audioVis\';

preprocAllFile = fullfile(preprocFolder,sprintf('all_ephys%s.mat',what2Load));
if exist(preprocAllFile,'file')
    load(preprocAllFile, 'db_bu', 'P_bu');
else
    % Load raw data
    recompute = 0;
    whichData = sprintf('ephys%s',what2Load);
    [db_bu,P_bu] = loadData_paper(whichData,recompute);

    % Process data (psths, etc.)  
    P_bu.proc.window = [-1,3.8];
    P_bu.proc.binSize = .030; % almost video resolution, so relationship between the two will be a bit jittered at the single bin scale
    P_bu.smooth = 3; % smooth psths with causal filter
    [db_bu,P_bu] = processData_paper(db_bu,P_bu,recompute);

    % Save it
    save(preprocAllFile, 'db_bu', 'P_bu', '-v7.3')
end

%% Subselect cells

% Requires the allenCCF repo (https://github.com/cortex-lab/allenCCF)
% Can bypass it easily by modifying the code.

% Cell subselection parameters--subselect a specific brain area
P.focusOnArea = 'VIS';
% P.focusOnArea = 'HPF';

P_bu.path2CCF = '\\zserver\Lab\Atlas\allenCCF'; % Path to the allen CCF
[db,P,Nc] = subselectCells(db_bu,P_bu,P.focusOnArea,1);
sideTransec = [db.side]; % Side of the transectomy: 1 if ipsi, -1 if contra/no transectomy

%% Define general parameters

% Processing params
P.reductionMethod = 'SVD';
P.demeanTimeCourse = 1;
P.zscore = 1;
P.foldtype = 'oddeven';

%% Generate further processed data

%%% Following sections will generate the data needed to plot each figure.

%% Main figure

mainFigP.focusOn = 'sound';
mainFigP.pc2Keep = 1:4; 
mainFigP.plotfold = 0; % visualize each fold
mainFigP.stimsubset = [27 9 18 25 0];
% mice and cells to plot
switch P.focusOnArea
    case 'VIS'
        % VIS
        switch P.whichData
            case {'ephys','ephysNormal'}
                mainFigP.k2plt = find(cellfun(@(x) ~isempty(x), strfind(P.mouseRef,'AL032')));  % main example mouse is the 1st one
                mainFigP.ii2plt = 100; % main example cell is the 100th one (AL032)
            case 'ephysTransec'
                mainFigP.k2plt = 1;  % main example mouse is the 1st one
                mainFigP.ii2plt = 1; % not really used
        end
    case 'HPF'
        % HPF
        mainFigP.k2plt = find(cellfun(@(x) ~isempty(x), strfind(P.mouseRef,'CB005')));  % main example mouse is the 1st one
        mainFigP.ii2plt = 1; % not really used
end
mainFigP.nGroupsSnd = P.nGroupsSnd;
mainFigP.labelsGroupsSnd = P.labelsGroupsSnd;
mainFigP.SndColors = P.SndColors;
mainFigP.VidColors = P.VidColors;
mainFigP.bins = P.bins;
mainFigP.nMice = P.nMice;
mainFigP.mouseRef = P.mouseRef;

% Generate the files
focusOnL = {'sound','video'}; % Will save different files
for f = 1:numel(focusOnL)
    mainFigP.focusOn = focusOnL{f};
    [dataMPC,M,dataMPC_odd,dataMPC_even,dataA] = getMainFigData(db,P,mainFigP);

    mainFigFile = fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephys%s_%s.mat',mainFigP.focusOn,what2Load,P.focusOnArea));
    mkdir_CB(fileparts(mainFigFile))
    save(mainFigFile,'dataMPC','M','dataMPC_odd','dataMPC_even','dataA','mainFigP','sideTransec')
end

%% Behavior

behFigP = mainFigP;
behFigP.pc2Keep = 1:4; 
behFigP.k2plt = 1;

% Generate the files
focusOnL = {'sound','video'}; % Will save different files
for f = 1:numel(focusOnL)
    behFigP.focusOn = focusOnL{f};
    dat = getBehFigData(db,P,behFigP);

    behFigFile = fullfile(preprocFolder,'Behfig',sprintf('allBehfig_%s_ephys%s_%s.mat',behFigP.focusOn,what2Load,P.focusOnArea));
    mkdir_CB(fileparts(behFigFile))
    save(behFigFile,'dat','behFigP')
end

%% Eigenspectra

compDimP.focusOnList = {'video','sound','interaction'};
compDimP.doitOn = 'spikeData';
compDimP.pcL = 1:max(Nc); % need to have many components for visual
if strcmp(what2Load,'Normal') && strcmp(P.focusOnArea,'VIS')
    compDimP.shuff = 1000;
else
    compDimP.shuff = 0; % no need here
end

% Generate the files
[r,r_shuff] = getEigenSpectra(db,P,compDimP);

eigenspectraFile = fullfile(preprocFolder,'Eigenspectra',sprintf('allEigenspectra_ephys%s_%s.mat',what2Load,P.focusOnArea));
mkdir_CB(fileparts(eigenspectraFile))
save(eigenspectraFile,'r','r_shuff','compDimP')

%% Decoding

decP.focusOnL = {'video','sound'};
decP.numshuf = 0;
switch P.focusOnArea
    case 'VIS'
        decP.pc2KeepSpk = num2cell([4 30]);
    case 'HPF'
        decP.pc2KeepSpk = num2cell(4);
end
decP.pc2KeepMot = num2cell(128);

% Generate the files
acc = getDecoding(db,P,decP);

decodingFile = fullfile(preprocFolder,'Decoding',sprintf('allDecoding_ephys%s_%s.mat',what2Load,P.focusOnArea));
mkdir_CB(fileparts(decodingFile))
save(decodingFile,'acc','decP')

%% Subspace Overlap

subOvP.motionCompNum = 128;
subOvP.lagbins = -ceil(0.1/mean(diff(db(1).videoTimestampsEphysTime_face))):ceil(0.2/mean(diff(db(1).videoTimestampsEphysTime_face)));
subOvP.shuff = 100;
subOvP.meth = 'RRR';
subOvP.nModelComp2keep = 4;

% Generate the files
[Cov_comp_vid,Cov_comp_snd,percExpl_vid,percExpl_snd] = getSubspaceOverlap(db,P,subOvP);

subspaceFile = fullfile(preprocFolder,'SubspaceOverlap',sprintf('allSubspace_ephys%s_%s.mat',what2Load,P.focusOnArea));
mkdir_CB(fileparts(subspaceFile))
save(subspaceFile,'Cov_comp_vid','Cov_comp_snd','percExpl_vid','percExpl_snd','subOvP')

%% Encoding models

% Requires glmnet-matlab (https://github.com/junyangq/glmnet-matlab)

predP.motpc2keepL = 2^7; %2.^(0:8);
predP.motpc2keepfix = 2^7; %%% remove?
predP.model.type = 'ridge';
predP.model.nfold = 3;
predP.model.NL = 0; % whether to add a non-linearity %%% remove?
predP.lagbins = -ceil(0.1/P.proc.binSize):ceil(0.2/P.proc.binSize);
predP.stim = 1;
predP.stimlagbins = 0:sum(P.bins>0)-1;
predP.getQuant = @getCorr; % or @getCov
predP.foldtype = P.foldtype;
predP.tmpSaveFolder =  fullfile(preprocFolder,'Model','tmp');

% Generate the files
focusOnL = {'sound','video'}; % Will save different files
for f = 1:numel(focusOnL)
    predP.focusOn = focusOnL{f};
    switch predP.focusOn
        case 'sound'
            predP.pc2predL = 1:4; % number of neural PCs to predict
        case 'video'
            switch P.focusOnArea
                case 'VIS'
                    predP.pc2predL = 1:30;
                case 'HPF'
                    predP.pc2predL = 1:4;
            end
    end
    predP.model.RRRdimL = 1:max(predP.pc2predL); % list of dimensions for RRR
    [ypred_test, W_all] = getPrediction(db,P,predP);

    predFile = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat',predP.focusOn,predP.model.type,what2Load,P.focusOnArea));
    mkdir_CB(fileparts(predFile))
    save(predFile,'ypred_test', 'W_all', 'predP','-v7.3')
end

% Save the masks to reconstruct the weights in pixel space
facemapWeightsFile = fullfile(preprocFolder,'Model',sprintf('facemapWeights_ephys%s_%s.mat',what2Load,P.focusOnArea));
facemapWeights = cellfun(@(x) x(:,:,1:predP.motpc2keepfix), {db.uMotMask}, 'uni', 0);
save(facemapWeightsFile,'facemapWeights','-v7.3')

%% Lag between behavior and neural activity
lagP.focusOn = 'sound';
lagP.runModel = 1;
lagP.pc2predL = 1; % number of neural PCs to predict
lagP.motpc2keepfix = 2^7; 
lagP.model.type = 'ridge';
lagP.model.nfold = 3;

% Generate the files
lagBehNeurFile = fullfile(preprocFolder,'Lags',sprintf('lagBehNeur_ephys%s_%s.mat',what2Load,P.focusOnArea));
[xc,lags_xc,w,lags_w] = getLags(db,P,lagP);

mkdir_CB(fileparts(lagBehNeurFile))
save(lagBehNeurFile,'xc', 'lags_xc', 'w', 'lags_w')

%% Lag between ipsi and contra--only for the transectomy experiments

if strcmp(what2Load,'Transec')
    ipsiContraP.window = [-1 3.8];
    ipsiContraP.psthBinSize = 0.001;

    % Generate the files
    lagIpsiContraFile = fullfile(preprocFolder,'Lags',sprintf('lagsIpsiContra_ephys%s_%s.mat',what2Load,P.focusOnArea));
    [pc1_ipsi,pc1_contra,lags_xc] = getLagIpsiContra(db,P,ipsiContraP);

    mkdir_CB(fileparts(lagIpsiContraFile))
    save(lagIpsiContraFile,'pc1_ipsi', 'pc1_contra', 'lags_xc')
end