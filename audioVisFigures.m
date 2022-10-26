%% Define general file names

preprocFolder = 'D:\ephys_results\processedData\audioVis\';

% Main Figures
% Fig. 1
mainFigFile_fig1 = fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephys%s_%s.mat','sound','Normal','VIS'));
eigenspectraFile_fig1 = fullfile(preprocFolder,'Eigenspectra',sprintf('allEigenspectra_ephys%s_%s.mat','Normal','VIS'));
decodingFile_fig1 = fullfile(preprocFolder,'Decoding',sprintf('allDecoding_ephys%s_%s.mat','Normal','VIS'));

% Fig. 2
mainFigFile_fig2 = fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephys%s_%s.mat','sound','Normal','HPF'));
eigenspectraFile_fig2 = fullfile(preprocFolder,'Eigenspectra',sprintf('allEigenspectra_ephys%s_%s.mat','Normal','HPF'));
decodingFile_fig2 = fullfile(preprocFolder,'Decoding',sprintf('allDecoding_ephys%s_%s.mat','Normal','HPF'));

% Fig. 3
mainFigFile_fig3 = fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephys%s_%s.mat','sound','Transec','VIS'));
eigenspectraFile_fig3 = fullfile(preprocFolder,'Eigenspectra',sprintf('allEigenspectra_ephys%s_%s.mat','Transec','VIS'));
decodingFile_fig3 = fullfile(preprocFolder,'Decoding',sprintf('allDecoding_ephys%s_%s.mat','Transec','VIS'));

% Fig. 4
behFigFile_fig4 = fullfile(preprocFolder,'Behfig',sprintf('allBehfig_%s_ephys%s_%s.mat','sound','Normal','VIS'));
lagBehNeurFile_fig4 = fullfile(preprocFolder,'Lags',sprintf('lagBehNeur_ephys%s_%s.mat','Normal','VIS'));
subspaceFile_fig4 = fullfile(preprocFolder,'SubspaceOverlap',sprintf('allSubspace_ephys%s_%s.mat','Normal','VIS'));
predFile_fig4 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','sound','ridge','Normal','VIS'));

% Extended Data Figures

% Extended Data Fig. 6 (also used in 10)
lagBehNeurFile_hpf_EDfig6 = fullfile(preprocFolder,'Lags',sprintf('lagBehNeur_ephys%s_%s.mat','Normal','HPF'));
lagBehNeurFile_transec_EDfig6 = fullfile(preprocFolder,'Lags',sprintf('lagBehNeur_ephys%s_%s.mat','Transec','VIS'));
subspaceFile_hpf_EDfig6 = fullfile(preprocFolder,'SubspaceOverlap',sprintf('allSubspace_ephys%s_%s.mat','Normal','HPF'));
subspaceFile_transec_EDfig6 = fullfile(preprocFolder,'SubspaceOverlap',sprintf('allSubspace_ephys%s_%s.mat','Transec','VIS'));

% Extended Data 7
facemapWeightsFile_EDfig7 = fullfile(preprocFolder,'Model',sprintf('facemapWeights_ephys%s_%s.mat','Normal','VIS'));

% Extended Data 8
predFile_video_NormalVIS_EDfig8 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','video','ridge','Normal','VIS'));
predFile_sound_NormalHPF_EDfig8 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','sound','ridge','Normal','HPF'));
predFile_video_NormalHPF_EDfig8 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','video','ridge','Normal','HPF'));
predFile_sound_TransecVIS_EDfig8 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','sound','ridge','Transec','VIS'));
predFile_video_TransecVIS_EDfig8 = fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-%s_ephys%s_%s.mat','video','ridge','Transec','VIS'));

% Transectomy sides--mainly for plotting reasons
side_NormalVIS = load(mainFigFile_fig1,'sideTransec'); side_NormalVIS = side_NormalVIS.sideTransec;
side_NormalHPF = load(mainFigFile_fig2,'sideTransec'); side_NormalHPF = side_NormalHPF.sideTransec;
side_TransecVIS = load(mainFigFile_fig3,'sideTransec'); side_TransecVIS = side_TransecVIS.sideTransec;

%% Fig. 1 or 2--Auditory responses in V1 or HPF

% Requires the matlab version of 'rastermap' (https://github.com/MouseLand/rastermap)

figNum = 1; % 1 or 2

% Load data
switch figNum
    case 1
        eigDat = load(eigenspectraFile_fig1);
        mainFigDat = load(mainFigFile_fig1);
        decDat = load(decodingFile_fig1);
    case 2
        eigDat = load(eigenspectraFile_fig2);
        mainFigDat = load(mainFigFile_fig2);
        decDat = load(decodingFile_fig2);
end

% Main figure
[c_audPC1Corr,p_audPC1Corr] = plotMainFig(mainFigDat,eigDat);

% Decoding
[accVid,accSnd,p_dec] = plotDecoding(decDat.acc.spk,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0);
% close up
plotDecoding(decDat.acc.spk,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0); 
xlim([80 100]); ylim([5 25]);

% --- Quantif ---
% Decoding
fprintf('Video decoding accuracy: %d +- %d (p-value: %d).\n', nanmean(accVid),nansem(accVid),p_dec(1))
fprintf('Sound decoding accuracy: %d +- %d (p-value: %d).\n', nanmean(accSnd),nansem(accSnd),p_dec(2))

% Percentage of sound-related variance explained by audPC1
audPC1_sndCov = squeeze(sum(nanmean(eigDat.r.Cov_comp(1,:,:,2),2),1)./nansum(nanmean(eigDat.r.Cov_comp(:,:,:,2),2), [1 4] ))*100;
fprintf('Percentage sound-related variance explained by audPC1: %d +- %d %% (example mouse: %d).\n', ...
    nanmean(audPC1_sndCov),nansem(audPC1_sndCov),audPC1_sndCov(mainFigDat.mainFigP.k2plt))

% Percentage of total variance explained by visPC1
visPC1_totCov = squeeze(sum(nanmean(eigDat.r.Cov_comp(1,:,:,1),2),1)./nansum(nanmean(eigDat.r.Cov_comp,2), [1 4] ))*100;
fprintf('Percentage total variance explained by visPC1: %d +- %d %% (example mouse: %d).\n', ...
    nanmean(visPC1_totCov),nansem(visPC1_totCov),visPC1_totCov(mainFigDat.mainFigP.k2plt))

% Percentage of total variance explained by audPC1
audPC1_totCov = squeeze(sum(nanmean(eigDat.r.Cov_comp(1,:,:,2),2),1)./nansum(nanmean(eigDat.r.Cov_comp,2), [1 4] ))*100;
fprintf('Percentage total variance explained by audPC1: %d +- %d %% (example mouse: %d).\n', ...
    nanmean(audPC1_totCov),nansem(audPC1_totCov),audPC1_totCov(mainFigDat.mainFigP.k2plt))

% Bias of audPC1 weights
bias_audPC1 = cell2mat(cellfun(@(x) mean(x(:,1)), mainFigDat.M, 'uni', 0));
bias_audPC1_p = signrank(bias_audPC1);
fprintf('Bias of weight distribution: %d +- %d (p-value: %d).\n', ...
    nanmean(bias_audPC1),nansem(bias_audPC1),bias_audPC1_p)

% Correlations across and within mice
s = size(mainFigDat.dataMPC_even);
c = corr(reshape(mainFigDat.dataMPC_even(:,:,1,:),[s(1)*s(2),s(4)]), ...
        reshape(mainFigDat.dataMPC_odd(:,:,1,:),[s(1)*s(2),s(4)]));
c_across = ztransform(c(triu(true(size(c)),1)),'mean'); % Z-transform
c_within = ztransform(diag(c),'mean');
fprintf('Correlation of audPC1 across mice: %d (test-retest: %d).\n',c_across,c_within)

% Correlation with Fig1
fprintf('Correlation with VIS audPC1 timecourse: %d (p-value: %d).\n',c_audPC1Corr,p_audPC1Corr)

% Number of components that are significant
if isfield(eigDat,'r_shuff')
    thre = 1; % 99% confidence interval
    dim = nan(numel(eigDat.compDimP.focusOnList),size(eigDat.r.Cov_comp,3));
    for f = 1:numel(eigDat.compDimP.focusOnList)
        for k = 1:size(eigDat.r.Cov_comp,3)
            uplim = prctile(squeeze(nanmean(eigDat.r_shuff.Cov_comp(:,:,k,strcmp(eigDat.compDimP.focusOnList,eigDat.compDimP.focusOnList{f}),:),2))',100-thre)';
            d = eigDat.compDimP.pcL(find(nanmean(eigDat.r.Cov_comp(:,:,k,strcmp(eigDat.compDimP.focusOnList,eigDat.compDimP.focusOnList{f})),2)-uplim<0,1))-1;
            if isempty(d)
                d = eigDat.compDimP.pcL(find(isnan(uplim),1)); % lower bound
            end
            dim(f,k) = d;
        end
    end
    fprintf('Number of visual components: %d +- %d.\n',nanmean(dim(1,:)),nansem(dim(1,:)))
    fprintf('Number of auditory components: %d +- %d.\n',nanmean(dim(2,:)),nansem(dim(2,:)))
    fprintf('Number of interaction components: %d +- %d.\n',nanmean(dim(3,:)),nansem(dim(3,:)))
end

%% Fig. 3--Transectomy experiments

% Load data
mainFigDatTransec = load(mainFigFile_fig3);
eigDatTransec = load(eigenspectraFile_fig3);
decDatTransec = load(decodingFile_fig3);
eigDatNormal = load(eigenspectraFile_fig1);
decDatNormal = load(decodingFile_fig1);
transecAnatFolder = 'C:\Users\Hamish\OneDrive - University College London\Documents\GitHub\Bimbard2022\transecAnat';

[c_ContravsIpsi,p_ContravsIpsi,c_Contra,p_Contra,c_Ipsi,p_Ipsi] = plotTransecFig(mainFigDatTransec,eigDatTransec);

figure('Position',[595   685   763   293]);
% Auditory input
subplot(151)
audIn = load(fullfile(transecAnatFolder,'TransectomyAnat','auditoryInput.mat'));
maxIn = nan(1,numel(audIn.A));
audInContra = nan(1,numel(audIn.A));
audInIpsi = nan(1,numel(audIn.A));
for ani = 1:numel(audIn.A)
    maxIn(ani) = nansum(audIn.A(ani).fibers2IpsiVIS_tot +  audIn.A(ani).fibers2ContraVIS_tot);
    audInContra(ani) = nansum(audIn.A(ani).fibers2IpsiVIS_tot-audIn.A(ani).fibers2ContraVIS_cut_fromContra + audIn.A(ani).fibers2ContraVIS_tot-audIn.A(ani).fibers2ContraVIS_cut);
    audInIpsi(ani) = nansum(audIn.A(ani).fibers2IpsiVIS_tot-audIn.A(ani).fibers2IpsiVIS_cut + audIn.A(ani).fibers2ContraVIS_tot-audIn.A(ani).fibers2IpsiVIS_cut_fromContra);
end
plotVersus(audInContra./maxIn,audInIpsi./maxIn,-ones(1,numel(audInContra)),[],0,[0 1],1)
xlabel('To contra VIS')
ylabel('To ipsi VIS');
scatter(1,1,30,[9 133 140]/256,'s','filled')

% Sound-related variance explained
subplot(152)
audPC1_sndCov = squeeze(sum(nanmean(eigDatTransec.r.Cov_comp(1:4,:,:,2),2),1)./nansum(nanmean(eigDatTransec.r.Cov_comp,2), [1 4] ))*100;
snd_cov_Normal = nanmean(squeeze(sum(nanmean(eigDatNormal.r.Cov_comp(1:4,:,:,2),2),1)./nansum(nanmean(eigDatNormal.r.Cov_comp,2), [1 4] ))*100);
plotVersus(audPC1_sndCov(mainFigDatTransec.sideTransec == -1),audPC1_sndCov(mainFigDatTransec.sideTransec == 1),-ones(1,sum(mainFigDatTransec.sideTransec == -1)),[],0,[0 10],1)
scatter(snd_cov_Normal,snd_cov_Normal,30,[9 133 140]/256,'s','filled')
xlabel('In contra VIS')
ylabel('In ipsi VIS');

% Decoding
% Video
subplot(153)
cmp2keepneur = 30;
accVid_ipsi = decDatTransec.acc.spk(strcmp(decDatTransec.decP.focusOnL, 'video'),mainFigDatTransec.sideTransec == 1,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatTransec.decP.pc2KeepSpk),1),1);
accVid_contra = decDatTransec.acc.spk(strcmp(decDatTransec.decP.focusOnL, 'video'),mainFigDatTransec.sideTransec == -1,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatTransec.decP.pc2KeepSpk),1));
accVid_normal = nanmean(decDatNormal.acc.spk(strcmp(decDatNormal.decP.focusOnL, 'video'),:,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatNormal.decP.pc2KeepSpk),1)));
plotVersus(accVid_contra,accVid_ipsi,-ones(1,numel(accVid_contra)),[1/12 1/12]*100,0,[0 100],1)
scatter(accVid_normal,accVid_normal,30,[9 133 140]/256,'s','filled')
xlabel('In contra VIS')
ylabel('In ipsi VIS');

% Audio
subplot(154)
cmp2keepneur = 4;
accSnd_ipsi = decDatTransec.acc.spk(strcmp(decDatTransec.decP.focusOnL, 'sound'),mainFigDatTransec.sideTransec == 1,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatTransec.decP.pc2KeepSpk),1),1);
accSnd_contra = decDatTransec.acc.spk(strcmp(decDatTransec.decP.focusOnL, 'sound'),mainFigDatTransec.sideTransec == -1,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatTransec.decP.pc2KeepSpk),1));
accSnd_normal = nanmean(decDatNormal.acc.spk(strcmp(decDatNormal.decP.focusOnL, 'sound'),:,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatNormal.decP.pc2KeepSpk),1)));
plotVersus(accSnd_contra,accSnd_ipsi,-ones(1,numel(accSnd_contra)),[1/12 1/12]*100,0,[0 100],1)
scatter(accSnd_normal,accSnd_normal,30,[9 133 140]/256,'s','filled')
%close up
subplot(155)
plotVersus(accSnd_contra,accSnd_ipsi,-ones(1,numel(accSnd_contra)),[1/12 1/12]*100,0,[0 100],1)
scatter(accSnd_normal,accSnd_normal,30,[9 133 140]/256,'s','filled')
xlim([10 40])
ylim([10 40])
offsetAxes

% --- Quantif ---
% Timecourses correlations
fprintf('Correlation of audPC1 ipsi vs. contra side: %d (p-value: %d).\n', c_ContravsIpsi,p_ContravsIpsi)
fprintf('Correlation of audPC1 ipsi vs. Fig1: %d (p-value: %d).\n', c_Ipsi,p_Ipsi)
fprintf('Correlation of audPC1 contra vs. Fig1: %d (p-value: %d).\n', c_Contra,p_Contra)

% Ratio
ratio = load(fullfile(transecAnatFolder,'TransectomyAnat','superprocessedData_ratioTransec.mat'));
fprintf('Ratio contra/ipsi: %s (average: %d)', num2str(1./ratio.ratio), mean(1./ratio.ratio))

% Expected explained variance
audPC1_sndCov_ipsi = audPC1_sndCov(mainFigDatTransec.sideTransec == 1);
audPC1_sndCov_contra = audPC1_sndCov(mainFigDatTransec.sideTransec == -1);
expected_audPC1_sndCov_ipsi = audPC1_sndCov_contra'.*[ones(1,2)*ratio.ratio(1) ratio.ratio(2) ones(1,3)*ratio.ratio(3)]; % hardcoded
p_exp = signrank(audPC1_sndCov_ipsi,expected_audPC1_sndCov_ipsi);
fprintf('P-value for expected vs. real: %d.\n', p_exp)

% Decoding
pVid_contra = signrank(accVid_contra,1/12*100,'tail','right');
pVid_ipsi = signrank(accVid_ipsi,1/12*100,'tail','right');
pVid = signrank(accVid_ipsi,accVid_contra);
fprintf(['Video decoding accuracy on the ipsi vs. contra side: %d +- %d (p-value: %d)' ...
    ' vs %d +- %d (p-value: %d) (comparison p-value: %d).\n'], ...
    nanmean(accVid_ipsi),nansem(accVid_ipsi),pVid_ipsi,nanmean(accVid_contra),nansem(accVid_contra),pVid_contra,pVid)

pSnd_contra = signrank(accSnd_contra,1/12*100,'tail','right');
pSnd_ipsi = signrank(accSnd_ipsi,1/12*100,'tail','right');
pSnd = signrank(accSnd_ipsi,accSnd_contra);
fprintf(['Sound decoding accuracy on the ipsi vs. contra side: %d +- %d (p-value: %d)' ...
    ' vs %d +- %d (p-value: %d) (comparison p-value: %d).\n'], ...
    nanmean(accSnd_ipsi),nansem(accSnd_ipsi),pSnd_ipsi,nanmean(accSnd_contra),nansem(accSnd_contra),pSnd_contra,pSnd)

%% Fig. 4--Prediction by movements

% Load data
mainFigDat = load(mainFigFile_fig1);
decDat = load(decodingFile_fig1);
decDatTransec = load(decodingFile_fig3);
behFigDat = load(behFigFile_fig4);
lagDat = load(lagBehNeurFile_fig4,'xc','lags_xc');
subOvDat = load(subspaceFile_fig4);
predDat = load(predFile_fig4);

% Behavior
[c_behPC1Corr,p_behPC1Corr] = plotBehFig(behFigDat.dat(find(strcmp({behFigDat.dat.name},'Body motion'))),behFigDat.behFigP);

% Decoding from motion
[accVid,accSnd,p_dec] = plotDecoding(decDat.acc.motonly,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0);

% Comparison decoding beh. and neur.
cmp2keepmot = 128;
cmp2keepneur = 4; 
accmotNormal = decDat.acc.motonly(strcmp(decDat.decP.focusOnL, 'sound'),:,find(cellfun(@(x) isequal(x,cmp2keepmot),decDat.decP.pc2KeepMot),1),1);
accspkNormal = decDat.acc.spk(strcmp(decDat.decP.focusOnL, 'sound'),:,find(cellfun(@(x) isequal(x,cmp2keepneur),decDat.decP.pc2KeepSpk),1));
[c_behvsNeur, p_behvsNeur] = plotVersus(accspkNormal,accmotNormal,side_NormalVIS,[1/12 1/12]*100,1);
accmotTransec = decDatTransec.acc.motonly(strcmp(decDatTransec.decP.focusOnL, 'sound'),:,find(cellfun(@(x) isequal(x,cmp2keepmot),decDatTransec.decP.pc2KeepMot),1),1);
accspkTransec = decDatTransec.acc.spk(strcmp(decDatTransec.decP.focusOnL, 'sound'),:,find(cellfun(@(x) isequal(x,cmp2keepneur),decDatTransec.decP.pc2KeepSpk),1));
scatter(accspkTransec(side_TransecVIS == -1),accmotTransec(side_TransecVIS == -1),20,[.7 .7 .7])
scatter(accspkTransec(side_TransecVIS == 1),accmotTransec(side_TransecVIS == 1),20,[.7 .7 .7],'x')
xlabel('From V1')
ylabel('From motion')
xlim([5 40])
ylim([5 40])

% Cross-correlation
[delays,valAtDelay,valAt0] = plotCrossCorrelogram(lagDat.xc,lagDat.lags_xc);

% Subspace overlap
[~,~,p_sub] = plotSubspaceOverlap(subOvDat,side_NormalVIS);

% Prediction
[~,p_pred] = plotPredResults(predDat,mainFigDat.sideTransec);

% --- Quantif ---
% Decoding
fprintf('Video decoding accuracy: %d +- %d (p-value: %d).\n', nanmean(accVid),nansem(accVid),p_dec(1))
fprintf('Sound decoding accuracy: %d +- %d (p-value: %d).\n', nanmean(accSnd),nansem(accSnd),p_dec(2))

% Comparison with decoding from neural data
[~,accSnd_spk,~] = plotDecoding(decDat.acc.spk,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0); close
p = signrank(accSnd,accSnd_spk);    
fprintf('Comparison of sound decoding beh. vs. neur.: p-value: %d.\n',p)

% Correlation in decoding beh vs neur
fprintf('Correlation between decoding from motion vs. from V1: %d (p-value: %s).\n',c_behvsNeur,p_behvsNeur)

% Correlation with Fig1
fprintf('Correlation with VIS audPC1 timecourse: %d (p-value: %d).\n',c_behPC1Corr,p_behPC1Corr)

% Correlation during spontaneous
fprintf('Correlation between behPC1 and neur. audPC1 during spontaneous: %d +- %d.\n', ...
    ztransform(valAtDelay,'mean'),ztransform(valAtDelay,'sem'))

% Delay during spontaneous
p_delay = signrank(delays);
fprintf('Delay between behPC1 and neur. audPC1 during spontaneous: %d +- %d (p-value: %s).\n', ...
    nanmean(delays),nansem(delays), p_delay)

% Overlap by movement subspace
fprintf('Overlap between movement-related subspace and video: %d +- %d /  and sound: %d +- %d (comparison p-value: %d).\n', ...
    nanmean(subOvDat.percExpl_vid.mov),nansem(subOvDat.percExpl_vid.mov), ...
    nanmean(subOvDat.percExpl_snd.mov),nansem(subOvDat.percExpl_snd.mov),p_sub)

% Goodness of the predictions--p-values only
fprintf('Auditory vs. behavioral model: p-value %d.\n',p_pred(1))
fprintf('Auditory vs. full model: p-value %d.\n',p_pred(2))
fprintf('Full vs. behavioral model: p-value %d.\n',p_pred(3))

%% Extended Data Fig. 1--Similarity across animals and regions

focusOn = 'video'; % or 'video'
spkVIS = load(fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephysNormal_VIS.mat',focusOn)));
spkHPF = load(fullfile(preprocFolder,'Mainfig',sprintf('allMainfig_%s_ephysNormal_HPF.mat',focusOn)));
spkPredVIS = load(fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-ridge_ephysNormal_VIS.mat',focusOn)));
spkPredHPF = load(fullfile(preprocFolder,'Model',sprintf('pred-%s_reg-ridge_ephysNormal_HPF.mat',focusOn)));
plotCompDetails(spkVIS,spkHPF,spkPredVIS,spkPredHPF,focusOn)

% Remove PC1
%%% hacky
load(fullfile(preprocFolder,'Decoding','hacky_noPC1.mat'))
decodingFile_fig1 = fullfile(preprocFolder,'Decoding',sprintf('allDecoding_ephys%s_%s.mat','Normal','VIS'));
load(decodingFile_fig1,'decP');
switch focusOn
    case 'sound'
        idx = 1;
        yl = [0 25];
    case 'video'
        idx = 2;
        yl = [75 100];
end
accNoPC1 = accNoPC1.spk(strcmp(decP.focusOnL, focusOn),:,idx);
accNormal = accNormal.spk(strcmp(decP.focusOnL, focusOn),:,idx);
figure('Position',[680   865   120   110]);
hold all
for k = 1:numel(accNoPC1)
    plot([1 2], [accNormal(k) accNoPC1(k)],'o-','color','k')
end
xlim([0 3])
xticks([1 2])
xticklabels({'All PCs','No PC1'})
xtickangle(45)
p = signrank(accNoPC1,accNormal);
sigstar({[1 2]},p)
offsetAxes
hline(1/12*100)
ylim(yl)

%% Extended Data Fig. 2--Eigenspectra

% VIS
eigDat = load(eigenspectraFile_fig1,'r','compDimP');
plotEigenspectra(eigDat.r,eigDat.compDimP,side_NormalVIS);

% HPF
eigDat = load(eigenspectraFile_fig2,'r','compDimP');
plotEigenspectra(eigDat.r,eigDat.compDimP,side_NormalHPF);

% VIS transectomy
eigDat = load(eigenspectraFile_fig3,'r','compDimP');
plotEigenspectra(eigDat.r,eigDat.compDimP,side_TransecVIS);

%% Extended Data Fig. 4--Timecourse comparison beh. and neur.

mainFigDat = load(mainFigFile_fig1);
behFigDat = load(behFigFile_fig4);

plotAllSounds(mainFigDat,behFigDat)

%% Extended Data Fig. 5--Pupil size and eye movements

behFigDat = load(behFigFile_fig4);
decDat = load(decodingFile_fig1);

% Pupil area
plotBehFig(behFigDat.dat(find(strcmp({behFigDat.dat.name},'Pupil area'))),behFigDat.behFigP);
% Decoding from pupil area
[accVid,accSnd,p_dec] = plotDecoding(decDat.acc.puparea,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0);

% Pupil motion
plotBehFig(behFigDat.dat(find(strcmp({behFigDat.dat.name},'Pupil motion'))),behFigDat.behFigP);
% Decoding from
[accVid,accSnd,p_dec] = plotDecoding(decDat.acc.pupmot,decDat.decP,mainFigDat.sideTransec,[1/12 1/12]*100,0);

%% Extended Data Fig. 6--Crosscorrelation and subspace overlap for HPF and transectomy experiments

lagDat = load(lagBehNeurFile_hpf_EDfig6,'xc','lags_xc');
subOvDat = load(subspaceFile_hpf_EDfig6);
plotCrossCorrelogram(lagDat.xc,lagDat.lags_xc)
plotSubspaceOverlap(subOvDat,side_NormalHPF);

lagDat = load(lagBehNeurFile_transec_EDfig6,'xc','lags_xc');
subOvDat = load(subspaceFile_transec_EDfig6);
plotCrossCorrelogram(lagDat.xc,lagDat.lags_xc)
plotSubspaceOverlap(subOvDat,side_TransecVIS);

%% Extended Data Fig. 7--Weights

predDat = load(predFile_fig4);
weights = load(facemapWeightsFile_EDfig7);
plotPredResults(predDat) % second plot is the one
plotPredFace(predDat,weights.facemapWeights)%,P.fileref.rootDir)

%% Extended Data Fig. 8--Predictions by models

% VIS
predDat = load(predFile_video_NormalVIS_EDfig8);
plotPredResults(predDat,side_NormalVIS)

% HPF
predDat = load(predFile_sound_NormalHPF_EDfig8);
plotPredResults(predDat,side_NormalHPF)

predDat = load(predFile_video_NormalHPF_EDfig8);
plotPredResults(predDat,side_NormalHPF)

% Transectomy
predDat = load(predFile_sound_TransecVIS_EDfig8);
plotPredResults(predDat,side_TransecVIS)

predDat = load(predFile_video_TransecVIS_EDfig8);
plotPredResults(predDat,side_TransecVIS)

%% Extended Data Fig. 9--Single trial fluctuations

predDat = load(predFile_fig4);
mainFigDat = load(mainFigFile_fig1,'mainFigP');
[~,z_ST,p_ST] = plotPredSingleTrial(predDat,mainFigDat.mainFigP);

% --- Quantif ---
fprintf('Correlation for single trials : %d (p-value: %s).\n',z_ST,p_ST)

%% Extended Data Fig. 10--Lag analysis

lagDat = load(lagBehNeurFile_fig4,'w','lags_w');
delays = plotLags(lagDat.w,lagDat.lags_w);

lagDat = load(lagBehNeurFile_hpf_EDfig6,'w','lags_w');
plotLags(lagDat.w,lagDat.lags_w)

lagDat = load(lagBehNeurFile_transec_EDfig6,'w','lags_w');
plotLags(lagDat.w,lagDat.lags_w)

% --- Quantif ---
p_delay = signrank(delays);
fprintf('Delay between behPC1 and neur. audPC1 during spontaneous: %d +- %d (p-value: %s).\n', ...
    nanmean(delays),nansem(delays), p_delay)

%% Suppl. Fig. XXX--Lag between ipsi and contra

lagICDat = load(lagIpsiContraFile);
xc = plotCrossCorrelogram(lagICDat.pc1_ipsi,lagICDat.pc1_contra,lagICDat.lags_xc); title('Ipsi x contra')
