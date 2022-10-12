clear all
close all

%% Load CCF atlas
allen_atlas_path = '\\zserver\Lab\Atlas\allenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);

%% get all voxels of visual cortex

% get the voxels. There should be a much easier way
allen_atlas_path = '\\zserver\Lab\Atlas\allenCCF';
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% get visual areas
visArea = {'VISp','VISpm','VISl'}; % VISp for primary visual cortex, VIS for all visual cortex
visID = [];
for a = 1:numel(visArea)
    struct_id = find(strcmpi(st.acronym,visArea{a}));
    % get all subareas
    visID = [visID find(contains(st.structure_id_path,st.structure_id_path{struct_id}))];
end    

% get auditory areas -- just for display
audArea = {'AUD'}; % just primary visual cortex
audID = [];
for a = 1:numel(audArea)
    struct_id = find(strcmpi(st.acronym,audArea{a}));
    % get all subareas
    audID = [audID find(contains(st.structure_id_path,st.structure_id_path{struct_id}))];
end 

midline = 5739;
subsamp = 500; % quick and dirty subsampling to avoid taking too much time for computation

% get their indices and coordinates - VIS
idx = find(ismember(av,visID));
[i,j,k] = ind2sub(size(av),unique(idx));
i = 10*i(1:subsamp:end); j = 10*j(1:subsamp:end); k = 10*k(1:subsamp:end);

clear vis
vis.vR(:,1) = i(k<midline); vis.vR(:,2) = j(k<midline); vis.vR(:,3) = k(k<midline);
vis.vL(:,1) = i(k>midline); vis.vL(:,2) = j(k>midline); vis.vL(:,3) = k(k>midline);
clear i j k

% get their indices and coordinates - AUD
idx = find(ismember(av,audID));
[i,j,k] = ind2sub(size(av),unique(idx));
i = 10*i(1:subsamp:end); j = 10*j(1:subsamp:end); k = 10*k(1:subsamp:end);

clear aud
aud.vR(:,1) = i(k<midline); aud.vR(:,2) = j(k<midline); aud.vR(:,3) = k(k<midline);
aud.vL(:,1) = i(k>midline); aud.vL(:,2) = j(k>midline); aud.vL(:,3) = k(k>midline);
clear i j k

%% Animal references
animals = {'CR_transectomy3','Transectomy7','Transectomy8'};

for ani = 1:numel(animals)
    animal = animals{ani};
    cut{ani} = readObj(fullfile('D:\BrainSaw',animal,'downsampled_stacks\025_micron\brainreg_output\manual_segmentation\standard_space\regions\cut.obj'));
end

%% Load axon tractographies
streamlines_path = 'C:\Users\Hamish\Desktop\InjectionsA1_AllenBrainAtlas';

% get info about experiments
T = readtable(fullfile(streamlines_path, 'projection_search_results.csv'));

% loop over all of them
d = dir(fullfile(streamlines_path,'*.json'));
% d = dir(fullfile(streamlines_path,'*606100558.json'));

% processing plotting and saving options
plt = 0;
pltAll = 0;
addInjSite = 1; % add injection site to trajectory
saveJSON = 0;
distlimCut = 50;
distlimVIS = 50;

% send warning as this isn't the most conservative way, yet might be needed
%%% discuss this
if addInjSite
    warning('Warning. Adding injection site at the beginning of each branche to compute intersect with the cut. Not conservative.')
end
    
%%% note that here might be some weird fibers coming out of nowhere that
%%% will still be counted -- but that's a first attempt. It can only
%%% decrease the amount we actually cut?

%%% also, take only the endpoint? -- otherwise, some axons are just passing by
%%% V1 but then really seem to project only to some other place?
%%% or remove the ones of which starting point is very far from injection
%%% site

%%% finally, note that some fibers will "start" quite away from injection
%%% point -- not sure how to deal with that, esp. if injection point is
%%% right at the cut, and fibers start on the other side!

A = struct();
tic;
for f = 1:numel(d)
    %% load the specific experiment
    filePath = fullfile(d(f).folder,d(f).name); % filename in JSON extension
    fid = fopen(filePath); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    data = jsondecode(str);
    
    %% get exp info
    expID = find(T.id==str2num(d(f).name(strfind(d(f).name,'_')+1:strfind(d(f).name,'.json')-1)));
    expNum(f) = T.id(expID);
    struct{f} = T.structure_abbrev{expID};
    injSite{f} = str2num(T.injection_coordinates{expID}); % not sure how good these injection site coordinates are
    injSiteStruct(f) = av(floor(injSite{f}(1)/10),floor(injSite{f}(2)/10),floor(injSite{f}(3)/10)); % maybe discrepancy between the two?
    
    %% mirror it if on the contra side
    % check side of the beginning, flip from midline if wrong side
    if injSite{f}(3) > midline
        for l = 1:numel(data.lines)
            for ii = 1:numel(data.lines{l})
                data.lines{l}(ii).z = 2*midline - data.lines{l}(ii).z;
            end
        end
        injSite{f}(3) = 2*midline - injSite{f}(3);
        for l = 1:numel(data.injection_sites)
            data.injection_sites(l).z = 2*midline - data.injection_sites(l).z;
        end
    end
    
    % get projections to VIS
    idxVIS = [];
    RorL = [];
    for l = 1:numel(data.lines)
        branch = [[data.lines{l}.x]; [data.lines{l}.y]; [data.lines{l}.z]];
        d2visR = nan(1,size(branch,2));
        d2visL = nan(1,size(branch,2));
        d2injsite = sqrt(sum((injSite{f} - branch(:,end)').^2,2));
        for p = 1:size(branch,2)
            d2visR(p) = min(sqrt(sum((vis.vR - branch(:,p)').^2,2)));
            d2visL(p) = min(sqrt(sum((vis.vL - branch(:,p)').^2,2)));
        end
        if (d2visR(1) < distlimVIS || d2visL(1) < distlimVIS) && d2injsite < 1500 % condition where only the end of the tract is going in V1
            %         if (any(d2visR < distlimVIS) || any(d2visL < distlimVIS)) && d2injsite < 1500 % condition where any piece of the tract is going through
            idxVIS = [idxVIS l]; % select the branches that pass by V1
            RorL = [RorL -(d2visL(1) < distlimVIS)+(d2visR(1) < distlimVIS)]; % -1: contra / 1: ipsi
            % RorL = [RorL -any(d2visL < distlimVIS)+any(d2visR < distlimVIS)]; % condition where any piece of the tract is going through
        end
        
        % quantify that better somehow? the length of the portion that is in V1?
    end
    
    %% quantify how much the cut removed
    for ani = 1:numel(animals)
        % get which of these fibers was actually cut
        idxCutIpsi = []; % if the injection is ipsi to the cut
        idxCutContra = []; % if the injection is contra to the cut
        for ll = 1:numel(idxVIS)
            l = idxVIS(ll);
            branch = [[data.lines{l}.x]; [data.lines{l}.y]; [data.lines{l}.z]];
            if addInjSite
                % not conservative
                branch = [branch ...
                    [linspace(branch(1,end),injSite{f}(1),10); ...
                    linspace(branch(2,end),injSite{f}(2),10); ...
                    linspace(branch(3,end),injSite{f}(3),10)]];
            end
            d2ipsicut = nan(1,size(branch,2));
            d2contracut = nan(1,size(branch,2));
            for p = 1:size(branch,2)
                d2ipsicut(p) = min(sqrt(sum((cut{ani}.v - branch(:,p)').^2,2)));
                d2contracut(p) = min(sqrt(sum((cut{ani}.v - ([0 0 2*midline] + [1 1 -1].* branch(:,p)')).^2,2))); % if injection on the contra side
            end
            if any(d2ipsicut < distlimCut)
                idxCutIpsi = [idxCutIpsi l];
            end
            if any(d2contracut < distlimCut)
                idxCutContra = [idxCutContra l];
            end
        end
        
        %% save all
        A(ani).fibers2IpsiVIS_tot(f) = sum(RorL > -1);
        A(ani).fibers2ContraVIS_tot(f) = sum(RorL < 1);
        A(ani).fibers2IpsiVIS_cut(f) = sum(RorL(ismember(idxVIS,idxCutIpsi)) > -1); % fibers going in ipsi V1
        A(ani).fibers2ContraVIS_cut(f) = sum(RorL(ismember(idxVIS,idxCutIpsi)) < 1); % fibers going in contra V1
        A(ani).fibers2ContraVIS_cut_fromContra(f) = sum(RorL(ismember(idxVIS,idxCutContra)) > -1); % this should be zero? / flipped because theoretically on the other side
        A(ani).fibers2IpsiVIS_cut_fromContra(f) = sum(RorL(ismember(idxVIS,idxCutContra)) < 1); % flipped because theoretically on the other side
        
        %% plot all
        if plt
            %%% plot only the case where injection is ipsi to cut
            
            figure('Name',['Exp. #' num2str(T.id(expID))]);
            subplot(1,2,1); hold all
            histogram(RorL,[-2,-0.5,0.5,2])
            histogram(RorL(ismember(idxVIS,idxCutIpsi)),[-2,-0.5,0.5,2])
            
            ax = subplot(1,2,2);
            % taken from plotBrainGrid
            brainGridData = readNPY(fullfile('brainGridData.npy')); % don't know how it works..?
            brainGridData = double(brainGridData)*10; % 10um
            bp = double(brainGridData);
            bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there wer
            h = plot3(bp(:,1), bp(:,3), bp(:,2), 'Color', [0 0 0 0.3]);
            set(ax, 'ZDir', 'reverse')
            axis(ax, 'equal');
            axis(ax, 'vis3d');
            axis(ax, 'off');
            hold all
            scatter3(injSite{f}(1),injSite{f}(2),injSite{f}(3),200,'k','filled'); % injection site
            text(injSite{f}(1),injSite{f}(2),injSite{f}(3),'INJ. SITE')
            scatter3(cut{ani}.v(:,1),cut{ani}.v(:,2),cut{ani}.v(:,3),'g'); % cut
            scatter3(vis.vL(:,1),vis.vL(:,2),vis.vL(:,3),50,[0 0 0.9]); % VIS
            scatter3(vis.vR(:,1),vis.vR(:,2),vis.vR(:,3),50,[0 0 0.9]); % VIS
            scatter3(aud.vL(:,1),aud.vL(:,2),aud.vL(:,3),50,[0 0 0.7]); % AUD
            scatter3(aud.vR(:,1),aud.vR(:,2),aud.vR(:,3),50,[0 0 0.7]); % AUD
            if pltAll
                branch2plt = 1:numel(data.lines);
            else
                branch2plt = idxVIS;
                if isempty(idxVIS)
                    branch2plt = 1:numel(data.lines);
                    fprintf('No fibers to V1 found on exp #%s. Plot all to check.\n',num2str(T.id(expID)))
                end
            end
            for l = branch2plt
                branch = [[data.lines{l}.x]; [data.lines{l}.y]; [data.lines{l}.z]];
                if addInjSite
                    % not conservative
                    branch = [branch ...
                        [linspace(branch(1,end),injSite{f}(1),10); ...
                        linspace(branch(2,end),injSite{f}(2),10); ...
                        linspace(branch(3,end),injSite{f}(3),10)]];
                end
                % recompute the cut if we want it all.
                d2cut = nan(1,size(branch,2));
                for p = 1:size(branch,2)
                    d2cut(p) = min(sqrt(sum((cut{ani}.v - branch(:,p)').^2,2)));
                end
                if any(d2cut < distlimCut)
                    plot3(branch(1,:),branch(2,:),branch(3,:),'r')
                else
                    plot3(branch(1,:),branch(2,:),branch(3,:),'k')
                end
            end
        end
        
        %% save new .json files with i) all V1 fibers and ii) only uncut V1 fibers
        if saveJSON
            data_allV1cut = data;
            data_allV1cut.lines = data_allV1cut.lines(idxCutIpsi);
            data_allV1uncut = data;
            data_allV1uncut.lines = data_allV1uncut.lines(idxVIS(~ismember(idxVIS,idxCutIpsi)));
            
            % all the V1 fibers
            JSONFILE_name = fullfile(d(f).folder,'processed',[animal '_CutV1fibers_' d(f).name]);
            fid = fopen(JSONFILE_name,'w');
            encodedJSON = jsonencode(data_allV1cut);
            fprintf(fid, encodedJSON);
            fclose(fid);
            
            % only the uncut V1 fibers
            JSONFILE_name = fullfile(d(f).folder,'processed',[animal '_UncutV1fibers_' d(f).name]);
            fid = fopen(JSONFILE_name,'w');
            encodedJSON = jsonencode(data_allV1uncut);
            fprintf(fid, encodedJSON);
            fclose(fid);
        end
    end
end

for ani = 1:numel(animals)
    % clean the ones where there was nothing going to V1
    A(ani).fibers2IpsiVIS_tot(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    A(ani).fibers2IpsiVIS_cut(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    A(ani).fibers2ContraVIS_tot(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    A(ani).fibers2ContraVIS_cut(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    A(ani).fibers2ContraVIS_cut_fromContra(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    A(ani).fibers2IpsiVIS_cut_fromContra(A(ani).fibers2IpsiVIS_tot == 0) = nan;
    
    injSite2Exlude = [169]; % no control over injection site on this one
    A(ani).fibers2IpsiVIS_tot(ismember(injSiteStruct,injSite2Exlude)) = nan;
    A(ani).fibers2IpsiVIS_cut(ismember(injSiteStruct,injSite2Exlude)) = nan;
    A(ani).fibers2ContraVIS_tot(ismember(injSiteStruct,injSite2Exlude)) = nan;
    A(ani).fibers2ContraVIS_cut(ismember(injSiteStruct,injSite2Exlude)) = nan;
    A(ani).fibers2ContraVIS_cut_fromContra(ismember(injSiteStruct,injSite2Exlude)) = nan;
    A(ani).fibers2IpsiVIS_cut_fromContra(ismember(injSiteStruct,injSite2Exlude)) = nan;
    
    ratio(ani) = nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut + A(ani).fibers2ContraVIS_tot-A(ani).fibers2IpsiVIS_cut_fromContra)/ ...
    nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2ContraVIS_cut_fromContra + A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut);
end
toc

%% plot final results -- for one animal

structList = unique(struct);
cmap = winter(numel(structList));
colorDot = zeros(numel(struct),3);
for s = 1:numel(struct)
    colorDot(s,:) = cmap(strcmp(structList,struct{s}),:);
end

figure('Position', [771   818   469   160]); 
for ani = 1:numel(animals)
    subplot(numel(animals),3,(ani-1)*3+1); hold all
    for s = 1:numel(structList) % convoluted way just to get the legend...
        sc(s) = scatter(A(ani).fibers2IpsiVIS_tot(strcmp(struct,structList{s})), ...
            A(ani).fibers2IpsiVIS_tot(strcmp(struct,structList{s}))-A(ani).fibers2IpsiVIS_cut(strcmp(struct,structList{s})), ...
            20,colorDot((strcmp(struct,structList{s})),:),'filled');
    end
    scatter(nanmean(A(ani).fibers2IpsiVIS_tot),nanmean(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut),40,'k','filled')
    plot([0 max(A(ani).fibers2IpsiVIS_tot)],[0,max(A(ani).fibers2IpsiVIS_tot)],'k--')
    axis equal tight
    xlabel('Initial # of fibers')
    ylabel('Remaining # of fibers')
    legend(sc,structList)
    title('Ipsi')
    subplot(numel(animals),3,(ani-1)*3+2); hold all
    scatter(A(ani).fibers2ContraVIS_tot,A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut,20,colorDot,'filled')
    scatter(nanmean(A(ani).fibers2ContraVIS_tot),nanmean(A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut),40,'k','filled')
    plot([0 max(A(ani).fibers2IpsiVIS_tot)],[0,max(A(ani).fibers2IpsiVIS_tot)],'k--')
    axis equal tight
    xlabel('Initial # of fibers')
    ylabel('Remaining # of fibers')
    title('Contra')
    
    subplot(numel(animals),3,(ani-1)*3+3); hold all
    % ipsi
    b(1) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44]);
    b(2) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0]);
    bar(2,nansum(A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44])
    bar(2,nansum(A(ani).fibers2ContraVIS_tot-A(ani).fibers2IpsiVIS_cut_fromContra),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0])
    
    % contra
    bar(4,nansum(A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44]);
    bar(4,nansum(A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0]);
    bar(5,nansum(A(ani).fibers2IpsiVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44])
    bar(5,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2ContraVIS_cut_fromContra),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0])
    
    ylabel('# of fibers')
    xticks([1 2 4 5])
    xticklabels({'In ipsiV1, from ipsi','In ipsiV1, from contra','In contraV1, from ipsi','In contraV1, from contra'})
    xtickangle(45)
    text(2, 800, 'total', 'color',[0.98 0.50 0.44])
    text(2, 600, 'remaining', 'color',[0 0 0])
end

figure('Position', [771   818   469   160]); 
for ani = 1:numel(animals)
    subplot(numel(animals),3,(ani-1)*3+1); hold all
    for s = 1:numel(structList) % convoluted way just to get the legend...
        sc(s) = scatter(nanmean(A(ani).fibers2IpsiVIS_tot(strcmp(struct,structList{s}))), ...
            nanmean(A(ani).fibers2IpsiVIS_tot(strcmp(struct,structList{s}))-A(ani).fibers2IpsiVIS_cut(strcmp(struct,structList{s}))), ...
            20,cmap(s,:),'filled');
    end
    scatter(nanmean(A(ani).fibers2IpsiVIS_tot),nanmean(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut),40,'k','filled')
    plot([0 max(A(ani).fibers2IpsiVIS_tot)],[0,max(A(ani).fibers2IpsiVIS_tot)],'k--')
    axis equal tight
    xlabel('Initial # of fibers')
    ylabel('Remaining # of fibers')
    legend(sc,structList)
    title('Ipsi')
    subplot(numel(animals),3,(ani-1)*3+2); hold all
    for s = 1:numel(structList) % convoluted way just to get the legend...
        sc(s) = scatter(nanmean(A(ani).fibers2ContraVIS_tot(strcmp(struct,structList{s}))), ...
            nanmean(A(ani).fibers2ContraVIS_tot(strcmp(struct,structList{s}))-A(ani).fibers2ContraVIS_cut(strcmp(struct,structList{s}))), ...
            20,cmap(s,:),'filled');
    end
    scatter(nanmean(A(ani).fibers2ContraVIS_tot),nanmean(A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut),40,'k','filled')
    plot([0 max(A(ani).fibers2IpsiVIS_tot)],[0,max(A(ani).fibers2IpsiVIS_tot)],'k--')
    axis equal tight
    xlabel('Initial # of fibers')
    ylabel('Remaining # of fibers')
    title('Contra')
    
    subplot(numel(animals),3,(ani-1)*3+3); hold all
    % ipsi
    b(1) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44]);
    b(2) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0]);
    bar(2,nansum(A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44])
    bar(2,nansum(A(ani).fibers2ContraVIS_tot-A(ani).fibers2IpsiVIS_cut_fromContra),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0])
    
    % contra
    bar(4,nansum(A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44]);
    bar(4,nansum(A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0]);
    bar(5,nansum(A(ani).fibers2IpsiVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 1.0,'EdgeColor',[0.98 0.50 0.44])
    bar(5,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2ContraVIS_cut_fromContra),'FaceColor',[0.6 0 0], 'Facealpha', 1.0,'EdgeColor',[0.6 0 0])
    
    ylabel('# of fibers')
    xticks([1 2 4 5])
    xticklabels({'In ipsiV1, from ipsi','In ipsiV1, from contra','In contraV1, from ipsi','In contraV1, from contra'})
    xtickangle(45)
    text(2, 800, 'total', 'color',[0.98 0.50 0.44])
    text(2, 600, 'remaining', 'color',[0 0 0])
end

% or just sum of all
figure('Position', [989   717   251   261]); hold all
b(1) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot+A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 0.4,'EdgeColor',[0.98 0.50 0.44]);
b(2) = bar(1,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut + A(ani).fibers2ContraVIS_tot-A(ani).fibers2IpsiVIS_cut_fromContra),'FaceColor',[0.6 0 0], 'Facealpha', 0.4,'EdgeColor',[.6 0 0]);
bar(2,nansum(A(ani).fibers2IpsiVIS_tot+A(ani).fibers2ContraVIS_tot),'FaceColor',[0.98 0.50 0.44], 'Facealpha', 0.4,'EdgeColor',[0.98 0.50 0.44])
bar(2,nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2ContraVIS_cut_fromContra + A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut),'FaceColor',[0.6 0 0], 'Facealpha', 0.4,'EdgeColor',[.6 0 0])
ylabel('# of fibers')
xticks([1 2])
xticklabels({'In ipsiV1','In contraV1'})
xtickangle(45)
text(3, 800, 'total', 'color',[0.98 0.50 0.44])
text(3, 750, 'remaining', 'color',[0 0 0])
title(sprintf('Expected ratio: %s',ratio(ani)))

%% plot for all animals

figure; 
subplot(1,3,1:2); hold all
for ani = 1:numel(animals)
    m(ani) = nansum(A(ani).fibers2IpsiVIS_tot +  A(ani).fibers2ContraVIS_tot);
    scatter(nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2ContraVIS_cut_fromContra + A(ani).fibers2ContraVIS_tot-A(ani).fibers2ContraVIS_cut), ...
        nansum(A(ani).fibers2IpsiVIS_tot-A(ani).fibers2IpsiVIS_cut + A(ani).fibers2ContraVIS_tot-A(ani).fibers2IpsiVIS_cut_fromContra),40,'k','filled')
end
axis equal tight
xlim([0 max(m)])
ylim([0 max(m)])
plot([0 max(m)],[0 max(m)],'k--')
xlabel('Aud. input to contra')
ylabel('Aud. input to ipsi');

subplot(1,3,3); hold all
scatter(1:3,ratio,40,'k','filled')
ylim([0 1])
ylabel('Ratio');
xlabel('Mouse')
