function [db,P,Nc] = subselectCells(db,P,focusOnArea,pltFig)
    %%% This function needs the allenCCF repo to work (https://github.com/cortex-lab/allenCCF),
    %%% or at least the 'structure_tree_safe_2017.csv'

    if ~exist('pltFig','var')
        pltFig = 0;
    end

    %% Subselect the cells according to various criteria
    % UNIT TYPE
    P.proc.lab2keep = [2 3.2];
    
    % SPIKE COUNT
    P.proc.spkNumThre = 300;
    
    % AREA
    % get areas name and color properties
    allen_atlas_path = P.path2CCF;
    st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
    cmap_filename = [allen_atlas_path filesep 'allen_ccf_colormap_2017.mat'];
    load(cmap_filename);

    % ISI VIOLATIONS
    P.proc.ISIViolThre = inf;
    P.proc.tauC = 0.0001;
    P.proc.tauR = 0.0015;

    switch focusOnArea
        case {'ephys','ephysTrans','VIS'} % all visual cortex
            P.proc.areas2keep = {'VIS'};
        case 'VISp'
            P.proc.areas2keep ={'VISp'}; % just primary visual cortex
        case {'ephysHPC','HPF'}
            P.proc.areas2keep ={'HPF'}; % hippocampal areas
        case {'all','root'}
            P.proc.areas2keep = {'root'};
            %         P.proc.areas2keep = {'grey'};
        case 'MG'
            P.proc.areas2keep = {'MG'};
        otherwise
            error('Check area name')
    end
    
    % Get all substructures of each
    areas2keep_idx = [];
    for a = 1:numel(P.proc.areas2keep)
        struct_id = find(strcmpi(st.acronym,P.proc.areas2keep{a}));
        % get all subareas
        areas2keep_idx = [areas2keep_idx find(contains(st.structure_id_path,st.structure_id_path{struct_id}))];
    end
    
    % Get subselection indices
    mouse2Keep = []; % in case no neuron left
    idx2keep = cell(1,P.nMice);
    for k = 1:P.nMice
        % subselect unit type
        idxUnitType = ismember(db(k).C.CluLab, P.proc.lab2keep); 
        
        % subselect spike number
        idxSpknum = db(k).C.CluSpknum > P.proc.spkNumThre;
        
        % subselect depth
        idxDepth = db(k).C.Depth > 0; % just for tests
        
        % subselect area
        idxArea = ismember(db(k).C.area,areas2keep_idx);
        
        % subselect ISI violations
        for c = 1:numel(db(k).C.CluID)
            clu = db(k).C.CluID(c);
            spikeTimes = db(k).sp.st(db(k).sp.clu == clu);
            [db(k).C.ISIViol(c), ~] = computeISIViolations(spikeTimes,P.proc.tauR,P.proc.tauC);
        end
        idxISIViol = (db(k).C.ISIViol < P.proc.ISIViolThre)';

        % make sure there's no duplicate
        tmp = reshape(permute(db(k).spikeData,[1 3 2]),[P.nBins*size(db(k).spikeData,3),size(db(k).spikeData,2)]);
        corrtmp = corr(tmp,tmp);
        idxDuplicate = ~(sum(triu(corrtmp)>0.99)>1)';
        
        % subselect combined
        idx2keep{k} = idxArea & idxDepth & idxUnitType & idxSpknum & idxISIViol & idxDuplicate;
        
        if ~all(~idx2keep{k})
            mouse2Keep = [mouse2Keep k];
        end
    end
    
    %% Plot figure if asked
    if pltFig
        % figure
        %%% will have to adapt it with several shanks....
        for k = 1:P.nMice
            % make nice figure for each animal
            figure;
            % plot areas on the background

            if ~isempty(P.fileref.anatDir{k})
                if contains(P.fileref.anatDir{k},'slices')
                    % get specific files
                    load(fullfile(P.fileref.anatDir{k},'probe_ccf.mat'))
                    trajDepth = 3840-probe_ccf.probe_depths;
                    trajArea = probe_ccf.trajectory_areas;

                    curr_probe = 1;
                    trajectory_area_boundaries = ...
                        [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0);length(probe_ccf(curr_probe).trajectory_areas)];
                    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2;
                    trajectory_area_labels = st.acronym(probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers)));
                end
                imagesc(10:60,trajDepth,repmat(trajArea,[1,2]))
                colormap(cmap);
                caxis([1,size(cmap,1)])
                set(gca,'YTick',trajDepth(round(trajectory_area_centers(end:-1:1))),'YTickLabels',trajectory_area_labels(end:-1:1));
                set(gca,'YDir','Normal')
            end
            hold all
            % plot clusters
            scatter(db(k).C.XPos,db(k).C.Depth,40,'k','filled') % all in black
            scatter(db(k).C.XPos(idx2keep{k}),db(k).C.Depth(idx2keep{k}),'r','filled') % the ones that are kept will be plotted in red
            ylim([min(db(k).C.Depth)-10,max(db(k).C.Depth)+10])
            xlim([min(db(k).C.XPos),max(db(k).C.XPos)])
        end
    end

    %% Save subselection
    for k = 1:P.nMice
        db(k).spikeData = db(k).spikeData(:,idx2keep{k},:);
        db(k).C.XPos = db(k).C.XPos(idx2keep{k});
        db(k).C.Depth = db(k).C.Depth(idx2keep{k});
        db(k).C.CluID = db(k).C.CluID(idx2keep{k});
        db(k).C.CluLab = db(k).C.CluLab(idx2keep{k});
        db(k).C.CluSpknum = db(k).C.CluSpknum(idx2keep{k});
        db(k).C.area = db(k).C.area(idx2keep{k});
        db(k).C.ISIViol = db(k).C.ISIViol(idx2keep{k});
    end
    
    db = db(mouse2Keep);
    P.nMice = numel(mouse2Keep);
    P.mouseColor = P.mouseColor(mouse2Keep,:);
    P.mouseRef = P.mouseRef(mouse2Keep);
    P.miceweye = P.miceweye(mouse2Keep);
    P.fileref.blockDir = P.fileref.blockDir(mouse2Keep);
    P.fileref.rootDir = P.fileref.rootDir(mouse2Keep);
    P.fileref.ksDir = P.fileref.ksDir(mouse2Keep);
    P.fileref.anatDir = P.fileref.anatDir(mouse2Keep);

    Nc = nan(1,P.nMice);
    for k = 1:P.nMice
        Nc(k) = size(db(k).spikeData,2);
    end

    P.focusOnArea = focusOnArea;
end