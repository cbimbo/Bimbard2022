function plotLFPandCSD(LFP,CSD,depthChan,timebins,cluXPos,cluDepths,anatDir,granBound)

    % Get Allen Brain Atlas references
    allen_atlas_path = '\\zserver\Lab\Atlas\allenCCF';
    st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
    cmap_filename = [allen_atlas_path filesep 'allen_ccf_colormap_2017.mat'];
    load(cmap_filename);

    figure('Position',[579   583   467   395]);
    clear ax
    numXidx = numel(LFP);

    % Plot LFP
    for xposIdx = 1:numXidx
        ax(xposIdx) = subplot(1,3*numXidx,xposIdx);
        imagesc(timebins,depthChan{xposIdx},interp2(LFP{xposIdx},1))
        set(gca,'YDir','normal')
        colormap(ax(xposIdx),'RedBlue'); caxis([-max(abs(LFP{xposIdx}(:))) max(abs(LFP{xposIdx}(:)))])
        if xposIdx == 1
            ylabel('Distance to tip of probe')
            xlabel('Time (s)')
        end

        if exist('granBound','var')
            hline(granBound{xposIdx},'k--')
        end
        set(gca,'box','off')
        xlim([-0.01 0.08])
    end

    % Plot CSD
    for xposIdx = 1:numXidx
        ax(numXidx+xposIdx) = subplot(1,3*numXidx,numXidx+xposIdx);
        imagesc(timebins,depthChan{xposIdx}(3:end-2),interp2(CSD{xposIdx},1))
        set(gca,'YDir','normal')
        colormap(ax(numXidx+xposIdx),'RedBlue'); caxis([-max(abs(CSD{xposIdx}(:))) max(abs(CSD{xposIdx}(:)))])

        if exist('granBound','var')
            hline(granBound{xposIdx},'k--')
        end
        set(gca,'box','off')
        xlim([-0.01 0.08])
    end
    
    % Plot neuron location
    ax(3*numXidx+1) = subplot(1,3*numXidx,2*numXidx+1:3*numXidx);
    scatter(cluXPos,cluDepths)
    if ~isempty(anatDir)
        if contains(anatDir,'slices')
            % get specific files
            load(fullfile(anatDir,'probe_ccf.mat'))
            trajDepth = 3840-probe_ccf.probe_depths;
            trajArea = probe_ccf.trajectory_areas;

            curr_probe = 1;
            trajectory_area_boundaries = ...
                [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0);length(probe_ccf(curr_probe).trajectory_areas)];
            trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2;
            trajectory_area_labels = st.acronym(probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers)));
        end
        hold all
        imagesc(10:60,trajDepth,repmat(trajArea,[1,2]))
        colormap(ax(3*numXidx+1),cmap);
        caxis([1,size(cmap,1)])
        set(gca,'YTick',trajDepth(round(trajectory_area_centers(end:-1:1))),'YTickLabels',trajectory_area_labels(end:-1:1));
        set(gca,'YDir','Normal')
        hline(trajDepth(trajectory_area_boundaries))

        if exist('granBound','var')
            hline(granBound{xposIdx},'k--')
        end
    end
    % plot clusters
    hold all
    scatter(cluXPos,cluDepths,40,'k','filled') % all in black
    ylim([min(cluDepths)-10,max(cluDepths)+10])
    xlim([min(cluXPos)-10,max(cluXPos)])
    linkaxes(ax,'y')
    linkaxes(ax(1:end-1),'x')