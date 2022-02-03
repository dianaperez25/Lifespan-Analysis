%% CreateVariantFiles.m
%This script makes the variants from spatial correlation maps, excluding
%regions with low signal. 
%Written by Brian Kraus. Edited by Diana

clear all

%% Paths
%change paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';
data_dir = '/Users/dianaperez/Desktop/'; %directory where the data is located
rest_file = [data_dir 'sub-LS05_spCorrMap.dtseries.nii']; % name of the spatial correlation map file
SNRpath = '/Users/dianaperez/Box/HCP_variants/bottomBrainMask.dtseries.nii';
outfilepath = '/Users/dianaperez/Desktop/';
subject = 'LS05';

%%
threshold = [5];  %% Thresholds used to calculate variants (lowest % or correlation values)
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR, 1 = exclude, 0 = don't exclude
ExcludeBySize = 1;
SNRmap = ft_read_cifti_mod(SNRpath);

%%
%     % reads file paths, sub numbers split-halves from txt files    
%     [task_files, subjects1, tasksplithalf] = textread([dirpath tasktxtname],'%s%s%s');
%     [rest_files, subjects2, restsplithalf] = textread([dirpath resttxtname],'%s%s%s');

% main for-loop, creates variant maps
for x = 1:numel(threshold)

        %subject = subjects2{x};
        
        %reads cifti files from the txt files
        cifti_rest = ft_read_cifti_mod(rest_file);
        %%
        % if you want to exclude low signal regions, this will exclude by SNR
        if SNRexclusion == 1
            SNRmap = ft_read_cifti_mod(SNRpath);
            
%             SNRmap.data = SNRmap.data(1:59412,:);
%             SNRexclude = find(SNRmap.data < 750);
            
            cifti_rest.data(SNRmap.data==1) = NaN;
            %cifti_task.data(SNRexclude,1) = NaN;
        end
        
        %% Makes variant maps
        %this part makes the variant mask by finding the vertices that are
        %below the threshold, making those vertices = 1 and all other
        %vertices = 0
        %cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
        cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold(x)));

        cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
        cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;

%         cifti_task_thresh_dat = zeros(size(cifti_task.data));
%         cifti_task_thresh_dat(cifti_task_threshold,1) = 1;
        
        %%
        % this does exactly what the six lines above do?
        cifti_rest_final_dat = zeros(size(cifti_rest.data));
        %cifti_task_final_dat = zeros(size(cifti_task.data));

        for w = 1:length(cifti_rest.data)

            if cifti_rest_thresh_dat(w) == 1

                cifti_rest_final_dat(w) = 1;

            end

%             if cifti_task_thresh_dat(w) == 1
% 
%                 cifti_task_final_dat(w) = 1;
% 
%             end
        end
        %%

        %This creates the output file names for rest and task 
        %outfilerest = strrep(rest_file{x}, 'vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold(x))]);
        %outfiletask = strrep(task_files{x}, 'cortex_vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold)]);
        outfilerest = [outfilepath '/' subject '_uniqueIDs_variants_sizeExcluded_thresh-' num2str(threshold(x)) '_smooth_2.55.dtseries.nii'];
        %strrest = ['SNRExclude_' restsplithalf{x}];
        %strtask = ['SNRExclude_' tasksplithalf{x}];
        %outfilewbtask = [outfilepath '/' subject '_matcheddata_Variant_Size_' strtask '_' num2str(threshold) '.dtseries.nii'];
        outfilewbrest = [outfilepath '/' subject '_uniqueIDs_matcheddata_REST_Variant_Size_thresh-' num2str(threshold(x)) '_smooth_2.55.dtseries.nii'];

        %this creates and writes the file in cifti format
        cifti_rest.data = cifti_rest_final_dat;
        %cifti_task.data = cifti_task_final_dat;
        
        ft_write_cifti_mod(outfilerest, cifti_rest)
        %ft_write_cifti_mod(outfiletask, cifti_task)
        
        % This numbers that variants, clustering adjacent vertices that =1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf])
        %system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        % writes file with numbered clusters in preparation for size
        % exclusion
        cifti_rest = ft_read_cifti_mod(outfilewbrest);
        %cifti_task = ft_read_cifti_mod(outfilewbtask);

        if ExcludeBySize == 1             
            % exclusion criteria is set to 15 vertices (any variant less
            % than 15 vertices big will be excluded
            [cifti_rest.data] = ExcludeVariantbySize(cifti_rest.data, subject, threshold(x), 50);
        end 
        
        % writes size-excluded variant masks
        ft_write_cifti_mod(outfilerest, cifti_rest)
        %ft_write_cifti_mod(outfiletask, cifti_task)
    end

function [cifti_rest_data cifti_task_data] = ExcludeVariantbySize(cifti_rest_data, subject, threshold, exclusion_criteria)

        rest_sizes = [];
        
        % same values as in cifti_rest.data but w/o repetitions, of
        % vertices?
        vars_rest = unique(cifti_rest_data);
        
        
        %counts vertices to determine if variant meet exclusion criteria
        for q = 1:length(vars_rest)
            vertcount = 0;
            for r = 1:length(cifti_rest_data)
                if cifti_rest_data(r) == vars_rest(q)
                    vertcount = vertcount + 1;
                end
            end
            rest_sizes = [rest_sizes; vars_rest(q) vertcount];
        end

        
        %removing vertices belonging to variants that are not at least 15
        %vertices big
        for i = 2:size(rest_sizes,1)            
            if rest_sizes(i,2) < exclusion_criteria                
                removeverts = find(cifti_rest_data == rest_sizes(i,1));                
                cifti_rest_data(removeverts,1) = 0;                
%             else                
%                 setverts = find(cifti_rest_data == rest_sizes(i,1));                
%                 cifti_rest_data(setverts,1) = 1;                
            end
        end
        

 end