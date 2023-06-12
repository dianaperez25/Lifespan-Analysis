%% Script to make average timecourse for parcels
clear all
addpath(genpath('/projects/b1081/Scripts/Resources/'))
addpath(genpath('/Users/dianaperez/Desktop/Dependencies'))

% PATHS
ind_parc_loc =  '/scratch/dcr8536/parcellations/subsample_10/';
surf_ts_loc = '/scratch/dcr8536/TimeB/Nifti/postFCproc_CIFTI/';
output_dir = '/scratch/dcr8536/parcellations/avg_timecourses/';
tmask_loc = '/scratch/dcr8536/TimeB/Nifti/preproc_fmriprep-20.2.0/fmriprep/';
% VARIABLES
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = 5;
max_runs = 11;

if ~exist(output_dir)
    mkdir(output_dir)
end

% step 1: load the individualized parcellation file
for sub = 1:numel(subject)
    % load the cifti file with individual parcels for each subject
    ind_parc_fname = sprintf('%s/sub-%s/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', ind_parc_loc, subject{sub}, subject{sub});
    ind_parc = ft_read_cifti_mod(ind_parc_fname);
    surf_ts_concat = [];
    avg_parc_ts = [];
    for ses = 1:sessions
        for run = 1:max_runs
            % step 2: load the timecourse by vertex file for each session
            % and run
            surf_ts_fname = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii', surf_ts_loc, subject{sub}, ses, subject{sub}, ses, run);
            if exist(surf_ts_fname)
                surf_ts = ft_read_cifti_mod(surf_ts_fname);
                tmask_fname = sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%d_desc-tmask_fFD.txt',tmask_loc,subject{sub},ses,subject{sub},ses,run);
                tmask = table2array(readtable(tmask_fname)); 
                surf_ts = surf_ts.data(:,logical(tmask));
                surf_ts_concat = [surf_ts_concat surf_ts];
            end
        end
    end                         
    % step 3: index the vertices for each parcel; average their timecourses
    % together
    unique_ind_parc = unique(ind_parc.data);
    if unique_ind_parc(1) == 0
        unique_ind_parc(1) = [];
    end
    count = 1;
    for parc = 1:length(unique_ind_parc)
        inds = find(ind_parc.data==unique_ind_parc(parc));
        parc_ts = surf_ts_concat(inds, :);
        avg_parc_ts(count,:) = mean(parc_ts);
        count = count + 1;
        clear inds parc_ts
    end
    
    % step 4: save the file
    save([output_dir 'sub-' subject{sub} '_individual_parcels_average_timecourses.mat'], 'avg_parc_ts', '-v7.3');
    clear avg_parc_ts 
end