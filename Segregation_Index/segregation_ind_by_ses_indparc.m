%% Network Segregation Script by Session - Individualized Parcels
% This script calculates a measure of system segregation (how much nodes in
% a system communicate with each other vs with nodes in other systems).
% Input: 300x300 functional connectivity correlation matrix.
% Output: mean within-system correlation, mean between-system correlation,
% a segregation index for each session.
% Based on the analysis describe in: 
% Chan, Micaela Y., et al. "Decreased segregation of brain systems across 
% the healthy adult lifespan." Proceedings of the National Academy of 
% Sciences 111.46 (2014): E4997-E5006.
%
% NOTES: in original analysis negative z-values were set to zero.
% Within-system connectivity was calculated as the mean node-to-node 
% z-value of all nodes of that system to each other. 
% Between-system connectivity was calculated as the mean node-to-node 
% z-value between each node of a system and all nodes of all other systems. 
% ------------------------------------------------------------------------
clear all
addpath(genpath('/Users/dianaperez/Documents/Dependencies'))
% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/parcellations/';
surf_ts_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/';
tmask_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
output_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_Analyses/indiv_parcs/';
if ~exist(output_dir)
    mkdir(output_dir)
end
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'}; %
parcellation = 'Gordon';
fname_stem = 'individual_parcels_edgethresh_0.5.dtseries';
%fname_stem = 'cMSHBM';
sessions = 5;
max_runs = 13;
% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
ses_seg_index = [];
for sub = 1:numel(subjects)
    %cifti containing parcel id for each vertex
%     parcel_id = ft_read_cifti_mod(sprintf('%s/%s_Parcellation/sub-%s/sub-%s_%s.dtseries.nii', data_dir, parcellation, subjects{sub}, subjects{sub}, fname_stem));
%     parcel_id = parcel_id.data;
%     parcels = unique(parcel_id);
%     if parcels(1) == 0
%         parcels(1) = [];
%     end
    parc_fname = sprintf('%s/%s_Parcellation/avg_timecourses/corrmats/sub-%s_indiv_parcels_corrmats_sorted.mat', data_dir, parcellation, subjects{sub});
    mat_struct = load(parc_fname); % load subjects data
    timeseries = mat_struct.sorted_ts;
    parcel_netid = mat_struct.parcel_netid;
    clear mat_struct %to save memory
    for ses = 1:sessions
        %sub_struct = []; % a structure that will contain segregation indices for each network, and all networks at the end
        all_within = []; % initialize variable that will contain all within-networks correlations across networks 
        all_between = []; % initialize variable that will contain all between-networks correlations across networks 
        ses_fname=(sprintf('%s/%s_Parcellation/avg_timecourses/sub-%s_ses-%d_indiv_parcels_timseries_sorted.mat', data_dir, parcellation, subjects{sub}, ses));

        
        if ~exist(ses_fname) %if data by session, doesn't exist, make the files
            ses_tmask = [];    
            %concat_tseries = [];
            for run =  1:max_runs     
                %we're going to use this variable just to see if the run exists             
                surf_ts_fname = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',surf_ts_dir, subjects{sub}, ses, subjects{sub}, ses, run);
                if exist(surf_ts_fname)
                    %surf_tseries = ft_read_cifti_mod(surf_ts_fname);
                    tmask_fname = sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%d_desc-tmask_fFD.txt',tmask_dir,subjects{sub},ses,subjects{sub},ses,run);
                    tmask = table2array(readtable(tmask_fname)); 
                    %surf_tseries = surf_tseries.data(:,logical(tmask));
                    %concat_tseries = [concat_tseries surf_tseries];
                    ses_tmask = [ses_tmask tmask'];
                end
            end
            
            
            ses_inds = 1:sum(ses_tmask);
            ses_tseries = timeseries(:,ses_inds);
            timeseries(:,ses_inds) = [];
            save(ses_fname, 'ses_tseries')
        else
            ses_tseries = load(ses_fname);
            ses_tseries = ses_tseries.ses_tseries;
            clear timeseries
        end
        
                
        %load the structure with the sorted parcel info        
        %sorted_parcels = mat_struct.sorted_parcels;
        
        % sort parcels by network
        %sorted_tseries = ses_avg_tseries(sorted_parcels,:);
        % now we can create a sortex corrmat
        sorted_matrix = paircorr_mod(ses_tseries');
        % and we Fisher transform it
        sorted_matrix = single(FisherTransform(sorted_matrix));% fisher transform r values
        % now we have to figure out the network labels...
        networks = unique(parcel_netid(:,2)); % collect network ids represented in subjects data
        num_rois = length(parcel_netid); % collect number of parcels for this subject
            
        count = 1; % initialize a count that will change according to number of parcels
    
        for net = 1:size(networks,1) % for each network...
            net_size = length(find(parcel_netid(:,2)==networks(net))); % we want the number of parcels assigned to the network...
            if net == 1
                rois = [count:net_size(net)]; % we want to know which parcels are assigned to the network...
            else
                rois = [count:(count-1) + net_size];
            end

            tmp_within = sorted_matrix(rois,rois); % we want the within-network correlations...
            maskmat = ones(size(tmp_within)); % but only the upper triangle...
            maskmat = logical(triu(maskmat,1));
            within = tmp_within(maskmat);
            within(within<0) = 0; % and we want to set negative correlations to zero...
            all_within = [all_within; within]; % then we put them together with other within-network correlations for other systems for later...

            % then we do the same for between-network correlations...
            tmp_between = sorted_matrix(rois(1):rois(end), 1:num_rois); % all network correlations
            maskmat = ones(size(tmp_between)); % mask out within-network correlations
            maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
            between = tmp_between(maskmat==1); % because we only want between-network correlations
            between(between<0) = 0; % set negatives to zero
            all_between = [all_between;between];

            count = count + net_size; % then we move on to the next network
            clear within between % making sure to clear some variables first
        end
        %% calculate the segregation index by network by session
        ses_seg_index(sub,ses) = (mean(all_within) - mean(all_between))/mean(all_within);

    end
end
save(sprintf('%s/Lifespan_seg_index_by_ses_%s.mat', output_dir, parcellation), 'ses_seg_index')

