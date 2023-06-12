%% Network Segregation Script - Individualized Parcels
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
addpath(genpath('/Users/dianaperez/Desktop/Dependencies'))
% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/parcellations/';
output_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_Analyses/indiv_parcs/';
if ~exist(output_dir)
    mkdir(output_dir)
end
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
parcellation = 'Gordon333';

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
subjects_seg_index = [];
for sub = 1:numel(subjects)
    sub_struct = []; % a structure that will contain segregation indices for each network, and all networks at the end
    all_within = []; % initialize variable that will contain all within-networks correlations across networks 
    all_between = []; % initialize variable that will contain all between-networks correlations across networks 

    fname = sprintf('%s/%s_Parcellation/avg_timecourses/corrmats/sub-%s_indiv_parcels_corrmats_sorted.mat', data_dir, parcellation, subjects{sub});
    mat_struct = load(fname); % load subjects data
    sorted_matrix = single(FisherTransform(mat_struct.sorted_corrmat));% fisher transform r values
    networks = unique(mat_struct.parcel_netid(:,2)); % collect network ids represented in subjects data
    num_rois = length(mat_struct.parcel_netid); % collect number of parcels for this subject
    
    count = 1; % initialize a count that will change according to number of parcels
    
    for net = 1:size(networks,1) % for each network...
        net_size = length(find(mat_struct.parcel_netid(:,2)==networks(net))); % we want the number of parcels assigned to the network...
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
        
        % we put the subjects data into a structure
        sub_struct(net,1) = networks(net); %this is the network id 
        sub_struct(net,2) = (mean(within) - mean(between))/mean(within); %and the networks segregation index
        count = count + net_size; % then we move on to the next network
        clear within between % making sure to clear some variables first
    end
        %% calculate the segregation index by network by session
        subjects_seg_index(sub) = (mean(all_within) - mean(all_between))/mean(all_within);
        sub_struct(end+1,1) = 'A'; %this is the network id 
        sub_struct(end,2) = subjects_seg_index(sub);
        save(sprintf('%s/sub-%s_seg_index_by_net_%s.mat',output_dir,subjects{sub}, parcellation), 'sub_struct')

end
save(sprintf('%s/Lifespan_seg_index_by_sub_%s.mat', output_dir, parcellation), 'subjects_seg_index')

