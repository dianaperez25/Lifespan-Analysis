%% SEGREGATION INDEX
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

% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
dataLoc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas = 'Parcels333';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 

% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = [5, 5, 5, 5, 5, 5, 5, 5];
num_nodes = 333;
network_z.within = [];
network_z.between = [];
network_z.networks = atlas_params.networks{2:end}; %copy names of each network
%network_z.networks{1} = [];
seg_ind = [];
within_corrs = [];
between_corrs = [];


% ------------------------------------------------------------------------
%% MAIN FOR-LOOP
% ------------------------------------------------------------------------
for sub = 1:numel(subs)
    for ses = 1:sessions(sub)
        %load correlation matrix file
        fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', dataLoc, subs{sub}, ses);
        %fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', dataLoc, subs{sub}, subs{sub}, ses);
        mat_struct = load(fname);
        %extracts only correlation matrix
        matrix_orig = mat_struct.parcel_corrmat;
        clear mat_struct
        % sort matrix into correct order (so that rois that belong to same system
        % are grouped together)
        matrix_sorted = matrix_orig(atlas_params.sorti,atlas_params.sorti); %% this might be an extra step since the matrix might already be stored after sorting in the correct order, but just to be sure, we'll do it again        
        matrix = single(FisherTransform(matrix_sorted));% fisher transform r values
        sub_corrmat(ses,:,:,:) = matrix;
        clear matrix_orig matrix_sorted
        
        % Separate ROIs by Networks   
        % initialize some variables
        within = {}; %this will contain all the correlation values for the rois in each system, only for checking, will delete later
        between = {}; %this will contain all the correlation values for the rois outside of each system, only for checking, will delete later
        node_within = zeros(num_nodes,1);
        node_between = zeros(num_nodes,1);
        node_seg_index = zeros(num_nodes,1);
        nodes_by_net = {};
        
        %num_networks = size(atlas_params.networks,1); % extract number of networks
        for net = 2:size(atlas_params.networks,2) % go through each network
            rois = atlas_params.mods{1,net}; %extract the rois belonging to system n
            %num_rois = length(rois); % number of rois in system n
            size_nets(net) = length(rois);
            weights(net) = size_nets(net)/num_nodes;
            %% GET WITHIN SYSTEM CORRELATIONS
            tmp = matrix(rois(1):rois(end),rois(1):rois(end)); % extract correlation values within system n
            %make a mask to only keep correlation values in upper triangle (because
            %corrmats are symmetrical)
            maskmat = ones(size_nets(net)); 
            maskmat = logical(triu(maskmat,1));
            within = tmp(maskmat); % mask out values in lower triangle
            within(within<0) = 0; % set negative z values to 0            
            % put those z values in a structure to check/use later
            %means_within(net) = mean(within); 
            weighted_means_within(net) = weights(net) * mean(within); 
            within_network{ses,net,sub} = within; %maybe delete this
            within_corrs = [within_corrs; within];
            % get within system correlations by node
            for roi = 1:size_nets(net)
                node = rois(roi);
                node_values = tmp(roi,:);
                node_values(roi) = []; 
                node_values(node_values<0) = 0;
                node_mean_within(roi) = mean(node_values);
                node_within(node) = node_mean_within(roi);
            end                    
            nodes_by_net{net,1} = node_mean_within;
            
            clear maskmat tmp % clear up some variables
            
            %% GET BETWEEN SYSTEM CORRELATIONS
            tmp = matrix(rois(1):rois(end), 1:num_nodes); % extract z values for nth network
            maskmat = ones(size(tmp)); % make another
            maskmat(:,rois(1):rois(end)) = 0; % mask out z values within system n
            between = tmp(maskmat==1); % put only between system z values in a variable
            between(between<0) = [];% set negative z values to 0
            between_corrs = [between_corrs; between];
            %put those z values in a structure to check/use later
            %means_between(net) = mean(between); 
            weighted_means_between(net) = weights(net) * mean(between);
            between_network{ses,net,sub} = between;
            
            % get between system correlations by node
            tmp(:,rois(1):rois(end)) = [];
            for roi = 1:size_nets(net)
                node = rois(roi);
                node_values = tmp(roi,:); 
                node_values(node_values<0) = 0;
                node_mean_between(roi) = mean(node_values);
                node_between(node) = node_mean_between(roi);
                node_SI(roi) = (node_mean_within(roi) - node_mean_between(roi))/node_mean_within(roi);
                node_seg_index(node) = node_SI(roi);
            end                    
            nodes_by_net{net,2} = node_mean_between;
            nodes_by_net{net,3} = node_SI;
            
            clear maskmat tmp node_mean_within node_mean_between node_SI% clear up some variables
            %% calculate the segregation index by network by session
            seg_index_by_net(sub,ses,net) = (mean(within) - mean(between))/mean(within);
        end
        
        %% calculate segregation index by session across all networks (maybe exclude unassigned network (1))
        seg_index_by_ses(sub,ses) = (mean(weighted_means_within) - mean(weighted_means_between))/mean(weighted_means_within);
        seg_index_by_ses2(sub,ses) = (mean(within_corrs) - mean(between_corrs))/mean(within_corrs);
        %network_z.within{sub,ses} = means_within; %put mean within system corrs in main structure
        %network_z.between{sub,ses} = means_between; %put mean between system corrs in main structure
        
        %calculate segregation index = (mean within system Z - mean between system Z)/mean within system Z
        %seg_ind(sub,ses) = (mean(means_within)-mean(means_between))/(mean(means_within));
%         out_mat = [node_within node_between node_seg_index];
%         out_dir = '/Users/dianaperez/Desktop/';
%         outfile1 = sprintf('%s/sub-%s_ses-%d_segregation_indices_300rois.mat', out_dir, subs{sub}, ses);
%         save(outfile1, 'out_mat');
%         outfile2 = sprintf('%s/sub-%s_ses-%d_node_segInds_byNet.mat', out_dir, subs{sub}, ses);
%         save(outfile2, 'nodes_by_net');
    end
    %% calculate segregation index per subject (across all sessions)
    % either average across all segregation indices across sessions
    sub_avg = squeeze(mean(sub_corrmat));
    % or average the matrices and calculate segregation index
end

% calculate seg_ind by network?
tmp = [];
segInd_byNet = [];
for sub = 1:numel(subs)
    for ses = 1:sessions(sub)
        for net = 1:size(atlas_params.networks,1)
            within_sys = network_z.within{sub,ses}(net);
            between_sys = network_z.between{sub,ses}(net);
            tmp(net,ses) = (within_sys - between_sys)/within_sys;
        end
    end
    segInd_byNet = cat(2,segInd_byNet,tmp);
    clear tmp
end