%% SEGREGATION INDEX FOR EACH PARCEL (from individualized parcellation)
% This script calculates a measure of system segregation.
% Input: individualized parcellation, parcel average timeseries, parcel
% network assignment
% Output: mean within-system correlation, mean between-system correlation,
% a segregation index for parcel and one segregation index per subject (averaged across parcels).
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
parcel_timecourse_dir = '/scratch/dcr8536/parcellations/avg_timecourses/';
network_map_dir = '/scratch/dcr8536/parcellations/'
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};

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
    
    %load parcels average timecourse
    %% DP: CHECK IF IT'S ALREADY MASKED
    parcel_tc_fname = sprintf('%s/sub-%s_individual_parcels_average_timecourses.mat', parcel_timecourse_dir, subs{sub});
    mat_struct = load(parcel_tc_fname);
        
    %% DP: DO I HAVE TO EXTRACT A SPECIFIC VARIABLE?
    
    %load network map
    netmap_fname = sprintf('%s/sub-%s/sub-%s_indiv_parcels_net_assigned.dtseries.nii', network_map_dir, subs{sub}, subs{sub});
    netmap = ft_read_cifti_mod(netmap_fname);
    
    % first sort by network
    unique_net_IDs = unique(netmap.data);
    if unique_net_IDs(1) == 0
        unique_net_IDs(1) = [];
    end
    net_size = [];
    parcels_sorted = [];
    for net = 1:length(unique_net_IDs)
        % fing the parcels that are assigned to each network
        net_size(net,1) = unique_net_IDs(net); 
        net_parcels = find(netmap.data==unique_net_IDs(net));
        net_size(net,2) = length(net_parcels);
        parcels_sorted = [parcels_sorted; net_parcels];
    end
    num_nodes = length(parcels_sorted);
    parc_tc_sorted = mat_struct(parcels_sorted);
    matrix_sorted = single(FisherTransform(paircorr_mod(parc_tc_sorted')));% fisher transform r values
        
        
        % Separate ROIs by Networks   
        % initialize some variables
        within = {}; %this will contain all the correlation values for the rois in each system, only for checking, will delete later
        between = {}; %this will contain all the correlation values for the rois outside of each system, only for checking, will delete later
        
        %num_networks = size(atlas_params.networks,1); % extract number of networks
        count = 1;
        for net = 2:size(net_size,2) % go through each network
            net_parcels = count:(count+net_size(net,2)-1); %extract the parcels belonging to system n            
            weights(net) = net_size(net,2)/length(parcels_sorted);
            %% GET WITHIN SYSTEM CORRELATIONS
            tmp = matrix_sorted(net_parcels(1):net_parcels(end),net_parcels(1):net_parcels(end)); % extract correlation values within system n
            %make a mask to only keep correlation values in upper triangle (because
            %corrmats are symmetrical)
            maskmat = ones(net_size(net,2)); 
            maskmat = logical(triu(maskmat,1));
            within = tmp(maskmat); % mask out values in lower triangle
            within(within<0) = 0; % set negative z values to 0            
            % put those z values in a structure to check/use later
            %means_within(net) = mean(within); 
            weighted_means_within(net) = weights(net) * mean(within); 
            within_network{net,sub} = within; %maybe delete this
            within_corrs = [within_corrs; within];
            clear maskmat tmp % clear up some variables
            
            %% GET BETWEEN SYSTEM CORRELATIONS
            tmp = matrix(net_parcels(1):net_parcels(end), 1:num_nodes); % extract z values for nth network
            maskmat = ones(size(tmp)); % make another
            maskmat(:,net_parcels(1):net_parcels(end)) = 0; % mask out z values within system n
            between = tmp(maskmat==1); % put only between system z values in a variable
            between(between<0) = 0;% set negative z values to 0
            between_corrs = [between_corrs; between];
            %put those z values in a structure to check/use later
            %means_between(net) = mean(between); 
            weighted_means_between(net) = weights(net) * mean(between);
            between_network{net,sub} = between;

            clear maskmat tmp node_mean_within node_mean_between node_SI% clear up some variables
            %% calculate the segregation index by network by session
            seg_index_by_net(sub,net) = (mean(within) - mean(between))/mean(within);
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
