%% New Network Segregation Script
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
dataLoc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
output_dir = '/Users/dianaperez/Desktop/Segregation_Analyses/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas = 'Seitzman300';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
for net = 1:14
    net_size(net) = length(atlas_params.mods{1,net});
end
weights = net_size./(300-net_size(1));
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = [5, 5, 5, 5, 5, 5, 5, 5];

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
for sub = 1:numel(subs)
    sub_struct = {};
    node_mean_within = [];
    node_mean_between = [];
    node_SI = [];
    for ses = 1:sessions(sub)
        fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', dataLoc, subs{sub}, subs{sub}, ses);
        mat_struct = load(fname);
        matrix = mat_struct.corrmat; clear mat_struct
        matrix = single(FisherTransform(matrix));% fisher transform r values
        %sub_corrmat(ses,:,:,:) = matrix;
        count = 1;
        for net = 1:size(atlas_params.networks,1)
            if net == 1
                rois = [count:net_size(net)];
            else
                rois = [count:(count-1) + net_size(net)]; %extract the rois belonging to system n
            end
            tmp_within = matrix(rois,rois); % within-network correlations
            tmp_between = matrix(rois(1):rois(end), 1:300); % all network correlations
            maskmat = ones(size(tmp_between)); % mask out within-network correlations
            maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
            between = tmp_between(maskmat==1); %between-network correlations
            between = reshape(between,length(rois), 300-length(rois)); %reshape matrix 
            means_within = []; % initialize variable containing 
            means_between = [];
            for roi = 1:length(rois)
                node_within = tmp_within(roi,:); % all within-network corrs for the node roi
                node_within(roi) = []; % delete correlation with itself
                node_within(node_within<0) = []; % delete negative correlations
                %within_corrs = [within_corrs; node_within']; % variable containing all within-network correlations for each session
                node_mean_within(count,ses) = mean(node_within); %mean within-network correlations for node roi
                %means_within = [means_within; node_mean_within(count,ses)]; %variable containing the mean within-network corrs for all nodes in a network
                node_between = between(roi,:); %all between-network corrs for node roi
                node_between(node_between<0) = []; % delete negative correlations
                %between_corrs = [between_corrs; node_between']; % variable containing all between-network correlations for each session
                node_mean_between(count,ses) = mean(node_between); %mean between-network corrs for node roi
                %means_between = [means_between; node_mean_between(count,ses)];%variable containing the mean between-network corrs for all nodes in a network
                node_SI(count,ses) = (node_mean_within(count,ses) - node_mean_between(count,ses))/node_mean_within(count,ses);
                count = count + 1;
            end
        
            %weighted_means_within(net) = weights(net) * mean(means_within);             
            %weighted_means_between(net) = weights(net) * mean(means_between);
            %means_between_all = [means_between_all; means_between];
            %means_within_all = [means_within_all; means_within];
            %% calculate the segregation index by network by session
            %seg_index_net(net,ses) = (mean(means_within) - mean(means_between))/mean(means_within);
            %within_net(net,ses) = mean(means_within);
            %between_net(net,ses) = mean(means_between);
        end
        %seg_index_ses(ses) = (mean(node_mean_within(:,ses))-mean(node_mean_between(:,ses)))/mean(node_mean_within(:,ses));
    end
         sub_struct.within = node_mean_within;
         sub_struct.between = node_mean_between;
         sub_struct.seg_ind = node_SI;
         save([output_dir subs{sub} '_seg_index_node.mat'], 'sub_struct')
%         sub_struct.network.within = within_net;
%         sub_struct.network.between = between_net;
%         sub_struct.network.seg_ind = seg_index_net;
%         sub_struct.session.within = mean(node_mean_within,1);
%         sub_struct.session.between = mean(node_mean_between,1);
%         sub_struct.session.seg_ind = seg_index_ses;
        %sub_struct.subject.within = 
        %sub_struct.subject.between = 
        %sub_struct.subject.seg_ind = 
end
