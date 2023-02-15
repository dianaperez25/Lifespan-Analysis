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
%dataLoc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
output_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';
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
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
%subs = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
sessions = [5, 5, 5, 5, 5, 5, 5, 5];
%sessions = [4,4,4,4,4,4,4,4,4,4];

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
for sub = 1:numel(subs)
    sub_struct = {};
    all_within = [];
    all_between = [];
    %sub_SI = [];
    clear subcorrmat
    clear mean_matrix
    for ses = 1:sessions(sub)
        fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', dataLoc, subs{sub}, subs{sub}, ses);
        mat_struct = load(fname);
        matrix = mat_struct.corrmat; clear mat_struct
        matrix = single(FisherTransform(matrix));% fisher transform r values
        sub_corrmat(ses,:,:,:) = matrix;
    end
        mean_matrix(1,:,:,:) = squeeze(mean(sub_corrmat));
        if ndims(mean_matrix)>2
            mean_matrix = squeeze(mean_matrix);
        end
        clear sub_corrmat
        count = 1;
        for net = 1:size(atlas_params.networks,1)
            if net == 1
                rois = [count:net_size(net)];
            else
                rois = [count:(count-1) + net_size(net)]; %extract the rois belonging to system n
            end
            
            tmp_within = mean_matrix(rois,rois); % within-network correlations
            maskmat = ones(size(tmp_within));
            maskmat = logical(triu(maskmat,1));
            within = tmp_within(maskmat);
            %within(within<0) = [];
            within(within<0) = 0;
            all_within = [all_within; within];
            %all_within(net,ses) = all_within;
            
            tmp_between = mean_matrix(rois(1):rois(end), 1:300); % all network correlations
            maskmat = ones(size(tmp_between)); % mask out within-network correlations
            maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
            between = tmp_between(maskmat==1); %between-network correlations
            %between(between<0) = [];
            between(between<0) = 0;
            all_between = [all_between;between];

            count = count + net_size(net);
        end
        %% calculate the segregation index by network by session
        sub_SI(sub) = (mean(all_within) - mean(all_between))/mean(all_within);                      
end
save([output_dir 'LS_allsubs_seg_index_sub_with0neg.mat'], 'sub_SI')

