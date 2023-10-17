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
dataLoc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/iNetworks/Nifti/FC_Parcels_333/';
%dataLoc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
output_dir = '/Users/dianaperez/Desktop/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas = 'Parcels333';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
num_nets = length(atlas_params.networks);
for net = 1:num_nets
    net_size(net) = length(atlas_params.mods{1, net});
end
num_rois = atlas_params.num_rois;
weights = net_size./(num_rois-net_size(1));
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
%subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
subs = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %
%sessions = [5, 5, 5, 5, 5, 5, 5, 5];
sessions = 4;

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
    for ses = 1:sessions
        fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', dataLoc, subs{sub}, ses);
        mat_struct = load(fname);
        matrix = mat_struct.parcel_corrmat; clear mat_struct
        matrix = single(FisherTransform(matrix));% fisher transform r values
        sub_corrmat(ses,:,:,:) = matrix;
    end
        mean_matrix(1,:,:,:) = squeeze(mean(sub_corrmat));
        if ndims(mean_matrix)>2
            mean_matrix = squeeze(mean_matrix);
        end
        clear sub_corrmat
        count = 1;
        for net = 1:size(atlas_params.networks,2)
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
            
            tmp_between = mean_matrix(rois(1):rois(end), 1:333); % all network correlations
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

