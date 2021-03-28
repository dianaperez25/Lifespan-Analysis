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
% PATHS
% ------------------------------------------------------------------------
dataLoc = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
atlas_dir = '/Users/diana/Desktop/Research/Atlases/';
atlas = 'Seitzman300';
% load atlas that contains roi info (including which rois belong to each network) 
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);

% ------------------------------------------------------------------------
% VARIABLES
% ------------------------------------------------------------------------

subs = {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS10'};
sessions = [3, 5, 1, 5, 2, 1];
network_z.within = [];
network_z.between = [];
network_z.networks = atlas_params.networks; %copy names of each network
seg_ind = [];

% ------------------------------------------------------------------------
% MAIN FOR-LOOP
% ------------------------------------------------------------------------

for sub = 1:numel(subs)
    for ses = 1:sessions(sub)
        % ------------------------------------------------------------------------
        % STEP 1: Load Correlation Matrix
        % ------------------------------------------------------------------------
        %loads entire structure saved by FCProc script
        fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', dataLoc, subs{sub}, subs{sub}, ses);
        mat_struct = load(fname);
        %extracts only correlation matrix
        matrix_orig = mat_struct.corrmat;
        clear mat_struct
        % sort matrix into correct order (so that rois that belong to same system
        % are grouped together)
        matrix_sorted = matrix_orig(atlas_params.sorti,atlas_params.sorti); %% this might be an extra step since the matrix might already be stored after sorting in the correct order, but just to be sure, we'll do it again
        % fisher transform r values
        matrix = single(FisherTransform(matrix_sorted));
        clear matrix_orig matrix_sorted
        
        % ------------------------------------------------------------------------
        % STEP 2: Separate ROIs by Networks
        % ------------------------------------------------------------------------
        
        % initialize some variables
        within = {}; %this will contain all the correlation values for the rois in each system, only for checking, will delete later
        between = {}; %this will contain all the correlation values for the rois outside of each system, only for checking, will delete later

        %num_networks = size(atlas_params.networks,1); % extract number of networks
        for n = 1:size(atlas_params.networks,1) % go through each network
            rois = atlas_params.mods{1,n}; %extract the rois belonging to system n
            num_rois = length(rois); % number of rois in system n
            
            %% GET WITHIN SYSTEM CORRELATIONS
            tmp = matrix(rois(1):rois(end),rois(1):rois(end)); % extract correlation values within system n
            
            %make a mask to only keep correlation values in upper triangle (because
            %corrmats are symmetrical)
            maskmat = ones(num_rois); 
            maskmat = logical(triu(maskmat,1));
            within = tmp(maskmat); % mask out values in lower triangle
            
            % set negative z values to 0
            for z = 1:length(within)
                if within(z) < 0
                    within(z) = 0;
                end
            end
            
            % put those z values in a structure to check/use later
            means_within(n) = mean(within); 
            within_network{ses,n,sub} = within;
            clear maskmat within tmp % clear up some variables

            %% GET BETWEEN SYSTEM CORRELATIONS
            tmp = matrix(rois(1):rois(end), 1:300); % extract z values for nth network
            maskmat = ones(size(tmp)); % make another
            maskmat(:,rois(1):rois(end)) = 0; % mask out z values within system n
            between = tmp(maskmat==1); % put only between system z values in a variable
            
            % set negative z values to 0
            for z = 1:length(between)
                if between(z) < 0
                    between(z) = 0;
                end
            end
            
            %put those z values in a structure to check/use later
            means_between(n) = mean(between); 
            between_network{ses,n,sub} = between;
            clear maskmat between tmp % clear up some variables
        end
        
        network_z.within{sub,ses} = means_within; %put mean within system corrs in main structure
        network_z.between{sub,ses} = means_between; %put mean between system corrs in main structure
        
        %calculate segregation index = (mean within system Z - mean between system Z)/mean within system Z
        seg_ind(sub,ses) = (mean(means_within)-mean(means_between))/(mean(means_within));
    end
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