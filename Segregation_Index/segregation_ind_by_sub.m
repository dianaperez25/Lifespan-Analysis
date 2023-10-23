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
%% OPTIONS
% ------------------------------------------------------------------------
dataset = 'lifespan';
atlas = 'Parcels333';
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 988; % if this is commented out or set to 0, then the script will calculate it

% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
data_dir = '/Volumes/Back_Up/Dissertation/parcellations/FC_Parcels_333';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
ls_subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
inet_subject = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %
ls_sessions = 5;
inet_sessions = 4;

% ------------------------------------------------------------------------
%% DATA MATCHING
% ------------------------------------------------------------------------
if match_data
    if amt_data == 0 || ~exist('amt_data')
        allSubs_amtData = [];
        % get minimum amt of data per session
        subject = [ls_subject inet_subject];
        allSubs_amtData = [];
        for sub = 1:numel(subject)
            if contains(subject{sub}, 'LS')
                sessions = ls_sessions;
            elseif contains(subject{sub}, 'INET')
                sessions = inet_sessions;
            else error('Invalid subject ID');
            end
            for ses = 1:sessions
                load([data_dir '/sub-' subject{sub} '_rest_ses-' num2str(ses) '_parcel_timecourse.mat'])
                masked_data = parcel_time(logical(tmask_concat),:)';
                allSubs_amtData = [allSubs_amtData; size(masked_data,2)];
                if size(masked_data,2) < 800
                    disp(sprintf('subject %s session %d has %d data points', subject{sub}, ses, size(masked_data,2)))
                end
            end
        end
        amt_data = min(min(allSubs_amtData));
    end
end

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
% get correct subject and session info
if strcmpi(dataset, 'lifespan')
    subs = ls_subject;
    sessions = ls_sessions;
elseif strcmpi(dataset, 'inetworks')
    subs = inet_subject;
    sessions = inet_sessions;
else
    error('Invalid Dataset: please select either lifespan or inetworks')
end

% get atlas information
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
for net = 1:13
    net_size(net) = length(atlas_params.mods{1,net});
end
num_parcs = atlas_params.num_rois; 
weights = net_size./(num_parcs-net_size(1));

% begin segregation analysis calculations
for sub = 1:numel(subs)
    sub_struct = {};
    all_within = [];
    all_between = [];

    clear subcorrmat mean_matrix

    for ses = 1:sessions
        fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_timecourse.mat', data_dir, subs{sub}, ses);
        timecourse_struct = load(fname);
        masked_data = timecourse_struct.parcel_time(logical(timecourse_struct.tmask_concat), :)';
        sorted_data = masked_data(atlas_params.sorti, :);
        clear timecourse_struct

        % ... match amount of data ...
        if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
            amt_data = size(sorted_data,2);
        end
        matched_data = datasample(sorted_data,amt_data,2,'Replace', false);

        % ... calculate the correlation matrix ...
        matrix_sorted = paircorr_mod(matched_data');
        matrix_sorted = single(FisherTransform(matrix_sorted));% fisher transform r values
        sub_corrmat(ses,:,:,:) = matrix_sorted;
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
        within(within<0) = 0;
        all_within = [all_within; within];
        
        tmp_between = mean_matrix(rois(1):rois(end), 1:num_parcs); % all network correlations
        maskmat = logical(ones(size(tmp_between))); % mask out within-network correlations
        maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
        between = tmp_between(maskmat==1); %between-network correlations
        between(between<0) = 0;
        all_between = [all_between;between];

        count = count + net_size(net);
    end
        
    %% calculate the segregation index by network by session
    sub_SI(sub) = (mean(all_within) - mean(all_between))/mean(all_within);                      
end
save(sprintf('%s/%s_allsubs_seg_index_sub_matched.mat', output_dir, dataset), 'sub_SI')

