%% Calculate Segregation Index by Session
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
%data_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/FC_Parcels_333/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/';
output_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';
if ~exist(output_dir)
    mkdir(output_dir)
end
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';

% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
LS_subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
iNet_subject = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
LS_sessions = 1;
iNet_sessions = 1;

% ------------------------------------------------------------------------
%% OPTIONS
% ------------------------------------------------------------------------
neg_corrs = 'zero'; % choose: 'nan', 'zero', 'asis'. How to deal with negative correlations
data_set = 'Lifespan'; %'Lifespan' or 'iNetworks' 
atlas = 'Parcels333';
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 968; % if this is commented out or set to 0, then the script will calculate it

atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
for net = 1:length(atlas_params.networks)
    net_size(net) = length(atlas_params.mods{1,net});
end
weights = net_size./(300-net_size(1));

%% DATA MATCHING
if match_data
    if amt_data == 0
        allSubs_amtData = [];
        % get minimum amt of data per session
        subject = [LS_subject iNet_subject];
        allSubs_amtData = [];
        for s = 1:numel(subject)
            if contains(subject{s}, 'LS')
                LS_sessions = LS_sessions;
            elseif contains(subject{s}, 'INET')
                LS_sessions = iNet_sessions;
            else error('Invalid subject ID');
            end

            for i = 1:LS_sessions
                if strcmpi(atlas, 'Seitzman300')            
                    load([data_dir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
                    masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
                elseif strcmpi(atlas, 'Parcels333')
                    load([data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
                    masked_data = parcel_time(logical(tmask_concat),:)';
                end
                if match_data
                    allSubs_amtData = [allSubs_amtData; size(masked_data,2)];
                end
            end
        end
        amt_data = min(min(allSubs_amtData));
    end
end

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
if strcmpi(data_set, 'Lifespan')
    subject = LS_subject;
    sessions = LS_sessions;
elseif strcmpi(data_set, 'iNetworks')
    subject = iNet_subject;
    sessions = iNet_sessions;
end

for s = 1:numel(subject)
    sub_struct = {};
    all_within = [];
    all_between = [];
    %ses_SI = [];
    for i = 1:sessions
        if strcmpi(atlas, 'Seitzman300')
            load([data_dir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
            masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
            ind = 1;
            num_nodes = 300;
        elseif strcmpi(atlas, 'Parcels333')
            load([data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
            masked_data = parcel_time(logical(tmask_concat), :)';
            ind = 2;
            num_nodes = 333;
        end
        if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
            amt_data = size(masked_data,2);
        end
        %datasample(masked_data,amt_data,2,'Replace', false);
        %disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(matched_data_1,1), size(matched_data_1,2)))
        % ... calculate the correlation matrix...
        
        amt_data = floor(amt_data/2);
        matched_data_sorted = masked_data(atlas_params.sorti,:);
        matched_data_1 = matched_data_sorted(:,1:amt_data);
        matched_data_2 = matched_data_sorted(:,end-amt_data+1:end);
        matrix_1 = paircorr_mod(matched_data_1');
        matrix_1 = single(FisherTransform(matrix_1));% fisher transform r values
        matrix_2 = paircorr_mod(matched_data_2');
        matrix_2 = single(FisherTransform(matrix_2));
        
        for half = 1:2
            count = 1;
                if half == 1
                    matrix = matrix_1;
                elseif half == 2
                    matrix = matrix_2;
                end
        for net = 1:size(atlas_params.networks,ind)
            if net == 1
                % we are ignoring ROI's labeled as 'Unassigned'
                rois = [count:net_size(net)];
                count = rois(end) + 1;
                continue;
            else
                rois = [count:(count-1) + net_size(net)]; %extract the rois belonging to system n
            end            
            tmp_within = matrix(rois,rois); % within-network correlations
            maskmat = ones(size(tmp_within));
            maskmat = logical(triu(maskmat,1));
            within = tmp_within(maskmat);
            if strcmpi(neg_corrs, 'nan')
                within(within<0) = [];
            elseif strcmpi(neg_corrs, 'zero')
                within(within<0) = 0;
            end
            all_within = [all_within; within];
            %all_within(net,ses) = all_within;
            
            tmp_between = matrix(rois(1):rois(end), 1:num_nodes); % all network correlations
            maskmat = ones(size(tmp_between)); % mask out within-network correlations
            maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
            between = tmp_between(maskmat==1); %between-network correlations
            if strcmpi(neg_corrs, 'nan')
                between(between<0) = [];
            elseif strcmpi(neg_corrs, 'zero')
                between(between<0) = 0;
            end
            all_between = [all_between;between];

            count = count + net_size(net);
        end
            %% calculate the segregation index by network by session
            ses_SI(s,half) = (mean(all_within) - mean(all_between))/mean(all_within); 
        end
        end
    end

save([output_dir data_set '_allsubs_seg_index_splitfirstses_negcorrs' neg_corrs '_' atlas '.mat'], 'ses_SI')

