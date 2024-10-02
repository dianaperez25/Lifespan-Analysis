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
datasets = {'Lifespan-NU', 'iNet-NU', 'Lifespan-FSU', 'iNet-FSU'};% 
exclude_subs = {'LS46', 'INET108'};

atlas = 'Parcels333';
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 6252; % if this is commented out or set to 0, then the script will calculate it

% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
data_dir = '/Users/dianaperez/Desktop/FC_Parcels_333';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
%% DATA MATCHING
% ------------------------------------------------------------------------
if match_data
    if amt_data == 0 || ~exist('amt_data')
        allSubs_amtData = []; count = 1;
        for d = 1:numel(datasets)             
            [subjects, sessions, N] = get_subjects(datasets{d}, exclude_subs);            
            for sub = 1:N
                num_data = 0;
                for ses = 1:sessions
                    data_file = [data_dir '/sub-' subjects{sub} '_rest_ses-' num2str(ses) '_parcel_timecourse.mat'];
                    if exist(data_file, 'file')
                        load(data_file)
                        masked_data = parcel_time(logical(tmask_concat),:)';
                        num_data = num_data + size(masked_data,2); clear masked_data parcel_time tmask_concat
                    end
                end
                allSubs_amtData(count) = num_data;
                count = count + 1;
            end
        end
        amt_data = min(min(allSubs_amtData));
    end
end

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
% get atlas information
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
for net = 1:13
    net_size(net) = length(atlas_params.mods{1,net});
end
num_parcs = atlas_params.num_rois; 
weights = net_size./(num_parcs-net_size(1));


for d = 1:numel(datasets)
    disp(['dataset: ' datasets{d}])
% get correct subject and session info
    [subjects, sessions, N] = get_subjects(datasets{d}, exclude_subs);


    % begin segregation analysis calculations
    net_SI = {};
    for sub = 1:numel(subjects)
        sub_struct = {};
        all_within = [];
        all_between = [];
        cat_data = [];
        clear subcorrmat mean_matrix

        for ses = 1:sessions
            fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_timecourse.mat', data_dir, subjects{sub}, ses);
            timecourse_struct = load(fname);
            masked_data = timecourse_struct.parcel_time(logical(timecourse_struct.tmask_concat), :)';
            cat_data = [cat_data masked_data]; clear masked_data timecourse_struct
        end

            sorted_data = cat_data(atlas_params.sorti, :);
    
            % ... match amount of data ...
            if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
                amt_data = size(sorted_data,2);
            end

            if size(sorted_data,2) < amt_data
                exclude_subs{end} = subjects{sub};
                N = N - 1;
                disp(sprintf('sub %s excluded due to not enough data', subjects{sub}))
            else
                matched_data = sorted_data(:,1:amt_data);
    
                % ... calculate the correlation matrix ...
                mean_matrix = paircorr_mod(matched_data');
                mean_matrix = single(FisherTransform(mean_matrix));% fisher transform r values
                

        if ndims(mean_matrix)>2
            mean_matrix = squeeze(mean_matrix);
        end
        count = 1;

        for net = 1:size(atlas_params.networks,2)
            if net == 1
                count = net_size(net)+1;
                continue;
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
            
            net_SI.unweighted(sub, net) = (mean(within)-mean(between))/mean(within);
            weight = net_size(net)/300;
            net_SI.weighted(sub,net) = ((mean(within)*weight)-(mean(between)*weight))/(mean(within)*weight);
            count = count + net_size(net);
        end
        
        %% calculate the segregation index by network by session
        sub_SI(sub) = (mean(all_within) - mean(all_between))/mean(all_within);
            end
    end
    save(sprintf('%s/%s_allsubs_seg_index_sub_matched.mat', output_dir, datasets{d}), 'sub_SI', 'net_SI')
end
%% THE END

function [subject, sessions, N] = get_subjects(dataset, exclude_subs)

if strcmpi(dataset, 'Lifespan-NU')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
elseif strcmpi(dataset, 'iNet-NU')
    subject = {'INET002', 'INET003', 'INET005', 'INET006','INET010',...
    'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
    'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET062', 'INET063',...
    'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %'INET061', 'INET001', 
    sessions = 4;
elseif strcmpi(dataset, 'iNet-FSU')
    subject = {'INET074', 'INET075', 'INET077', 'INET078', 'INET083', 'INET084',...
    'INET085', 'INET086', 'INET087', 'INET088',...
    'INET091', 'INET093', 'INET094', 'INET096', 'INET098', 'INET099', 'INET101',...
    'INET103', 'INET104', 'INET105', 'INET106', 'INET107', 'INET108', 'INET109',...
    'INET110', 'INET111', 'INET112', 'INET114', 'INET115', 'INET116', 'INET120',...
    'INET123', 'INET133', 'INET136', 'INET137', 'INET140', 'INET141', 'INET143',...
    'INET156', 'INET158', 'INET159', 'INET160', 'INET165', 'INET168'};
    sessions = 4;
elseif strcmpi(dataset, 'Lifespan-FSU')
    subject = {'LS31', 'LS32', 'LS33', 'LS39', 'LS43', 'LS44', 'LS45', 'LS46',...
    'LS47', 'LS54', 'LS61', 'LS62', 'LS63', 'LS68', 'LS70', 'LS71', 'LS72',...
    'LS76', 'LS77', 'LS79', 'LS85', 'LS89', 'LS94', 'LS108'};
    sessions = 5;
end
subject(:,find(contains(subject, exclude_subs))) = [];
N = numel(subject);

end