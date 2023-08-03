%% Script for Similarity analysis, but with option to force the same amount of data per subject 
% This script calculates the correlation (aka the similarity) across
% sessions and subjects. 
% Input: .mat files containing timecourses for eithe the volume
% ROI's (Seitzman300) or the surface parcels (Gordon333)
% Output: a similarity matrix figure that will be saved, and the average
% between- and within-subject correlations (will not be saved)
% Written to work with both Lifespan and iNetworks subjects assumming 5 and
% 4 sessions for each dataset, respectively
% This script will force the same amount of data across the two datasets

clear all
%% ------------------------------------------------------------------------------------------------
%% PATHS
%data_dir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/iNetworks/Nifti/FC_Parcels_333/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
output_dir = '/Users/diana/Desktop/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';

%% VARIABLES
LS_subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
LS_sessions = 5;
iNet_subject = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %
iNet_sessions = 4;

%% OPTIONS
data_set = 'Lifespan'; %'Lifespan' or 'iNetworks' 
atlas = 'Parcels333'; %Parcels333 for the Gordon surface parcellations or Seitzman300 for the volumetric ROI's
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 988; % if this is commented out or set to 0, then the script will calculate it

%% SEPARATION BY SYSTEMS
% here we specify the indices for the systems that we want to designate to
% each category
SM_systems = [3, 9, 10, 11]; %3: visual, 8: motor hand, 9: motor mouth, 11: auditory
control_systems = [8, 4, 5, 6, 7]; % 8: CON, 4: FPN, 5: DAN, 6: VAN, 7: Salience
control_related = [8,4,5]; %CON, FPN, DAN
memory_default = [6, 7, 2, 12, 13]; %VAN, Salience, DMN, PERN, RetroSpl

% This structure contains the system categories that will be analyzed
% (I made it this way so that we can look at two or three categories without
% changing the script too much)
system_divisions = {SM_systems, control_systems, control_related, memory_default};
output_str = {'sensorimotor', 'control', 'control-related', 'memory-default'}; % output strings for each of the categories being analyzed

% load atlas parameters
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);

%% DATA MATCHING
if match_data
    if amt_data == 0
        allSubs_amtData = [];
        % get minimum amt of data per session
        subject = [LS_subject iNet_subject];
        allSubs_amtData = [];
        for s = 1:numel(subject)
            if contains(subject{s}, 'LS')
                sessions = LS_sessions;
            elseif contains(subject{s}, 'INET')
                sessions = iNet_sessions;
            else error('Invalid subject ID');
            end

            for i = 1:sessions
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

%% SIMILARITY CALCULATIONS

% sets number of subjects and sessions depending on the dataset being
% analyzed
if strcmpi(data_set, 'Lifespan')
    subject = LS_subject;
    sessions = LS_sessions;
elseif strcmpi(data_set, 'iNetworks')
    subject = iNet_subject;
    sessions = iNet_sessions;
end

for sys = 1:numel(system_divisions)

    % get indices for parcels belonging to the systems of interest
    inds = [];
    for n = 1:numel(system_divisions{sys})
        inds = [inds; atlas_params.mods{system_divisions{sys}(n)}];
    end
    
    % main loop; starts analysis
    count = 1; 
    for s = 1:numel(subject)  
        for i = 1:sessions
            % for each session and each subject, load the timeseries data...
            if strcmpi(atlas, 'Seitzman300')            
                load([data_dir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
                masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
            elseif strcmpi(atlas, 'Parcels333')
                load([data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
                masked_data = parcel_time(logical(tmask_concat), :)';            
            end
            if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
                amt_data = size(masked_data,2);
            end

            % ... then sample the pre-defined amount of data from the timeseries data...
            matched_data = datasample(masked_data,amt_data,2,'Replace', false);
            disp(sprintf('Total number of sample points for subject %s session %d is %d by %d...', subject{s}, i, size(matched_data,1), size(matched_data,2)))
            % ... calculate the correlation matrix...
            
            systems_of_interest = matched_data(inds, :);
            corrmat_matched_data = paircorr_mod(systems_of_interest');
            % ... make it linear and store it in a variable...
            
            maskmat = ones(size(corrmat_matched_data,1));
            maskmat = logical(triu(maskmat, 1));
            matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));
            % ... then onto the next session.
            count = count + 1;
        end
    end

% then, calculate the correlation/similarity across all of those linear matrices
simmat = corr(matcheddata_corrlin');
clear matcheddata_corrlin

%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
tick_marks = [0:sessions:(5*numel(subject))]+0.5;
hline_new(tick_marks,'k',1);
vline_new(tick_marks,'k',1);
set(gca,'XTick',tick_marks(1:numel(subject))+(sessions/2), 'YTick', tick_marks(1:numel(subject))+(sessions/2), 'XTickLabel',...
    subject, 'YTickLabel', subject);
axis square;
colorbar;
title(['Correlation Matrix Similarity - ' output_str{sys}]);
if match_data
    saveas(gcf,[output_dir data_set '_' atlas '_' output_str{sys} '_similarityMat_matchedData_' num2str(amt_data) '.tiff'],'tiff');
else
    saveas(gcf,[output_dir data_set '_' atlas '_' output_str{sys} '_similarityMat_unMatchedData.tiff'],'tiff');
end
close('all');

%% CALCULATE average within- and between-subject correlations
count = 1;
within = [];
between = [];
for s = 1:numel(subject)
    lines = [count:(count+sessions-1)];
    sub_vals = simmat(lines,:);
    maskmat = ones(sessions,sessions);
    maskmat = logical(triu(maskmat, 1));
    within_sub = sub_vals(:,lines);
    within = [within; within_sub(maskmat)];
    maskmat = ones(size(sub_vals));
    maskmat(:,lines) = 0;
    between = [between; sub_vals(maskmat==1)];
    count = count+sessions;
end

mean_between = mean(between);
mean_within = mean(within);
disp(['The average similarity between subjects for ' output_str{sys} ' is ' num2str(mean(between))])
disp(['The average similarity within subjects for ' output_str{sys} ' is ' num2str(mean(within))])
systems = system_divisions{sys};
save([output_dir data_set '_' atlas '_' output_str{sys} '_similarityMat_MatchedData.mat'], 'simmat', 'amt_data', 'mean_between', 'mean_within', 'LS_subject', 'systems')
end
%% THE END
    
    
    
    
    
   



