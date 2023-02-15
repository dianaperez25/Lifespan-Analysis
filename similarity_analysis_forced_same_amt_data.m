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
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/';
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
output_dir = '/Users/dianaperez/Desktop/';

%% VARIABLES
LS_subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
LS_sessions = 5;
iNet_subject = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
iNet_sessions = 4;

%% OPTIONS
data_type = 'vol'; %'vol' or 'surf' ; are the files for volume ROI's or parcels on the surface?
data_set = 'iNetworks'; %'Lifespan' or 'iNetworks' 
match_data = 0; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 968; % if this is commented out or set to 0, then the script will calculate it

%% DATA MATCHING
if match_data
    if ~exist(amt_data) || amt_data == 0
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
                if strcmpi(data_type, 'vol')            
                    load([data_dir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
                    masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
                elseif strcmpi(data_type, 'surf')
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
% sets the size of the matrix depending on the parcellation/ROI's used
if strcmpi(data_type, 'vol') 
    maskmat = ones(300);
elseif strcmpi(data_type, 'surf')
    maskmat = ones(333);
end
maskmat = logical(triu(maskmat, 1));

% sets number of subjects and sessions depending on the dataset being
% analyzed
if strcmpi(data_set, 'Lifespan')
    subject = LS_subject;
    sessions = LS_sessions;
elseif strcmpi(data_set, 'iNetworks')
    subject = iNet_subject;
    sessions = iNet_sessions;
end

% main loop; starts analysis
count = 1; 
for s = 1:numel(subject)  
    for i = 1:sessions
        % for each session and each subject, load the timeseries data...
        if strcmpi(data_type, 'vol')            
            load([data_dir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
            masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        elseif strcmpi(data_type, 'surf')
            load([data_dir '/sub-' subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
            masked_data = parcel_time(logical(tmask_concat), :)';            
        end
        if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
            amt_data = size(masked_data,2);
        end
        
        % ... then sample the pre-defined amount of data from the timeseries data...
        matched_data = datasample(masked_data,amt_data,2,'Replace', false);
        disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(matched_data,1), size(matched_data,2)))
        % ... calculate the correlation matrix...
        corrmat_matched_data = paircorr_mod(matched_data');
        % ... make it linear and store it in a variable...
        matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));
        % ... then onto the next session.
        count = count + 1;
    end
end

% then, calculate the correlation/similarity across all of those linear matrices
simmat = corr(matcheddata_corrlin');

%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
tick_marks = [0:sessions:(5*numel(subject))]+0.5;
hline_new(tick_marks,'k',2);
vline_new(tick_marks,'k',2);
set(gca,'XTick',tick_marks(1:numel(subject))+(sessions/2), 'YTick', tick_marks(1:numel(subject))+(sessions/2), 'XTickLabel',...
    subject, 'YTickLabel', subject);
axis square;
colorbar;
title('Correlation Matrix Similarity');
if match_data
    saveas(gcf,[output_dir data_set '_' data_type '_similarityMat_matchedData_' num2str(amt_data) '.tiff'],'tiff');
else
    saveas(gcf,[output_dir data_set '_' data_type '_similarityMat_unMatchedData.tiff'],'tiff');
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

disp(['The average similarity between subjects is ' mean(between)])
disp(['The average similarity within subjects is ' mean(within)])
    

%% THE END
    
    
    
    
    
   



