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
if ~exist(output_dir)
    mkdir(output_dir)
end
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
data_type = 'surf'; %'vol' or 'surf' ; are the files for volume ROI's or parcels on the surface?
data_set = 'Lifespan'; %'Lifespan' or 'iNetworks' 
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 988; % if this is commented out or set to 0, then the script will calculate it

%% DATA MATCHING
if match_data
    if amt_data == 0
        allSubs_amtData = [];
        % get minimum amt of data per session
        iNet_subject = [LS_subject iNet_subject];
        allSubs_amtData = [];
        for s = 1:numel(iNet_subject)
            if contains(iNet_subject{s}, 'LS')
                sessions = LS_sessions;
            elseif contains(iNet_subject{s}, 'INET')
                sessions = iNet_sessions;
            else error('Invalid subject ID');
            end

            for i = 1:sessions
                if strcmpi(data_type, 'vol')            
                    load([data_dir '/sub-' iNet_subject{s} '/sub-' iNet_subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
                    masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
                elseif strcmpi(data_type, 'surf')
                    load([data_dir '/sub-' iNet_subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
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
    iNet_subject = LS_subject;
    sessions = LS_sessions;
elseif strcmpi(data_set, 'iNetworks')
    iNet_subject = iNet_subject;
    sessions = iNet_sessions;
end

% main loop; starts analysis
count = 1; 
for s = 1:numel(iNet_subject)  
    for i = 1:sessions
        % for each session and each subject, load the timeseries data...
        if strcmpi(data_type, 'vol')            
            load([data_dir '/sub-' iNet_subject{s} '/sub-' iNet_subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
            masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        elseif strcmpi(data_type, 'surf')
            load([data_dir '/sub-' iNet_subject{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
            masked_data = parcel_time(logical(tmask_concat), :)';            
        end
        if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
            amt_data = size(masked_data,2);
        end
        
        % ... then sample the pre-defined amount of data from the timeseries data...
        matched_data = datasample(masked_data,amt_data,2,'Replace', false);
        disp(sprintf('Total number of sample points for subject %s is %d by %d...', iNet_subject{s}, size(matched_data,1), size(matched_data,2)))
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
tick_marks = [0:sessions:(5*numel(iNet_subject))]+0.5;
hline_new(tick_marks,'k',2);
vline_new(tick_marks,'k',2);
set(gca,'XTick',tick_marks(1:numel(iNet_subject))+(sessions/2), 'YTick', tick_marks(1:numel(iNet_subject))+(sessions/2), 'XTickLabel',...
    iNet_subject, 'YTickLabel', iNet_subject);
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
for s = 1:numel(iNet_subject)
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
disp(['The average similarity between subjects is ' num2str(mean(between))])
disp(['The average similarity within subjects is ' num2str(mean(within))])
save([output_dir data_set '_' data_type '_similarityMat_matchedData.mat'], 'simmat', 'amt_data', 'mean_between', 'mean_within', 'iNet_subject')


%% THE END
    
    
    
    
    
   



