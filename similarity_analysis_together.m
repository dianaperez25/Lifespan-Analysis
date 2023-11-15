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
data_dir = '/Volumes/Back_Up/Dissertation/parcellations/FC_Parcels_333/';
%data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

%% OPTIONS
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 968;%1600; % if this is commented out or set to 0, then the script will calculate it
separate_class = 1;
atlas = 'Parcels333';
%% VARIABLES
subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17',...
'INET001', 'INET002', 'INET003', 'INET005', 'INET050','INET010',...
'INET018','INET056','INET053', 'INET055', 'INET057', 'INET058',...
'INET019', 'INET026', 'INET030',  'INET032', 'INET033','INET006',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET049', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET042', 'INET051', 'INET052',  'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071','INET072', 'INET073'};
%sessions = [5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4];%,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4];

%% SEPARATION BY SYSTEMS
% here we specify the indices for the systems that we want to designate to
% each category
all_systems = 1:13;
SM_systems = [3, 9, 10, 11]; %3: visual, 8: motor hand, 9: motor mouth, 11: auditory
assoc_systems = [2, 8, 4, 5, 6, 7]; % 2: DMN, 8: CON, 4: FPN, 5: DAN, 6: VAN, 7: Salience
control_related = [8,4,5]; %CON, FPN, DAN
memory_default = [6, 7, 2, 12, 13]; %VAN, Salience, DMN, PERN, RetroSpl

% This structure contains the system categories that will be analyzed
% (I made it this way so that we can look at two or three categories without
% changing the script too much)
system_divisions = {all_systems SM_systems assoc_systems control_related memory_default};
output_str = {'all_nets', 'sensorimotor', 'control', 'control_related', 'memory-default'}; % output strings for each of the categories being analyzed


%% DATA MATCHING
if match_data
    if amt_data == 0
        allSubs_amtData = [];
        % get minimum amt of data per session
        allSubs_amtData = [];
        for s = 1:numel(subjects)
            
            for i = 1:sessions(s)
                    load([data_dir '/sub-' subjects{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
                    masked_data = parcel_time(logical(tmask_concat),:)';
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
maskmat = ones(333);
maskmat = logical(triu(maskmat, 1));
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);
for sys = 1:numel(system_divisions)

    % get indices for parcels belonging to the systems of interest
    inds = [];
    for n = 1:numel(system_divisions{sys})
        inds = [inds; atlas_params.mods{system_divisions{sys}(n)}];
    end
    
    % main loop; starts analysis
    count = 1; 
    new_ses_info = [];
    sub_count = 0;
    for s = 1:numel(subjects)  
        ses_count = 0;
        if contains(subjects{s}, 'INET')
            sessions = 4;
        elseif contains(subjects{s}, 'LS')
            sessions = 5;
        end
        for i = 1:sessions
            % for each session and each subject, load the timeseries data...
    
            load([data_dir '/sub-' subjects{s} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'])
            masked_data = parcel_time(logical(tmask_concat), :)';            
            if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
                amt_data = size(masked_data,2);
            end
            
            
            % ... then sample the pre-defined amount of data from the timeseries data...
            if size(masked_data,2)>= amt_data
                ses_count = ses_count + 1;
                matched_data = masked_data(:,1:amt_data);
                disp(sprintf('Total number of sample points for subject %s is %d by %d...', subjects{s}, size(matched_data,1), size(matched_data,2)))
                % ... calculate the correlation matrix...
                systems_of_interest = matched_data(inds, :);
                corrmat_matched_data = paircorr_mod(systems_of_interest');
                %corrmat_matched_data = paircorr_mod(matched_data');
                % ... make it linear and store it in a variable...
                 maskmat = ones(size(corrmat_matched_data,1));
                 maskmat = logical(triu(maskmat, 1));
                matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));
                % ... then onto the next session.
                count = count + 1;
            end
        end
        if ses_count > 0 
            new_ses_info = [new_ses_info, ses_count];
            sub_count = sub_count + 1;
        end
    end
% then, calculate the correlation/similarity across all of those linear matrices
simmat = corr(matcheddata_corrlin');
clear matcheddata_corrlin

%% MAKE FIGURE
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
tick_marks = [0:5:40]+0.5;
tick_marks = [tick_marks [44:4:232]+0.5];
hline_new(tick_marks,'k',.5);
vline_new(tick_marks,'k',.5);
hline_new(40.5,'k',1.5);
vline_new(40.5,'k',1.5);

colorbar;
saveas(gcf,[output_dir 'all_subs_simmat_' output_str{sys} '_matcheddata_' num2str(amt_data) '.tiff'],'tiff');

close('all');

%% CALCULATE average within- and between-subject correlations
count = 1;
within = [];
between = [];
sub_averages = [];
for s = 1:sub_count
    ses = new_ses_info(s);
    lines = [count:(count+ses-1)];
    sub_vals = simmat(lines,:);
    maskmat = ones(ses,ses);
    maskmat = logical(triu(maskmat, 1));
    within_sub = sub_vals(:,lines);
    sub_averages(s,1) = mean(within_sub(maskmat==1));
    within = [within; within_sub(maskmat)];
    maskmat = ones(size(sub_vals));
    maskmat(:,lines) = 0;
    maskmat(:, 41:end) = 2;
    sub_averages(s,2) = mean(sub_vals(maskmat>0));
    if s<=7
        sub_averages(s,3) = mean(sub_vals(maskmat==1));
        sub_averages(s,4) = mean(sub_vals(maskmat==2));
    elseif s>7
        sub_averages(s,3) = mean(sub_vals(maskmat==2));
        sub_averages(s,4) = mean(sub_vals(maskmat==1));
    end
    between = [between; sub_vals(maskmat>0)];
    count = count+ses;
end
save([output_dir 'all_subs_similarity_' output_str{sys} '_matcheddata_' num2str(amt_data) '.mat'],'simmat', 'sub_averages');

disp(['The average similarity between subjects for ' output_str{sys} ' is ' num2str(mean(between))])
disp(['The average similarity within subjects for ' output_str{sys} ' is ' num2str(mean(within))])

end
%% THE END
    
    
    
    
    
   



