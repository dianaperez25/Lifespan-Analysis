%% Script for Reliability of rs-FC Analysis

clear all

%% PATHS
data_dir = '/Volumes/Back_Up/Dissertation/Lifespan/derivatives/postFCproc_CIFTI/FC_Parcels_333';
output_dir = '/Volumes/Back_Up/Dissertation/Reliability_Analyses';

%% OPTIONS
dataset = 'Lifespan'; %'iNet' or 'Lifespan'
% How many points to sample for "true" half
pts2sample = 3808; %3808 -> ~70 min, 5454 -> ~100 min, 8181 -> ~150 min
% How much data to add at each step
samp_step = 136;% 136 -? ~2.5 min, 272 -> ~5 min, will add this number of frames each time it subsamples data
% How many iterations to run
iterations = 1;%000;

%% VARIABLES

if strcmpi(dataset, 'Lifespan')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
    %rgb colors for plotting results
    rgb_colors = [1 0 0;%LS02
                0, 1, 0;%LS03
                0, 0, 1;%LS05
                0, 1, 1;%LS08
                1, 0, 1;%LS11
                0.4660 0.6740 0.188;%LS14
                0.9290 0.6940 0.1250;%LS16
                0.4940 0.1840 0.5560];%LS17
elseif strcmpi(dataset, 'iNet')
    subject = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
    'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
    'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
    'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %
    sessions = 4;
end
        
for sub = 1:numel(subject)    
    
    data_struct = cell(1,sessions);
    num_pts = 0;

    for i = 1:sessions        
        %load mat file
        FC_fname = [data_dir '/sub-' subject{sub} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'];
        if exist(FC_fname)
            load(FC_fname)
            masked_data = parcel_time(logical(tmask_concat'),:);
            data_struct{i} = masked_data';
            num_pts = num_pts + size(masked_data,1);
        end
        %apply tmask
        %masked_data = parcel_time(logical(tmask_concat'),:);
        %concatenate data
        %cat_data = [cat_data masked_data'];
        
        %catTmask = [catTmask tmask_concat'];
    end

    disp(sprintf('Total number of sample points for subject %s is %d...', subject{sub}, num_pts))
    
    indices = [1:1:num_pts]; % create a matrix of indices for each data point
    num_sets = floor((length(indices))/samp_step); % calculate how many chunks of data are possible given the total number of data points
    clear indices
    count = 1; 
    
    % create a matrix with indices for data points belonging to each chunk
    % of data. Each row corresponds to a chunk of sampStep data points
    for set = 1:num_sets
        indices_for_data(set,:) = count:(count+samp_step-1);
        count = count + samp_step;
    end
    
    for p = 1:iterations
        % randomize order of sessions
        rng('shuffle');
        ses_ind = randperm(sessions);
    
        % concatenate data for all sessions in random order - note, the order
        % of the runs stays the same
        cat_data = zeros(1,num_pts);
        first = 1;
        for i = 1:sessions
            data = data_struct{ses_ind(i)};
            last = first + size(data,2) - 1;
            cat_data(first:last) = data;
            first = last + 1;
        end

        % now randomly pick a set of contigous chunks that add up to the
        % amount of data needed for the true half
        count = datasample(1:num_sets,1);
        true_half_inds = [];
        num_true_half_sets = 1;
        
        while length(true_half_inds) < pts2sample
            
            if count > num_sets % if we reach the end of contigous chunks, circle back to the beginning
                count = 1;            
            end
        
            true_half_inds = [true_half_inds; indices_for_data(count,:)']; % all the indices for the data points that will be in the true half of this iteration
            count = count + 1;
            num_true_half_sets = num_true_half_sets + 1; % count the number of sets so we can delete them later
        end
        
        %make corrmats for true half
        truehalf = cat_data(:,true_half_inds);
        corrmat_truehalf = paircorr_mod(truehalf');
        maskmat = ones(333);
        maskmat = logical(triu(maskmat, 1));
        truehalf_corrlin(1,:) = corrmat_truehalf(maskmat);
        
        % now let's chunk the rest of the data, after excluding the true half
        rest_of_data = cat_data;
        rest_of_data(:,true_half_inds) = []; % delete the chunks of data that were used for true half because independent samples
        indices_for_rest_of_data = indices_for_data; % let's copy this matrix, so we can use same strategy for sampling increasing amount of data
        indices_for_rest_of_data((end-(num_true_half_sets-1)):end,:) = []; % but let's delete the chunks that equate to the amount of data that was deleted so we don't get an error
        times = [2.5:2.5:((size(indices_for_rest_of_data,1))*2.5)]; % calculate the data steps

        %Let's start sampling data, make corrmats for each sample, make linear matrix with
        %corrmats
        
        for t = 1:numel(times) % for each of the time steps that we calculated before...
            inds_this_chunk = []; sets = []; %initialize some vars                 
            count = datasample(1:size(indices_for_rest_of_data,1),1); % randomly select a starting point from which we'll start a count            
            while length(sets) < t % while we are still sampling sets
                if count > size(indices_for_rest_of_data,1) %if we run out of sets
                    count = 1; % circle back to beginning
                end
                sets = [sets; count]; % else, keep adding a contigous set at a time
                count = count + 1;
            end
            inds_this_chunk = [inds_this_chunk; indices_for_rest_of_data(sets,:)'];
            sampledDatas{t} = rest_of_data(:,inds_this_chunk);
            corrmat = paircorr_mod(sampledDatas{t}');
            corrs{t} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
            corrlins(t,:) = corrmat(maskmat);
        end

    %run corrs between true half corrmats and each samples corrmats
    for j = 1:size(corrlins,1)
        tmpcorr = corrcoef(truehalf_corrlin', corrlins(j,:)', 'rows', 'complete');
        corr(p,j) = tmpcorr(2,1);
        clear tmpcorr
    end
    end
    allsubs_corrs{sub} = corr;
    means{sub} = mean(corr);
    %variable with amount of time sampled in each sample for each subject
    times_all(sub,1:size(times,2)) = times;
    
    clear cat_data
    clear cat_tmask
    clear masked_data
    clear num_sets
    clear indices_for_data
    clear true_half
    clear true_half_corrlin
    clear rest_of_data
    clear indices_for_rest_of_data
    clear corr
    clear times
    clear sampledDatas
    clear corrlins
    clear corrmat
    clear corrmat_truehalf
    clear corrs
    
    clear sampledDatas
    clear tmask
    clear tmask_concat
    clear truehalf
    clear truehalf_corrlin
    clear sess_roi_timeseries
    clear sess_roi_timeseries_concat
    
end
 
% Take mean across subjects of the mean of the correlation values across
% iterations

figure;
for s = 1:numel(subject)
   % plot(times_all(s,1:size(means{1,s},2)),means{1,s},'Color',rgb_colors(s,:),'LineWidth', 3)
   plot(times_all(s,1:size(means{1,s},2)),means{1,s})
    hold on
end


for t = 1:48
    tmp = [];
    for s = 1:numel(subject)
        if size(means{1,s},2)>=t
            tmp = [tmp;means{1,s}(t)];
        else
            continue;
        end
    end
    mean_of_means(t) = mean(tmp);
end

plot(times_all(1,1:35),mean_of_means(1:35), ':', 'Color', [0,0,0], 'LineWidth',3) %average

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

%m = findobj(gca,'Type','line');

%hleg1 = legend(m(1:9), 'Mean', 'LS17', 'LS16', 'LS14', 'LS11', 'LS08', 'LS05', 'LS03', 'LS02', 'Location', 'SouthEast');
%hleg1.FontSize = 20;
%ax = gca;
%ax.FontSize = 24;
% clear set
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);

print(gcf,[output_dir 'Lifespan_Reliability_truehalf_' num2str(pts2sample) '.jpg'],'-dpng','-r300');





