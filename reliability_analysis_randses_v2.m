% %% Test-Retest Reliability of rs-FC in young and older adults 
% This script conducts a test-retest reliability analysis by:
% 1) concatenating all the [subject's data,' ...
% 2) splitting up the data into a "true half" (a large chunk of the [subject's' ...
%    data) and an independent data pool
% 3) selecting subsamples from the independent data pool with increasing amounts
%    of data
% 4) comparing the correlation matrix resulting from the true half to the matrix
%    from each independent subsample of data using correlation
% 5) plotting the reliability curves for each subject + a group mean
% 
% User selects:
% - how much data to include in the "true half"
% - how many iterations of the analysis to conduct per subject
% - how much data does each subsample increase by

clear all

%% PATHS
data_dir = '/Users/dianaperez/Desktop/FC_Parcels_333'; % the timecourses for the different parcels
output_dir = '/Users/dianaperez/Desktop/'; % where to store the results

%% OPTIONS
datasets = {'Lifespan-FSU', 'iNet-FSU'};%'Lifespan-NU', 'iNet-NU', 

% How many points to sample for "true" half
truehalf_datapts = 3808; %3808 -> ~70 min, 5454 -> ~100 min, 8181 -> ~150 min

% How much data to add at each step in minutes (how many data points that
% corresponds to is calculated below)
step_time = 2.5;
step_datapts = round((60*step_time)/1.1);% 136 -? ~2.5 min, 272 -> ~5 min, will add this number of frames each time it subsamples data

% How many iterations to run
iterations = 1000;

for d = 1:numel(datasets)
    dataset = datasets{d};
    t0=tic;
    %% VARIABLES
if strcmpi(dataset, 'Lifespan-NU')
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
sw_subs=[]; sw_iters=[];

for sub = 1:numel(subject)    
    data_struct = cell(1,sessions);
    total_datapts = 0;
    tStart=tic;
    for i = 1:sessions        
        %load mat file
        FC_fname = [data_dir '/sub-' subject{sub} '_rest_ses-' num2str(i) '_parcel_timecourse.mat'];
        if exist(FC_fname)
            parcel_time=[];
            load(FC_fname)
            masked_data = parcel_time(logical(tmask_concat'),:);
            data_struct{i} = masked_data';
            total_datapts = total_datapts + size(masked_data,1);            
        end
    end

    disp(sprintf('Begin Processing for subject %s : Total number of sample points is %d...', subject{sub}, total_datapts))
    %disp(sprintf('data loading and storage finished -- time elapsed: %d', toc(t0)))
    clear masked_data i tmask_concat parcel_time FC_fname
        
    % create a matrix with indices for data points belonging to each chunk
    % of data. Each row corresponds to a chunk of sampStep data points
    num_sets = floor(total_datapts/step_datapts); % calculate how many chunks of data are possible given the total number of data points
    indices = [1:1:num_sets*step_datapts]; % create a matrix of indices for each data point
    indices_for_data = reshape(indices, step_datapts, num_sets)';
    clear indices

    fprintf('Iteration: 0000');
    for p = 1:iterations
        fprintf('\b\b\b\b%4d', p);
        %disp(sprintf('begin iteration %d', p))
        t1=tic;
        % randomize order of sessions
        rng('shuffle');
        ses_ind = randperm(sessions);
    
        % concatenate data for all sessions in random order - note, the order
        % of the runs stays the same
        cat_data = zeros(333,total_datapts);
        first = 1;
        for i = 1:sessions
            data = data_struct{ses_ind(i)};
            last = first + size(data,2) - 1;
            cat_data(:,first:last) = data(:,:);
            first = last + 1;
        end

        clear first last data i
        %disp(sprintf('data concatenation finished: time elapsed %d minutes (%d seconds since start of iteration)', toc(t0)/60, toc(t1)))

        % now randomly pick a starting point to begin sampling
        %t2=tic;
        count = datasample(1:num_sets,1);
        true_half_inds = []; true_half_sets = [];
        
        while length(true_half_inds) < truehalf_datapts
            
            if count > num_sets % if we reach the end of contigous chunks, circle back to the beginning
                count = 1;            
            end
        
            true_half_inds = [true_half_inds; indices_for_data(count,:)']; % all the indices for the data points that will be in the true half of this iteration
            count = count + 1;
            true_half_sets = [true_half_sets; count]; % count the number of sets so we can delete them later
        end
                
        %disp(sprintf('sampling of true half finished: time elapsed %d minutes (%d seconds since last step)', toc(t0)/60, toc(t2)))

        %make corrmats for true half
        %t3=tic;
        truehalf = cat_data(:,true_half_inds);
        corrmat_truehalf = paircorr_mod(truehalf');
        maskmat = ones(333);
        maskmat = logical(triu(maskmat, 1));
        truehalf_corrlin(1,:) = corrmat_truehalf(maskmat);
        
        %disp(sprintf('correlation of true half finished: time elapsed %d minutes (%d seconds since last step)', toc(t0)/60, toc(t3)))

        % now let's chunk the rest of the data, after excluding the true half
        %t4=tic;
        rest_of_data = cat_data;
        rest_of_data(:,true_half_inds) = []; % delete the chunks of data that were used for true half because independent samples
        indices_for_rest_of_data = indices_for_data; % let's copy this matrix, so we can use same strategy for sampling increasing amount of data
        indices_for_rest_of_data(end-(length(true_half_sets)-1):end,:) = []; % but let's delete the chunks that equate to the amount of data that was deleted so we don't get an error
        times = [2.5:2.5:size(indices_for_rest_of_data,1)*2.5]; % calculate the data steps
        clear cat_data count true_half_inds

        %Let's start sampling data, make corrmats for each sample, make linear matrix with
        %corrmats
        inds_this_chunk = [];
        set = datasample(1:size(indices_for_rest_of_data,1),1);

        for t = 1:numel(times) % for each of the time steps that we calculated before...
            %sets = [];               
             % randomly select a starting point from which we'll start a count            
            %while length(sets) < t % while we are still sampling sets
                if set > size(indices_for_rest_of_data,1) %if we run out of sets
                    set = 1; % circle back to beginning
                end
                %sets = [sets; count]; % else, keep adding a contigous set at a time
                %count = count + 1;
            %end
            inds_this_chunk = [inds_this_chunk; indices_for_rest_of_data(set,:)'];
            sampledDatas{t} = rest_of_data(:,inds_this_chunk);
            corrmat = paircorr_mod(sampledDatas{t}');
            corrs{t} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
            corrlins(t,:) = corrmat(maskmat);
            set = set + 1;
        end
        
        %disp(sprintf('sampling of data subsamples finished: time elapsed %d minutes (%d seconds since last step)', toc(t0)/60, toc(t4)))

    %run corrs between true half corrmats and each samples corrmats
    %t5=tic;
    for j = 1:size(corrlins,1)
        tmpcorr = corrcoef(truehalf_corrlin', corrlins(j,:)', 'rows', 'complete');
        corr(p,j) = tmpcorr(2,1);
        clear tmpcorr
    end
    sw_iters(sub,p)=toc(t1);
    end
    fprintf('\n');

    allsubs_corrs{sub} = corr;
    means{sub} = mean(corr,1);
    %disp(sprintf('correlation of true half and subsamples finished: time elapsed %d minutes (%d since last step)', toc(t0)/60, toc(t5)))
    %variable with amount of time sampled in each sample for each subject
    times_all(sub,1:size(times,2)) = times;
    
    clear cat_data
    clear cat_tmask
    clear masked_data
    clear num_sets
    clear indices_for_data
    clear truehalf
    clear true_half_corrlin
    clear rest_of_data
    clear indices_for_rest_of_data
    clear corr
    clear times
    clear sampled_datas
    clear corrlins
    clear this_corrmat
    clear truehalf_corrmat
    clear corrmats
    
    clear sampled_datas
    clear tmask
    clear tmask_concat
    clear truehalf
    clear truehalf_corrlin
    clear sess_roi_timeseries
    clear sess_roi_timeseries_concat
   disp(sprintf('Finished Processing for sub-%s: time elapsed %d minutes', subject{sub}, toc(tStart)/60)) 
   sw_subs=toc(tStart);
end
disp(sprintf('Finished Processing for Dataset %s', dataset))
disp('STATS:')
disp(sprintf('Total time elapsed - %d hours', toc(t0)/360))
disp(sprintf('Total number of subjects - %d', numel(subject)))
disp(sprintf('Average time elapsed per iteration - %d seconds', mean(mean(sw_iters))))
disp(sprintf('Average time elapsed per subject - %d minutes', mean(mean(sw_subs))/60))

% Take mean across subjects of the mean of the correlation values across
% iterations
save([output_dir dataset '_Reliability_truehalf_' num2str(truehalf_datapts) '_corrdata.mat'], 'allsubs_corrs', 'means', 'times_all', '-v7.3')
figure;

if strcmpi(dataset, 'Lifespan-NU')

    for s = 1:numel(subject)
       plot(times_all(s,1:size(means{1,s},2)),means{1,s},'Color',rgb_colors(s,:),'LineWidth', 3)
       hold on
    end
else
    for s = 1:numel(subject)
       plot(times_all(s,1:size(means{1,s},2)),means{1,s}, 'LineWidth', 2)
       hold on
    end
end

for t = 1:size(times_all,2)
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

times_all = step_time:step_time:size(times_all,2)*step_time;
plot(times_all,mean_of_means(1:size(times_all,2)), ':', 'Color', [0,0,0], 'LineWidth',3) %average

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

print(gcf,[output_dir dataset '_Reliability_truehalf_' num2str(truehalf_datapts) '.jpg'],'-dpng','-r300');
close all

end


