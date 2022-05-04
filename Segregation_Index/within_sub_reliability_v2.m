%% Segregation index with increasing amounts of data
% add 1000 iterations to smooth out data? 
% also figure out how to randomly sample contigous data points

clear all
%272 frames = 4.99 min, 
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
%dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
%subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
subject = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
sessions = 4;
output_dir = '/Users/dianaperez/Desktop/Segregation_Analyses/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas = 'Seitzman300';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
%runs = [9,9,11,8,9;8,8,8,9,9];
pts2sample = 3818; %roughly 70 mins%5454; %number of frames to sample for true half;roughly equals 100 minutes
sampStep=136;% roughly 2.5 mins -- 272 = 5 minutes, will add this number of frames each time it subsamples data
iterations = 1000;
rgb_colors = [1 0 0;%LS02
            0, 1, 0;%LS03
            0, 0, 1;%LS05
            0, 1, 1;%LS08
            1, 0, 1;%LS11
            0.4660 0.6740 0.188;%LS14
            0.9290 0.6940 0.1250;%LS16
            0.4940 0.1840 0.5560];%LS17
        
for s = 1:numel(subject)  
    catData = [];
    catTmask = [];    
    for i = 1:sessions        
        %load mat file
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        %apply tmask
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        %concatenate data
        catData = [catData masked_data];
        catTmask = [catTmask tmask_concat'];
    end    
    %I think this should be 10,816
    disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(catData,1), size(catData,2)))
    indices = [1:1:size(catData,2)]; % create a matrix of indices for each data point
    num_sets = floor((length(indices))/sampStep); % calculate how many chunks of data are possible given the total number of data points
    clear indices
    count = 1; 
    % create a matrix with indices for data points belonging to each chunk
    % of data. Each row corresponds to a chunk of sampStep data points
    for set = 1:num_sets
        indices_for_data(set,:) = count:(count+sampStep-1);
        count = count + sampStep;
    end
    
    for p = 1:iterations
        rng('shuffle');
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
        
        truehalf = catData(:,true_half_inds);
    corrmat_truehalf = paircorr_mod(truehalf');
    [seg_index_truehalf(p,s)] = calculate_seg_index(corrmat_truehalf, atlas_params);

    %calculate how many samples we have data for
    rest_of_data = catData;
        rest_of_data(:,true_half_inds) = []; % delete the chunks of data that were used for true half because independent samples
        indices_for_rest_of_data = indices_for_data; % let's copy this matrix, so we can use same strategy for sampling increasing amount of data
        indices_for_rest_of_data((end-(num_true_half_sets-1)):end,:) = []; % but let's delete the chunks that equate to the amount of data that was deleted so we don't get an error
        times = [2.5:2.5:((size(indices_for_rest_of_data,1))*2.5)]; % calculate the data steps

    
    %sample data, make corrmats for each sample, make linear matrix with
    %corrmats
    
    for t = 1:numel(times)
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
        [sub_seg_index(p,t)] = calculate_seg_index(corrmat, atlas_params);
    end

    
    %run corrs between true half corrmats and each samples corrmats
    for j = 1:size(sub_seg_index,2)
        diffs(p,j) = seg_index_truehalf(p,s) - sub_seg_index(p,j);
        abs_diffs(p,j) = abs(diffs(p,j))/seg_index_truehalf(p,s);
    end
    
    %variable with amount of time sampled in each sample for each subject
    end
    allsubs_seg_ind{s} = sub_seg_index;
    times_all(s,1:size(times,2)) = times;
    allsubs_abs_diffs{s} = abs_diffs;
    
    clear abs_diffs
    clear sub_seg_index
    clear times
    clear catData
    clear catTmask
    clear num_sets
    clear indices_for_data
    clear true_half
    clear true_half_corrlin
    clear rest_of_data
    clear indices_for_rest_of_data
    clear corrmat
    clear corrmat_truehalf
    clear masked_data
    clear sampledDatas
    clear tmask
    clear tmask_concat
    clear truehalf
    clear sess_roi_timeseries
    clear sess_roi_timeseries_concat
    
end

% corrs_for_mean = [corr(1:3,1:20); corr(5:8,1:20)];
% mean = mean(corr,1);%mean(corrs_for_mean);
%plot reliability curves
times =[2.5:2.5:100];
figure;
for s = 1:numel(subject) 
    plot(times_all(s,1:size(allsubs_seg_ind{s},2)),mean(allsubs_abs_diffs{1,s}(:,1:size(allsubs_seg_ind{s},2)),1),'Color',rgb_colors(s,:),'LineWidth', 3)
    hold on
    means_diff{s} = mean(allsubs_abs_diffs{1,s}(:,1:size(allsubs_seg_ind{s},2)),1); 
end
mean_diff_all = [0.451054915	0.31958749	0.250078328	0.205589374	0.173270193	0.149834654	0.130898036	0.116096011	0.103579279	0.09357542	0.085252215	0.077510628	0.071893195	0.067765274	0.063014982	0.0609474	0.057319998	0.055589204	0.053587631	0.051688755	0.049855826	0.049189074	0.047716066	0.047485058	0.046371319	0.046092851	0.045026644	0.049051036	0.048147125	0.04798667	0.047803287	0.048079218	0.047927852	0.047976362	0.048676489	0.048785449	0.048605224	0.048932758	0.049221036	0.049052662	0.049501161	0.049540971	0.049734258	0.04993298	0.050447245	0.050401381	0.0504621	0.050506165	0.048969311	0.042753616	0.042789566	0.043138017	0.043480223	0.026972895];
plot(times_all(1,1:48),mean_diff_all(1:48), ':', 'Color', [0,0,0], 'LineWidth',3)

axis([0 120 0 .6])
ylabel('% Difference');
xlabel('Time (Minutes)');
m = findobj(gca,'Type','line');
hleg1 = legend(m(1:9), 'Mean', 'LS17', 'LS16', 'LS14', 'LS11', 'LS08', 'LS05', 'LS03', 'LS02', 'Location', 'NorthEast');
hleg1.FontSize = 20;
ax = gca;
ax.FontSize = 24;
% 
print(gcf,[output_dir '/ReliabilityLifespanSegIndex_Diff.jpg'],'-dpng','-r300');

figure;
for s = 1:numel(subject)
    plot(times_all(s,1:size(allsubs_seg_ind{s},2)),mean(allsubs_seg_ind{1,s}(:,1:size(allsubs_seg_ind{s},2)),1),'Color',rgb_colors(s,:),'LineWidth', 3)
    hold on
end

means = [0.228121374	0.282854956	0.311825923	0.330394855	0.344075468	0.3540615	0.362085299	0.368580906	0.374353669	0.378886496	0.38303392	0.386997298	0.39018919	0.392779018	0.395671349	0.3921087	0.394652817	0.396349383	0.398223	0.39998822	0.401678269	0.402921549	0.404130596	0.405287927	0.406779144	0.407857626	0.408819803	0.40779906	0.4087539	0.409494202	0.41042813	0.411371703	0.4121878	0.412706817	0.413360333	0.414032785	0.41465492	0.419794232	0.420634614	0.420983464	0.421506702	0.421954376	0.422355276	0.42278862	0.423175282	0.42350195	0.423740618	0.42402081	0.434813233	0.42669529	0.42691267	0.42712188	0.42740375];
plot(times_all(1,1:48),means(1:48), ':', 'Color', [0,0,0], 'LineWidth',3)
axis([0 120 0.1 .5])
ylabel('Segregation Index');
xlabel('Time (Minutes)');
%m = findobj(gca,'Type','line');
%hleg1 = legend(m(1:9), 'Mean', 'LS17', 'LS16', 'LS14', 'LS11', 'LS08', 'LS05', 'LS03', 'LS02', 'Location', 'SouthEast');
%hleg1.FontSize = 20;
ax = gca;
ax.FontSize = 20;

print(gcf,[output_dir '/ReliabilityLifespanSegIndex.jpg'],'-dpng','-r300');

% figure;
% for s = 1:numel(subject)
%     plot(times_all(s,1:size(allsubs_seg_ind{s},2)),smooth(allsubs_seg_ind{1,s}),'Color',rgb_colors(s,:),'LineWidth', 3)
%     hold on
% end

function [seg_index] = calculate_seg_index(matrix, atlas_params)
all_within = [];
all_between = [];
matrix = single(FisherTransform(matrix));% fisher transform r values
for net = 1:14
    net_size(net) = length(atlas_params.mods{1,net});
end    
count = 1;
    for net = 1:size(atlas_params.networks,1)
        if net == 1
            rois = [count:net_size(net)];
        else
            rois = [count:(count-1) + net_size(net)]; %extract the rois belonging to system n
        end
        tmp_within = matrix(rois,rois); % within-network correlations
        maskmat = ones(size(tmp_within));
        maskmat = logical(triu(maskmat,1));
        within = tmp_within(maskmat);
        within(within<0) = [];
        all_within = [all_within; within];
        %all_within(net,ses) = all_within;

        tmp_between = matrix(rois(1):rois(end), 1:300); % all network correlations
        maskmat = ones(size(tmp_between)); % mask out within-network correlations
        maskmat(:,rois(1):rois(end)) = 0; % mask out within-network correlations
        between = tmp_between(maskmat==1); %between-network correlations
        between(between<0) = [];
        all_between = [all_between;between];

        count = count + net_size(net);
    end
        %% calculate the segregation index by network by session
        seg_index = (mean(all_within) - mean(all_between))/mean(all_within);                      
end