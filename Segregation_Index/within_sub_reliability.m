%% Segregation index with increasing amounts of data
% add 1000 iterations to smooth out data? 
% also figure out how to randomly sample contigous data points

clear all
%272 frames = 4.99 min, 
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = 5;
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
    %make corrmats for true half
    %true-half is 150min=8181 samp points (TR = 1.1; (8181*1.1)/60=149.99)
    for p = 1:iterations
        rng('shuffle');
    ind = randperm(size(catData,2));
    truehalf = catData(:,ind(1:pts2sample));
    corrmat_truehalf = paircorr_mod(truehalf');
    [seg_index_truehalf(p,s)] = calculate_seg_index(corrmat_truehalf, atlas_params);

    %calculate how many samples we have data for
    rest_of_data = catData(:,ind(pts2sample+1:end));
    times = floor((size(rest_of_data,2))/sampStep);
    times = [2.5:2.5:(times*2.5)];
    
    %sample data, make corrmats for each sample, make linear matrix with
    %corrmats
    
    for t = 1:numel(times)
        sampledDatas{t} = datasample(rest_of_data,sampStep*t,2);
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
end

ylabel('% Difference');
xlabel('Time (Minutes)');
m = findobj(gca,'Type','line');
hleg1 = legend(m(1:8), 'LS17', 'LS16', 'LS14', 'LS11', 'LS08', 'LS05', 'LS03', 'LS02', 'Location', 'NorthEast');
hleg1.FontSize = 20;
ax = gca;
ax.FontSize = 24;
% 
print(gcf,['/Volumes/RESEARCH_HD/Lifespan/CNS_Analyses/ReliabilityLifespanSegIndex_Diff.jpg'],'-dpng','-r300');

figure;
for s = 1:numel(subject)
    plot(times_all(s,1:size(allsubs_seg_ind{s},2)),mean(allsubs_seg_ind{1,s}(:,1:size(allsubs_seg_ind{s},2)),1),'Color',rgb_colors(s,:),'LineWidth', 3)
    hold on
end

ylabel('Expectation Value');
xlabel('Time (Minutes)');
m = findobj(gca,'Type','line');
ax = gca;
ax.FontSize = 24;
print(gcf,['/Volumes/RESEARCH_HD/Lifespan/CNS_Analyses/ReliabilityLifespanSegIndex.jpg'],'-dpng','-r300');

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