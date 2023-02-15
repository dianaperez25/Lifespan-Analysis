%% Script for Reliability Analysis
% This script will concatenate each subject's timeseries across runs and sessions, 
% then split them into a true half and a separata data pool
% Then create connectivity matrices for the true half and a subset of the data pool 
% to calculate the similarity and test-retest reliability
clear all
% ------------------------------------------------------------------------

%272 frames = 4.99 min, 
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
<<<<<<< Updated upstream
output_dir = '/Users/dianaperez/Desktop/Segregation_Analyses/';

=======
>>>>>>> Stashed changes
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = 5;
%runs = [9,9,11,8,9;8,8,8,9,9];
%pts2sample = 8181; %8181 roughly equivalent to 150 minutes
<<<<<<< Updated upstream
%pts2sample = 5454; %number of frames to sample for true half;roughly equals 100 minutes
pts2sample = 3818;
%sampStep=272; %5 minutes, will add this number of frames each time it subsamples data
sampStep=136;% roughly 2.5 mins -- 272 = 5 minutes, will add this number of frames each time it subsamples data
iterations = 10;
rgb_colors = [1 0 0;%LS02
            0, 1, 0;%LS03
            0, 0, 1;%LS05
            0, 1, 1;%LS08
            1, 0, 1;%LS11
            0.4660 0.6740 0.188;%LS14
            0.9290 0.6940 0.1250;%LS16
            0.4940 0.1840 0.5560];%LS17
        
=======
pts2sample = 5454; %number of frames to sample for true half;roughly equals 100 minutes
sampStep=272; %5 minutes, will add this number of frames each time it subsamples data
perms = 1000;

>>>>>>> Stashed changes
for s = 1:numel(subject)
    
    catData = [];
    catTmask = [];
    
    for i = 1:sessions
        
        %load mat file
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        %get number of runs
        num_runs = size(sess_roi_timeseries,2);
        %randomize runs
        runs_rand = randperm(num_runs);
        for 1:num_runs            
            %apply tmask
            masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));

        %concatenate data
    %        catData = [catData masked_data];
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
    maskmat = ones(300);
    maskmat = logical(triu(maskmat, 1));
    truehalf_corrlin(1,:) = corrmat_truehalf(maskmat);

    rest_of_data = catData(:,ind(pts2sample+1:end));
    times = floor((size(rest_of_data,2))/sampStep);
    times = [2.5:2.5:(times*2.5)];
    %calculate how many samples we have data for
   
    
    %sample data, make corrmats for each sample, make linear matrix with
    %corrmats
    for t = 1:numel(times)
        sampledDatas{t} = datasample(rest_of_data,sampStep*t,2);
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
    allsubs_corrs{s} = corr;
    %variable with amount of time sampled in each sample for each subject
    times_all(s,1:size(times,2)) = times;
    
    clear corr
    clear times
    clear catData
    clear catTmask
    clear corrlins
    clear corrmat
    clear corrmat_truehalf
    clear corrs
    clear masked_data
    clear sampledDatas
    clear tmask
    clear tmask_concat
    clear truehalf
    clear truehalf_corrlin
    clear sess_roi_timeseries
    clear sess_roi_timeseries_concat
    
end

<<<<<<< Updated upstream
corrs_for_mean = [corr(1:3,1:20); corr(5:8,1:20)];
mean = mean(corr,1);%mean(corrs_for_mean);
=======
%corrs_for_mean = [corr(1:3,1:20; corr(5:6,1:20)];
mean = mean(corrs_for_mean);
>>>>>>> Stashed changes
%plot reliability curves
times =[5:5:100];
figure;
plot(times(1:20),corr(1,1:20),'Color',[1, 0, 0],'LineWidth', 3) %LS02
hold on
plot(times(1:20),corr(2,1:20),'Color',[0, 1, 0],'LineWidth', 3) %LS03
hold on
plot(times(1:20),corr(3,1:20),'Color',[0, 0, 1],'LineWidth', 3)%LS05
hold on
plot(times(1:3),corr(4,1:3),'Color',[0, 1, 1],'LineWidth', 3)%LS08
hold on
plot(times(1:20),corr(5,1:20),'Color',[1, 0, 1],'LineWidth', 3)%LS11
hold on
plot(times(1:20),corr(6,1:20),'Color',[0.4660 0.6740 0.1880],'LineWidth', 3)%LS14
hold on
plot(times(1:9),corr(7,1:9),'Color',[0.9290 0.6940 0.1250],'LineWidth', 3)%LS14
hold on
plot(times(1:14),corr(8,1:14),'Color',[0.4940 0.1840 0.5560],'LineWidth', 3)%LS14
hold on
plot(times(1:20),mean(1:20), ':', 'Color', [0,0,0], 'LineWidth',3) %average

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

m = findobj(gca,'Type','line');

hleg1 = legend(m(1:9), 'Mean', 'LS17', 'LS16', 'LS14', 'LS11', 'LS08', 'LS05', 'LS03', 'LS02', 'Location', 'SouthEast');
hleg1.FontSize = 20;
ax = gca;
ax.FontSize = 24;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);

print(gcf,['/Volumes/RESEARCH_HD/Lifespan/CNS_Analyses/ReliabilityLifespanRestDatatruhalf' num2str(pts2sample) '.jpg'],'-dpng','-r300');





