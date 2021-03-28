%% Script for Reliability Analysis

clear all
%272 frames = 4.99 min, 
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
subject = {'LS03', 'LS05'};
sessions = 5;
%runs = [9,9,11,8,9;8,8,8,9,9];
%pts2sample = 8181; %8181 roughly equivalent to 150 minutes
pts2sample = 5454; %number of frames to sample for true half;roughly equals 100 minutes
sampStep=272; %5 minutes, will add this number of frames each time it subsamples data

for s = 1:numel(subject)
    
    catData = [];
    catTmask = [];
    
    for i = 1:sessions
        
        %load mat file
        load([dataDir 'sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
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
    truehalf = catData(:,1:pts2sample);
    corrmat_truehalf = paircorr_mod(truehalf');
    maskmat = ones(300);
    maskmat = logical(triu(maskmat, 1));
    truehalf_corrlin(1,:) = corrmat_truehalf(maskmat);

    %calculate how many samples we have data for
    times = floor((size(catData,2)-pts2sample)/sampStep);
    times = [5:5:(times*5)];
    
    %sample data, make corrmats for each sample, make linear matrix with
    %corrmats
    for t = 1:numel(times)
        sampledDatas{t} = catData(:, (pts2sample+1):(pts2sample+(sampStep*t)));
        corrmat = paircorr_mod(sampledDatas{t}');
        corrs{t} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
        corrlins(t,:) = corrmat(maskmat);
    end

    % run corr between true half and remaining data points (will be added
    % at the end as another sample, but will be less than 5 mins)
    sampledDatas{t+1} = catData(:, (pts2sample+1):end);
    corrmat = paircorr_mod(sampledDatas{t+1}');
    corrs{t+1} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
    corrlins((t+1),:) = corrmat(maskmat);
    times = [5:5:((numel(times)+1)*5)];
    
    %run corrs between true half corrmats and each samples corrmats
    for j = 1:size(corrlins,1)
        tmpcorr = corrcoef(truehalf_corrlin', corrlins(j,:)', 'rows', 'complete');
        corr(s,j) = tmpcorr(2,1);
        clear tmpcorr
    end
    
    %variable with amount of time sampled in each sample for each subject
    times_all(s,1:size(times,2)) = times;
    
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

%plot reliability curves
times =[5:5:100];
figure;
plot(times,corr(1,:),'Color',[1, 0.5, 0],'LineWidth', 3)
ylim([0 1]);
hold on
plot(times(1:18),corr(2,1:18),'Color',[0, 0, 1],'LineWidth', 3)
hold on

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

m = findobj(gca,'Type','line');

hleg1 = legend(m(1:2), 'LS05', 'LS03', 'Location', 'SouthEast');
hleg1.FontSize = 14;
ax = gca;
ax.FontSize = 17;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);

print(gcf,['/projects/p31161/ReliabilityLifespanRestDatatruhalf' num2str(pts2sample) '.jpg'],'-dpng','-r300');





