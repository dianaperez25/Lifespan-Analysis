%% Script for Similarity analysis, but with forcing the same amount of data per subject

clear all
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
outDir = '/Users/dianaperez/Desktop/';
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14'};
sessions = 5;
%runs = [9,9,11,8,9;8,8,8,9,9];
allSubs_amtData = [];

%% add code to separate by session, get minimum amt of data per session
for s = 1:numel(subject)  
    for i = 1:sessions
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        allSubs_amtData = [allSubs_amtData; size(masked_data,2)];
    end
end

min_data = min(allSubs_amtData);
maskmat = ones(300);
maskmat = logical(triu(maskmat, 1));

count = 1;
for s = 1:numel(subject)  
    catData = [];
    for i = 1:sessions
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        matched_data = masked_data(:,1:min_data);
%         catData = [catData masked_data(1:min_data)];
        disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(matched_data,1), size(matched_data,2)))
        corrmat_matched_data = paircorr_mod(matched_data');
        matcheddata_corrlin(count,:) = single(FisherTransform(corrmat_matched_data(maskmat)));
        count = count + 1;
    end
end

simmat = corr(matcheddata_corrlin');

%% Make Figure --- edit tick marks!!!!
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([5,10,15,20,25,30]+0.5,'k',2);
vline_new([5,10,15,20,25,30]+0.5,'k',2);
set(gca,'XTick',[4,13.5,19,24.5,30.5,34,37,40,45], 'YTick', [4,13.5,19,24.5,30.5,34,37,40,45], 'XTickLabel',...
     {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14'}, 'YTickLabel', {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outDir 'SimilarityMat_rest_matched_data.tiff'],'tiff');
close('all');

    
    
    
   



