%% Script for Similarity analysis, but with forcing the same amount of data per subject

clear all
%dataDir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
dataDir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
outDir = '/Users/dianaperez/Desktop/';
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
sessions = 5;
%runs = [9,9,11,8,9;8,8,8,9,9];
match_data = 1;
if match_data
    allSubs_amtData = [];
end
%% add code to separate by session, get minimum amt of data per session
for s = 1:numel(subject)  
    for i = 1:sessions
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        if match_data
        allSubs_amtData(i,s) = size(masked_data,2);
        end
    end
end

maskmat = ones(300);
maskmat = logical(triu(maskmat, 1));

count = 1;
for s = 1:numel(subject)  
    catData = [];
    for i = 1:sessions
        load([dataDir '/sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        if match_data
            amt_data = min(min(allSubs_amtData));
        else
            amt_data = size(masked_data,2);
        end
        matched_data = datasample(masked_data,amt_data,2,'Replace', false);
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
hline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
vline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
set(gca,'XTick',[3,8,13,18,23,28,33,38], 'YTick', [3,8,13,18,23,28,33,38], 'XTickLabel',...
     {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'}, 'YTickLabel', {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outDir 'SimilarityMat_rest.tiff'],'tiff');
close('all');

% find average within- and between-subject correlations
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

mean(between)
mean(within)
    
    
    
    
    
    
   



