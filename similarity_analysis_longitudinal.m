% Look at relationship between correlation matrices across sessions and
% subjects
clear all

%%PATHS
% QUEST
%datadir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
% outDir = '/home/dcr8536/';
% atlas_dir = '/projects/b1081/Atlases/';
% addpath(genpath('/home/dcr8536/Repositories/GrattonLab-General-Repo/FCPROCESS'));
% addpath(genpath('/projects/b1081/Scripts'))
% addpath(genpath('/projects/b1081/Darmouth_MIND2'))

% LOCAL
datadir = '/Volumes/RESEARCH_HD/Lifespan/CNS_Analyses/FC_Parcels_333/';
outDir = '/Volumes/RESEARCH_HD/Lifespan/CNS_Analyses/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
addpath(genpath('/Users/dianaperez/Documents/GitHub/GrattonLab-General-Repo/FCPROCESS'));
%addpath(genpath('/Users/dianaperez/Box/Scripts'))
addpath(genpath('/Users/dianaperez/Documents/Dependencies/cifti-matlab-master'))
%addpath(genpath('/Users/dianaperez/Box/Quest_Backup/Darmouth_MIND2'))

%% VARIABLES
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
%subs = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006', 'INET010', 'INET016',...
%    'INET018', 'INET019', 'INET030'};
%sessions = [13, 4, 4, 4, 4, 4, 4, 4, 4, 4];
sessions = [8, 10, 10, 5, 5, 5, 5, 5];
avg_sessions = 0;
fisher = 0;
atlas_params = atlas_parameters_GrattonLab('Parcels333',atlas_dir);

% count = 1;
% for sub = 1:numel(subs)
%     for ses = 1:sessions(sub)
%         data = load(sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', datadir, subs{sub}, ses));
%         %data = load(sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', datadir, subs{sub}, subs{sub}, ses));
%         corrmat(count,:,:,:) =  data.parcel_corrmat;
%         %corrmat(count,:,:,:) =  data.corrmat;
%         sub_corrmat(ses,:,:,:) = data.parcel_corrmat;
%         count = count + 1;
%     end
%     if ses > 1
%         mean_corrmat(sub,:,:,:) = squeeze(mean(sub_corrmat));
%         clear sub_corrmat
%     else
%         mean_corrmat(sub,:,:,:)= squeeze(sub_corrmat);
%     end
% end


% %% Make similarity matrices
% maskmat = ones(atlas_params.num_rois);
% maskmat = logical(triu(maskmat,1));
% 
% count = 1;
% for sub = 1:numel(subs)
%     for ses = 1:sessions(sub)
%         tmp = corrmat(count,:,:,:);
%         if fisher
%             corrlin(count,:) = single(FisherTransform(tmp(maskmat)));
%         else
%             corrlin(count,:) = single(tmp(maskmat));
%         end
%         count = count + 1;
%     end
% end
% 
% simmat = corr(corrlin');
% 
% %% Make Figure --- edit tick marks!!!!
% % figure('Position',[1 1 1000 800]);
% % imagesc(simmat,[0 1]); colormap('jet');
% % hline_new([8,18,19,29,31,36,37,42,47,52]+0.5,'k',2);
% % hline_new([3,13,24]+0.5,'k',.5);
% % vline_new([8,18,19,29,31,36,37,42, 47,52]+0.5,'k',2);
% % vline_new([3,13,24]+0.5,'k',.5);
% % set(gca,'XTick',[2,13.5,19,24.5,30.5,34,37,40,45, 50, 55], 'YTick', [4,13.5,19,24.5,30.5,34,37,40,45,50, 55], 'XTickLabel',...
% %     {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14', 'LS16', 'LS17'}, 'YTickLabel', {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14', 'LS16', 'LS17'});
% % axis square;
% % colorbar;
% % title('Correlation Matrix Similarity');
% % saveas(gcf,[outDir 'SimilarityMat_rest.tiff'],'tiff');
% % close('all');
% 
% %averaged across sessions
% 
% if avg_sessions  
%     clear corrlin
%     maskmat = ones(atlas_params.num_rois);
%     maskmat = logical(triu(maskmat,1));
%     count = 1;
%     for s = 1:9
%             tmp = mean_corrmat(s,:,:,:);
%             corrlin(count,:) = tmp(maskmat);
%             count = count+1;
%     end
% 
%     simmat = corr(corrlin');
%     figure('Position',[1 1 1000 800]);
%     imagesc(simmat,[0 1]); colormap('jet');
%     hline_new([1:8]+0.5,'k',2);
%     vline_new([1:8]+0.5,'k',2);
%     set(gca,'XTick',[1:9], 'YTick', [1:9], 'XTickLabel',...
%         {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14'}, 'YTickLabel', {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS08', 'LS10', 'LS11', 'LS14'});
%     axis square;
%     colorbar;
%     title('Correlation Matrix Similarity');
%     saveas(gcf,[outDir 'SimilarityMat_averaged.tiff'],'tiff');
%     close('all');
% end


% longitudinal analysis
subs = {'LS02', 'LS03', 'LS05'};
sessions = [8,10,10];
clear corrmat
count = 1;
for sub = 1:numel(subs)
    for ses = 1:sessions(sub)
        data = load(sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', datadir, subs{sub}, ses));
        %data = load(sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', datadir, subs{sub}, subs{sub}, ses));
        corrmat(count,:,:,:) =  data.parcel_corrmat;
        %corrmat(count,:,:,:) =  data.corrmat;
        sub_corrmat(ses,:,:,:) = data.parcel_corrmat;
       count = count + 1;
    end
end
clear corrlin
maskmat = ones(atlas_params.num_rois);
maskmat = logical(triu(maskmat,1));

count = 1;
for sub = 1:numel(subs)
    for ses = 1:sessions(sub)
        tmp = corrmat(count,:,:,:);
        if fisher
            corrlin(count,:) = single(FisherTransform(tmp(maskmat)));
        else
            corrlin(count,:) = single(tmp(maskmat));
        end
        count = count + 1;
    end
end

simmat = corr(corrlin');

figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([0,8,18,28]+0.5,'k',2);
hline_new([0,3,13,23]+0.5,'k',.5);
vline_new([0,8,18,28]+0.5,'k',2);
vline_new([0,3,13,23]+0.5,'k',.5);
set(gca,'XTick',[4,13.5,23.5], 'YTick', [4,13.5,23.5], 'XTickLabel',...
    {'LS02', 'LS03', 'LS05'}, 'YTickLabel', {'LS02', 'LS03', 'LS05'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outDir 'SimilarityMat_longitudinalSubs.tiff'],'tiff');
close('all');



ses_pre = [3,5,5];
ses_post = [5,5,5];
count = 1;
within = [];
between = [];
for s = 1:numel(subject)
    lines = [count:count+(ses_pre(s)+ses_post(s))-1];
    sub_vals = simmat(lines,:);
    maskmat = ones(length(lines),length(lines));
    maskmat = logical(triu(maskmat, 1));
    within_sub = sub_vals(:,lines);
    within = [within; within_sub(maskmat)];
    maskmat = ones(size(sub_vals));
    maskmat(:,lines) = 0;
    between = [between; sub_vals(maskmat==1)];
    count = lines(end)+1;
end

mean(between)
mean(within)

count = 1;
clear corrmat
for sub = 1:numel(subs)
    for ses = 1:ses_pre(sub)
        data = load(sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', datadir, subs{sub}, subs{sub}, ses));
        corrmat(count,:,:,:) =  data.corrmat;
        sub_corrmat(ses,:,:,:) = data.corrmat;
        count = count + 1;
    end
    if ses > 1
        mean_pre(sub,:,:,:) = squeeze(mean(sub_corrmat));
        clear sub_corrmat
    else
        mean_pre(sub,:,:,:)= squeeze(sub_corrmat);
    end
end

count = 1;
for sub = 1:numel(subs)
    for ses = ses_pre(sub)+1:ses_pre(sub)+ses_post(sub)
        data = load(sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', datadir, subs{sub}, subs{sub}, ses));
        corrmat(count,:,:,:) =  data.corrmat;
        sub_corrmat(ses,:,:,:) = data.corrmat;
        count = count + 1;
    end
    if ses > 1
        mean_post(sub,:,:,:) = squeeze(mean(sub_corrmat));
        clear sub_corrmat
    else
        mean_post(sub,:,:,:)= squeeze(sub_corrmat);
    end
end
maskmat = ones(atlas_params.num_rois);
maskmat = logical(triu(maskmat,1));

clear corrlin
count = 1;
for sub = 1:numel(subs)
    tmp = mean_pre(count,:,:,:);
    corrlin_pre(count,:) = single(tmp(maskmat));
    tmp = mean_post(count,:,:,:);
    corrlin_post(count,:) = single(tmp(maskmat));
    count = count + 1;
end

simmat = corr(corrlin_pre', corrlin_post');
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([1:2]+0.5,'k',2);
vline_new([1:2]+0.5,'k',2);
set(gca,'XTick',[1:3], 'YTick', [1:3], 'XTickLabel',...
    {'LS02', 'LS03', 'LS05'}, 'YTickLabel', {'LS02', 'LS03', 'LS05'});
axis square;
colorbar;
title('Correlation Matrix Similarity - 1 year apart');
saveas(gcf,[outDir 'SimilarityMat_longitudinal_avg.tiff'],'tiff');

