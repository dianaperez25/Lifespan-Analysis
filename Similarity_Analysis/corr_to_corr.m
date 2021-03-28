% Look at relationship between correlation matrices across sessions and
% subjects

%datadir = '/projects/b1081/Lifespan/derivatives/preproc_FCProc/corrmats_Seitzman300/';
% outDir = '/home/dcr8536/';
% atlas_dir = '/projects/b1081/Atlases/';
% addpath(genpath('/home/dcr8536/Repositories/GrattonLab-General-Repo/FCPROCESS'));
% addpath(genpath('/projects/b1081/Scripts'))
% addpath(genpath('/projects/b1081/Darmouth_MIND2'))
datadir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
outDir = '/Users/dianaperez/Desktop/Research/Lifespan/';
atlas_dir = '/Users/dianaperez/Box/Quest_Backup/Atlases/';
addpath(genpath('/Users/dianaperez/Documents/GitHub/GrattonLab-General-Repo/FCPROCESS'));
addpath(genpath('/Users/dianaperez/Box/Scripts'))
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
addpath(genpath('/Users/dianaperez/Box/Quest_Backup/Darmouth_MIND2'))
subs = {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS10'};
atlas_params = atlas_parameters_GrattonLab('Seitzman300',atlas_dir);


LS02_1 = load([datadir 'sub-LS02/sub-LS02_sess-1_task-rest_corrmat_Seitzman300.mat']);
LS02_2 = load([datadir 'sub-LS02/sub-LS02_sess-2_task-rest_corrmat_Seitzman300.mat']);
LS02_3 = load([datadir 'sub-LS02/sub-LS02_sess-3_task-rest_corrmat_Seitzman300.mat']);
LS03_1 = load([datadir 'sub-LS03/sub-LS03_sess-1_task-rest_corrmat_Seitzman300.mat']);
LS03_2 = load([datadir 'sub-LS03/sub-LS03_sess-2_task-rest_corrmat_Seitzman300.mat']);
LS03_3 = load([datadir 'sub-LS03/sub-LS03_sess-3_task-rest_corrmat_Seitzman300.mat']);
LS03_4 = load([datadir 'sub-LS03/sub-LS03_sess-4_task-rest_corrmat_Seitzman300.mat']);
LS03_5 = load([datadir 'sub-LS03/sub-LS03_sess-5_task-rest_corrmat_Seitzman300.mat']);
LS04_1 = load([datadir 'sub-LS04/sub-LS04_sess-1_task-rest_corrmat_Seitzman300.mat']);
LS05_1 = load([datadir 'sub-LS05/sub-LS05_sess-1_task-rest_corrmat_Seitzman300.mat']);
LS05_2 = load([datadir 'sub-LS05/sub-LS05_sess-2_task-rest_corrmat_Seitzman300.mat']);
LS05_3 = load([datadir 'sub-LS05/sub-LS05_sess-3_task-rest_corrmat_Seitzman300.mat']);
LS05_4 = load([datadir 'sub-LS05/sub-LS05_sess-4_task-rest_corrmat_Seitzman300.mat']);
LS05_5 = load([datadir 'sub-LS05/sub-LS05_sess-5_task-rest_corrmat_Seitzman300.mat']);
LS07_1 = load([datadir 'sub-LS07/sub-LS07_sess-1_task-rest_corrmat_Seitzman300.mat']);
LS07_2 = load([datadir 'sub-LS07/sub-LS07_sess-2_task-rest_corrmat_Seitzman300.mat']);
LS10_1 = load([datadir 'sub-LS10/sub-LS10_sess-1_task-rest_corrmat_Seitzman300.mat']);

corrmat_2(1,:,:,:) = LS02_1.corrmat;
corrmat_2(2,:,:,:) = LS02_2.corrmat;
corrmat_2(3,:,:,:) = LS02_3.corrmat;

mean_LS02 = squeeze(mean(corrmat_2,1));

corrmat_3(1,:,:,:) = LS03_1.corrmat;
corrmat_3(2,:,:,:) = LS03_2.corrmat;
corrmat_3(3,:,:,:) = LS03_3.corrmat;
corrmat_3(4,:,:,:) = LS03_4.corrmat;
corrmat_3(5,:,:,:) = LS03_5.corrmat;

mean_LS03 = squeeze(mean(corrmat_3,1));

mean_LS04 = LS04_1.corrmat;

corrmat_5(1,:,:,:) = LS05_1.corrmat;
corrmat_5(2,:,:,:) = LS05_2.corrmat;
corrmat_5(3,:,:,:) = LS05_3.corrmat;
corrmat_5(4,:,:,:) = LS05_4.corrmat;
corrmat_5(5,:,:,:) = LS05_5.corrmat;

mean_LS05 = squeeze(mean(corrmat_5,1));

corrmat_7(1,:,:,:) = LS07_1.corrmat;
corrmat_7(2,:,:,:) = LS07_2.corrmat;

mean_LS07 = squeeze(mean(corrmat_7,1));

mean_LS10 = LS10_1.corrmat;

group_corrmat(1,:,:,:) = mean_LS02;
group_corrmat(2,:,:,:) = mean_LS03;
group_corrmat(3,:,:,:) = mean_LS04;
group_corrmat(4,:,:,:) = mean_LS05;
group_corrmat(5,:,:,:) = mean_LS07;
group_corrmat(6,:,:,:) = mean_LS10;

% for s = 1:4 %different subjects, avg over sessions
%     figure_corrmat_GrattonLab(squeeze(group_corrmat(s,:,:,:)),atlas_params,-0.4,1);
%     title(['AllSessAvg, Subject ' subs{s}]);
%     colormap('jet');
%     saveas(gcf, [outDir 'Corrmat_' subs{s} '_sessmean.tiff'],'tiff');
% end
% close('all');

% Make similarity matrices
maskmat = ones(atlas_params.num_rois);
maskmat = logical(triu(maskmat,1));
%not averaged across sessions
tmp = corrmat_2(1,:,:,:);
corrlin(1,:) = tmp(maskmat);
tmp = corrmat_2(2,:,:,:);
corrlin(2,:) = tmp(maskmat);
tmp = corrmat_2(3,:,:,:);
corrlin(3,:) = tmp(maskmat);
tmp = corrmat_3(1,:,:,:);
corrlin(4,:) = tmp(maskmat);
tmp = corrmat_3(2,:,:,:);
corrlin(5,:) = tmp(maskmat);
tmp = corrmat_3(3,:,:,:);
corrlin(6,:) = tmp(maskmat);
tmp = corrmat_3(4,:,:,:);
corrlin(7,:) = tmp(maskmat);
tmp = corrmat_3(5,:,:,:);
corrlin(8,:) = tmp(maskmat);
% tmp = mean_LS04;
% corrlin(9,:) = tmp(maskmat);
tmp = corrmat_5(1,:,:,:);
corrlin(9,:) = tmp(maskmat);
tmp = corrmat_5(2,:,:,:);
corrlin(10,:) = tmp(maskmat);
tmp = corrmat_5(3,:,:,:);
corrlin(11,:) = tmp(maskmat);
tmp = corrmat_5(4,:,:,:);
corrlin(12,:) = tmp(maskmat);
tmp = corrmat_5(5,:,:,:);
corrlin(13,:) = tmp(maskmat);
tmp = corrmat_7(1,:,:,:);
corrlin(14,:) = tmp(maskmat);
tmp = corrmat_7(2,:,:,:);
corrlin(15,:) = tmp(maskmat);
% tmp = mean_LS10;
% corrlin(17,:) = tmp(maskmat);

simmat = corr(corrlin');
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([3,8,13,15]+0.5,'k',2);
vline_new([3,8,13,15]+0.5,'k',2);
set(gca,'XTick',[2,6,11,14.5], 'YTick', [2,6,11,14.5], 'XTickLabel',...
    {'LS02', 'LS03', 'LS05', 'LS07'}, 'YTickLabel', {'LS02', 'LS03', 'LS05', 'LS07'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outDir 'SimilarityMat_rest.tiff'],'tiff');
close('all');

%averaged across sessions
maskmat = ones(atlas_params.num_rois);
maskmat = logical(triu(maskmat,1));
count = 1;
for s = 1:6
        tmp = mean_corrmat(s,:,:,:);
        corrlin(count,:) = tmp(maskmat);
        count = count+1;
end
% end
simmat = corr(corrlin');
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([3,8,9,14,16]+0.5,'k',2);
vline_new([3,8,9,14,16]+0.5,'k',2);
set(gca,'XTick',[2,6,9,12,15.5,17], 'YTick', [2,6,9,12,15.5,17], 'XTickLabel',...
    {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS10'}, 'YTickLabel', {'LS02', 'LS03', 'LS04', 'LS05', 'LS07', 'LS10'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outDir 'SimilarityMat_averaged.tiff'],'tiff');
close('all');

% Talk about doing these measures with individual ROIs

% for more information, relevant references:
%   Power et al. (2011). Functional Network Organization of the Human Brain. Neuron, 72, 4, 665-678
%   Infomap: http://www.mapequation.org/code.html
%       Rosvall & Bergstrom (2008). Maps of invofmation flow reveal
%       community structure in complex networks. PNAS, 105, 1118


%% Spring-embedding plots

% Start with a single correlation matrix per subject and network assignment
% [For the group and for the individual separately?]
submats = squeeze(mean(corrmat_5,2));

% and colormat variable needed for function:
% information about prespecified networks
colors = zeros([size(corrmat_5,1),3]);
for m = 1:length(Parcel_params.mods)
    colors(Parcel_params.mods{m},:) = repmat(Parcel_params.colors(m,:),[length(Parcel_params.mods{m}),1]);
end

% Create spring embedding plots at one threshold for the group
%   1. Discuss different types of thresholding; discuss weighted vs. not
%   2. Discuss issues of network size and density for feasibility
make_spring_fig_MIND(groupmat,0.02,colors,Parcel_params.networks,Parcel_params.colors);

% Now create across a range of thresholds (1% - 10%, 20%, 30%) for the
% group -- IF TIME
%   1. Discuss consequences of different thresholds
thresholds = [0.01:0.01:0.04]; %[0.01:0.01:0.1,0.2,0.3]; %these take a while, esp for higher thresholds. Only do a few.
for t = 1:length(thresholds)
    %tic;
    make_spring_fig_MIND(groupmat,thresholds(t),colors,Parcel_params.networks,Parcel_params.colors);
    saveas(gcf,sprintf('%sSpring_group_t%.02f.tiff',outdir,thresholds(t)),'tiff');
    %toc
end
close all;

% Pick your favorite threshold, and do the same for all subjects -- IF TIME
%   1. Discuss differences
%   2. Discuss consequences of group vs. individual ROIs and network
%   assignments
for s = 1:10
    make_spring_fig_MIND(squeeze(submats(s,:,:)),0.02,colors,Parcel_params.networks,Parcel_params.colors);
    saveas(gcf,sprintf('%sSpring_sub%02d_thresh0.02.tiff',outdir,s),'tiff');
end
close all;

% for more information, relevant references:
%   Bullmore & Sporns (2009). Complex brain networks: graph theoretical analysis of structural and functional systems. Nature Reviews Neuroscience, 10 (3), 186-198
%   Sporns (2010). Networks of the Brain. MIT Press.
%   Spring Embedding Methods (Kamada-Kawai): https://arxiv.org/pdf/1201.3011.pdf
    


%% Hub measures [If time]

thresholds = [0.01:0.01:0.10];

% Start with a set of thresholded correlation matrices and network assignments for each threshold
%   [In the interest of time/ease, I am pre-computing network assignments
%   and providing them here]
infomapcomm = load([datadir 'Allsubavg_333parcels_infomapassn.mat']);
% code for plotting infomap output:
colors = distinguishable_colors(max(unique(infomapcomm.clrs)));
colors(1,:) = [1 1 1];
figure;
imagesc(infomapcomm.clrs(Parcel_params.sorti,:),[1 max(unique(infomapcomm.clrs))]);
hline_new(Parcel_params.transitions,'k',2);
set(gca,'YTick',Parcel_params.centers,'YTickLabels',Parcel_params.networks);
set(gca,'XTick',1:10,'XTickLabels',thresholds);
title('Assignments across thresholds');
colormap(colors);

% Compute hub measures - degree, PC, and WD - at one threshold in one
% subject
%   1. Look at formula, practice computing measures at one thershold in group
%       degree = total number of edges to a node
%       PC = 1 - for all mods (the # of edges of node to mod m/total # of edges of node)
%       WD = z-score normalized measure of degree within an network
%       See formula from Guimera & Amaral (2005). Functional Cartography of
%       Complex Metabolic Networks. Nature, 433, 895-900
%       http://www.nature.com/nature/journal/v433/n7028/full/nature03288.html?foxtrotcallback=true
%       Look at article end for formulas
%   2. Discuss thresholding, distance exclusion
[pc, wd, degree] = module_metrics_Dartmouth(groupmat,infomapcomm.clrs(:,1),0.02,Parcel_params.dist_thresh,Parcel_params.dmat);

% Now do this across thresholds
%   0. Plot image of hub measures across thresholds with network boundaries
%   1. Discuss how they change across thresholds
%   2. Discuss how they change across networks
%   3. Discuss potential issues with each measure
for t = 1:length(thresholds)
    [pc(:,t), wd(:,t), degree(:,t)] = module_metrics_Dartmouth(groupmat,infomapcomm.clrs(:,t),thresholds(t),Parcel_params.dist_thresh,Parcel_params.dmat);
end
figure_hubs(Parcel_params,thresholds,degree,wd,pc);
saveas(gcf,sprintf('%sHubmeasures_group.tiff',outdir),'tiff');

% [If time, do this across other subjects]
for s = 1:10
    infomapcomm = load(sprintf('data/MSC%02d_333parcels_infomapassn.mat',s));
    for t = 1:length(thresholds)
        [pc_sub(:,t,s), wd_sub(:,t,s), degree_sub(:,t,s)] = module_metrics_Dartmouth(squeeze(submats(s,:,:)),infomapcomm.clrs(:,t),thresholds(t),Parcel_params.dist_thresh,Parcel_params.dmat);
    end
    figure_hubs(Parcel_params,thresholds,degree_sub(:,:,s),wd_sub(:,:,s),pc_sub(:,:,s));
    saveas(gcf,sprintf('%sHubmeasures_MSC%02d.tiff',outdir,s),'tiff');
    clear infomapcomm;
end
close('all');

% [If time] Make a spring embedded plot, colored by hub measures rather than networks
%   1. Start with group and favorite threshold. Do other versions if time.
t = 2;
hub_colors = hub_colormap(pc(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_PC_t%.02f.tiff',outdir,thresholds(t)),'tiff');

hub_colors = hub_colormap(wd(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_WD_t%.02f.tiff',outdir,thresholds(t)),'tiff');

hub_colors = hub_colormap(degree(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_degree_t%.02f.tiff',outdir,thresholds(t)),'tiff');


% And for work on hubs and their importance in brain function, in addition to the references above, see:
%   Gratton, C., et al., (2012). Focal brain lesions to critical locations cause widespread disruption of the modular organization of the brain. Journal of Cognitive Neuroscience, 24 (6), 1275-1285
%   Power, J.D. et al. (2013). Evidence for hubs in human functional brain networks. Neuron, 79 (4), 798-813
%   Warren, D.E., et al. (2014). Network measures predict neuropsychological outcome after brain injury. PNAS, 111 (39), 14247-14252
%
% For a general intro on using graph theory in brain networks see: 
%   Bullmore & Sporns (2009). Complex brain networks: graph theoretical
%       analysis of structural and functional systems. Nature Reviews
%       Neuroscience 10(3), 186.
%   Sporns (2010). Networks of the brain
% For extensions of many graph metrics to weighted networks see also:
%   Rubinov & Sporns, 2011. Weight-conserving characterization of
%       complex functional brain networks. Neuroimage, 56(4) 2068-2079.
%   Rubinov & Sporns, 2010. Complex network measures of brain
%       connectivity: uses and interpretations, 52(3) 1059-1069
%
% The following packages contain tools for graph theoretical analyses:
%   Brain Connectivity Toolbox (Sporns, Matlab/Python/C++): https://sites.google.com/site/bctnet/
%   NetworkX (Python): https://networkx.github.io/ 
%       see also brainx extension: https://github.com/nipy/brainx
