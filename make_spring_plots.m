
clear all
%% Spring-embedding plots
load('/Users/dianaperez/Desktop/Neurohackademy_Tutorial-master/data/Lifespan_average_correlation_matrix.mat')
corrmat = mean_corrmat;
clear mean_corrmat
% Start with a single correlation matrix per subject and network assignment
% [For the group and for the individual separately?]
%submats = squeeze(mean(corrmat,2));

% and colormat variable needed for function:
% information about prespecified networks
colors = zeros([size(corrmat,1),3]);
atlas_params = atlas_parameters_GrattonLab('Seitzman300','/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/');

for m = 1:length(atlas_params.mods)
    colors(atlas_params.mods{m},:) = repmat(atlas_params.colors(m,:),[length(atlas_params.mods{m}),1]);
end

% Now create across a range of thresholds (1% - 10%, 20%, 30%) for the
% group -- IF TIME
%   1. Discuss consequences of different thresholds
thresholds = 0.05;%[0.01:0.01:0.04]; %[0.01:0.01:0.1,0.2,0.3]; %these take a while, esp for higher thresholds. Only do a few.
outdir = '/Users/dianaperez/Desktop/';
for t = 1:length(thresholds)
    %tic;
    make_spring_fig_MIND(corrmat,thresholds(t),colors,atlas_params.networks,atlas_params.colors);
    saveas(gcf,sprintf('%sSpring_group_t%.02f.tiff',outdir,thresholds(t)),'tiff');
    %toc
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
