%% New Network Segregation Script
% This script calculates a measure of system segregation (how much nodes in
% a system communicate with each other vs with nodes in other systems).
% Input: 300x300 functional connectivity correlation matrix.
% Output: mean within-system correlation, mean between-system correlation,
% a segregation index for each session.
% Based on the analysis describe in: 
% Chan, Micaela Y., et al. "Decreased segregation of brain systems across 
% the healthy adult lifespan." Proceedings of the National Academy of 
% Sciences 111.46 (2014): E4997-E5006.
%
% NOTES: in original analysis negative z-values were set to zero.
% Within-system connectivity was calculated as the mean node-to-node 
% z-value of all nodes of that system to each other. 
% Between-system connectivity was calculated as the mean node-to-node 
% z-value between each node of a system and all nodes of all other systems. 
% ------------------------------------------------------------------------
clear all
% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
dataLoc = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/corrmats/';
output_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
atlas = 'Seitzman300';
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);% load atlas that contains roi info (including which rois belong to each network) 
for net = 1:14
    net_size(net) = length(atlas_params.mods{1,net});
end
weights = net_size./(300-net_size(1));
% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
subs = {'LS02', 'LS03', 'LS05'};
sessions = [8, 10, 10];
%ses_post = [5,5,5];

% ------------------------------------------------------------------------
%% BEGIN ANALYSIS
% ------------------------------------------------------------------------
for sub = 1:numel(subs)
    sub_struct = {};
    all_within = [];
    all_between = [];
    %ses_SI = [];
    for ses = 1:sessions(sub)
        fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', dataLoc, subs{sub}, subs{sub}, ses);
        mat_struct = load(fname);
        matrix = mat_struct.corrmat; clear mat_struct
        matrix = single(FisherTransform(matrix));% fisher transform r values
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
            ses_SI(sub,ses) = (mean(all_within) - mean(all_between))/mean(all_within);                      
    end
end
save([output_dir 'allsubs_seg_index_ses_longitudinal.mat'], 'ses_SI')

LS02_pre = ses_SI(1,1:3);
LS03_pre = ses_SI(2,1:5);
LS05_pre = ses_SI(3,1:5);
LS02_post = ses_SI(1,4:8);
LS03_post = ses_SI(2,6:10);
LS05_post = ses_SI(3,6:10);

vals = [LS02_pre LS02_post LS03_pre LS03_post LS05_pre LS05_post];

for val1 = 1:length(vals)
    for val2 = 1:length(vals)
        diff(val1,val2) = abs(vals(val1)-vals(val2));
    end
end


figure('Position',[1 1 1000 800]);
imagesc(diff,[0 .12]); colormap('jet');
hline_new([0,8,18,28]+0.5,'k',2);
hline_new([3,13,23]+0.5,'k',.5);
vline_new([0,8,18,28]+0.5,'k',2);
vline_new([3,13,23]+0.5,'k',.5);
set(gca,'XTick',[4, 13, 23], 'YTick',[4, 13, 23], 'XTickLabel',...
     {'LS02', 'LS03', 'LS05'}, 'YTickLabel', {'LS02', 'LS03', 'LS05'});
axis square;
colorbar;
%title('Correlation Matrix Similarity');
saveas(gcf,[output_dir 'SimilarityMat_SegregationIndex_Longitudinal.tiff'],'tiff');
