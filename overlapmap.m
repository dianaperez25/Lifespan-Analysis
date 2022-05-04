%% Overlap Map Making Script
% input: variant maps
% output: a map of proportion of subjects that have overlapping variants at given vertex

clear all
%--------------------------------------------------------------------------
%% PATHS
addpath(genpath('/Users/dianaperez/Documents/Dependencies/cifti-matlab-master'));
root_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/';

%varmap_loc = [root_dir '/Variant_Maps/new_split_vars/reassigned/'];
varmap_str = '_allsess_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii';
varmap_loc = [root_dir '/Variant_Maps/']; %location of variant maps 
 
%template_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/100206_border1ectopic2.dtseries.nii';
template_loc = [varmap_loc 'LS02_allsess_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii'];
out_dir = [root_dir '/figures/'];
if ~exist(out_dir, 'file')
    mkdir(out_dir)
end
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'};

for g = 1:numel(subs)
    cifti = ft_read_cifti_mod([varmap_loc subs{g} varmap_str]);               
    group_map(:,g) = logical(cifti.data);
end
group_overlap = sum(group_map,2)/numel(subs);
template = ft_read_cifti_mod(template_loc);
overlap_map = template;
overlap_map.data = group_overlap;            
ft_write_cifti_mod(out_cifti, overlap_map);
