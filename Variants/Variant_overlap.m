%% Calculate variant overlap

clear all

% add dependencies
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
varmap1_loc = '/Users/dianaperez/Desktop/LS03_matcheddata_REST_Variant_Size_thresh-5_smooth_2.55.dtseries.nii';
varmap2_loc = '/Users/dianaperez/Desktop/Lifespan/LS03_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii';
outfile = '/Users/dianaperez/Desktop/LS03_longitudinal_variant_overlap.dtseries.nii';

% load the variant maps
varmap_1 = ft_read_cifti_mod(varmap1_loc);
varmap_2 = ft_read_cifti_mod(varmap2_loc);
overlap_map = [];


for v = 1:59412
    if varmap_1.data(v) > 0 && varmap_2.data(v) == 0
        overlap_map(v) = 1;
    elseif varmap_2.data(v) > 0 && varmap_1.data(v) == 0
        overlap_map(v) = 2;
    elseif varmap_2.data(v) > 0 && varmap_1.data(v) > 0
        overlap_map(v) = 3;
    else
        overlap_map(v) = 0;
    end
end

template = varmap_1;
template.data = overlap_map;

ft_write_cifti_mod(outfile, template);