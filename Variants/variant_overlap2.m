%% variant overlap

subs = {'LS02', 'LS03', 'LS05'}%, 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
data_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/Variant_Maps/';
out_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/';
corrs = [];
for sub = 1:numel(subs)
    vars1 = ft_read_cifti_mod([data_dir subs{sub} '_allsess_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii']);
    %vars2 = ft_read_cifti_mod(['/Users/dianaperez/Desktop/' subs{sub} '_second-half_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii']);
    vars2 = ft_read_cifti_mod(['/Users/dianaperez/Desktop/' subs{sub} '_allsess_timeA_uniqueIDs_variants_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii']);
    out_fname = [out_dir subs{sub} '_variant_overlap_longitudinal.dtseries.nii']
    overlap = zeros(59412,1);
    for vert = 1:59412
        if vars1.data(vert) >0 && vars2.data(vert) >0
            overlap(vert) = 3;
        elseif vars1.data(vert) >0
            overlap(vert) = 1;
        elseif vars2.data(vert) >0
            overlap(vert) = 2;
        end
    end
    template = vars1;
    template.data = overlap;
    corrs(sub) = corr(logical(vars1.data), logical(vars2.data));
    ft_write_cifti_mod(out_fname, template)
end