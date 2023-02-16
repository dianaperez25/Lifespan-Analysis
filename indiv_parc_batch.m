



subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};

for s = 1:numel(subs)
    cd /projects/b1081/member_directories/dperez/parcellations
    dconnfile = ['/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/Lifespan_Dconns/sub-' subs{s} '_allsess_tmasked.dconn.nii'];
    surfdir = ['/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-' subs{s} '/NativeVol/fsaverage_LR32k/'];
    outputdir = ['/projects/p31161/parcellations/sub-' subs{s}];
    if ~exists(outputdir)
        mkdir(outputdir);
    end
    surface_parcellation_singlesub('sub-LS03', dconnfile, surfdir, 100, 0, outputdir);
end