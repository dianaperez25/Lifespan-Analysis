% make SNR maps

subject = 'LS03';
sessions = 3;
data_dir = '/scratch/dcr8536/surf_testing/Nifti/derivatives/preproc_FCProc';
for ses = sessions
    volfile = [data_dir '/sub-' subject '/ses-' num2str(ses) '/func/sub-' subject '_ses-' num2str(ses) '_task-rest_desc-mode1000_mean'];
    map_vol_to_surface_Lifespan_specific(volfile,subject)
    %make into a cifti
    cd /projects/b1081/member_directories/dperez/gifti-master/@gifti/private/
    gifti_tmp = []; gifti_tmp.metaData = []; gifti_tmp.label = []; gifti_tmp.data = [];
    gifti_L = gifti_read([volfile '_L_dil10_32k_fs_LR.func.gii'], gifti_tmp);
    gifti_R = gifti_read([volfile '_R_dil10_32k_fs_LR.func.gii'], gifti_tmp);
    template = ft_read_cifti_mod('/projects/b1081/member_directories/dperez/template.dtseries.nii');
    tmp_cifti = [gifti_L.data{1,1}.data; gifti_R.data{1,1}.data];
    tmp_cifti = tmp_cifti(template.brainstructure>0);
    template.data = tmp_cifti;
    ft_write_cifti_mod([volfile '.dtseries.nii'], template)
end

