addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master'))
data_dir = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';
out_dir = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/thickness/'
subject = 'LS02';
hems = {'L','R'};
for h = 1:numel(hems)
    pial_surf = ft_read_cifti_mod(sprintf('%s/sub-%s/MNI/fsaverage_LR32k/sub-%s.%s.pial.32k_fs_LR.surf.gii',data_dir, subject, subject, hems{h}));
    wm_surf = ft_read_cifti_mod(sprintf('%s/sub-%s/MNI/fsaverage_LR32k/sub-%s.%s.white.32k_fs_LR.surf.gii',data_dir, subject, subject, hems{h}));
    cort_thickness = wm_surf.data-pial_surf.data;
    wm_surf.data = cort_thickness;
    ft_write_cifti_mod(sprintf('%s/sub-%s.%s.calculated_thickness.32k_fs_LR.surf.gii',out_dir, subject, hems{h}), wm_surf);
end