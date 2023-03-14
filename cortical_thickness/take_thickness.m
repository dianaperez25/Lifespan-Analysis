addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master
read_file_pial=ft_read_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/thickness/sub-114/MNI/fsaverage_LR32k/sub-114.L.pial.32k_fs_LR.surf.gii']);
read_file_white=ft_read_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/thickness/sub-114/MNI/fsaverage_LR32k/sub-114.L.white.32k_fs_LR.surf.gii']);
difference=read_file_white.data-read_file_pial.data;
read_file_white.data=difference;
ft_write_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/thickness/sub-114/MNI/fsaverage_LR32k/sub-114.L.calculated_thickness.32k_fs_LR.surf.gii'],read_file_post);