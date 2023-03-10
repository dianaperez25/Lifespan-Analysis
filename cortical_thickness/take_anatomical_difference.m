addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master
read_file_pre=ft_read_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-114/NativeVol/fsaverage_LR32k/sub-114.sulc.32k_fs_LR.dtseries.nii']);
read_file_post=ft_read_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-INET000/NativeVol/fsaverage_LR32k/sub-INET000.sulc.32k_fs_LR.dtseries.nii']);
difference=read_file_post.data-read_file_pre.data;
read_file_post.data=difference;
ft_write_cifti_mod(['/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-114/NativeVol/fsaverage_LR32k/sub-INET000_difference.sulc.32k_fs_LR.dtseries.nii'],read_file_post);