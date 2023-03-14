addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master/'))
L_thickness_old=gifti('/projects/b1081/Longitudinal_iNetworks/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-INET006post/NativeVol/fsaverage_LR32k/sub-INET006post.L.thickness.32k_fs_LR.shape.gii');
L_thickness_new=L_thickness_old;
L_thickness_new.cdata=L_thickness_old.cdata .*-1;

save(L_thickness_new,'/projects/b1081/Longitudinal_iNetworks/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-INET006post/NativeVol/fsaverage_LR32k/sub-INET006post.L.thickness_manuallyinverted.32k_fs_LR.shape.gii');
