addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master/'));
L_thickness_old=gifti('/projects/b1081/Longitudinal_Berkeley/Nifti/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-114/NativeVol/fsaverage_LR32k/sub-114.L.thickness_manuallyinverted.32k_fs_LR.shape.gii');
L_thickness_new=gifti('/projects/b1081/Longitudinal_Berkeley/NIFTI/derivatives/freesurfer-6.0.1/subjects/FREESURFER_fs_LR/sub-INET000/NativeVol/fsaverage_LR32k/sub-INET000.L.thickness_manuallyinverted.32k_fs_LR.shape.gii');
corr(L_thickness_old.cdata, L_thickness_new.cdata);