addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources
addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master
subject_baseline_path=ft_read_cifti_mod('/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI_20.2.0/Variants/spCorr/sub-INET003_preCOVID_spCorr_vs_WU120.dtseries.nii');
subject_followup_path=ft_read_cifti_mod('/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI_20.2.0/Variants/spCorr/sub-INET003_postCOVID_spCorr_vs_WU120.dtseries.nii');
template=ft_read_cifti_mod('/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI_20.2.0/Variants/spCorr/sub-INET003_preCOVID_spCorr_vs_WU120.dtseries.nii');
pre=subject_baseline_path.data;
post=subject_followup_path.data;
difference=post-pre;
[correlation,pval]=corrcoef(pre, post);
%template.data=difference;
template.data=correlation;
ft_write_cifti_mod('/projects/b1081/member_directories/aporter/sub-INET003_corr_spCorr_vs_WU120.dtseries.nii',template);
