addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master
read_file_pre=ft_read_cifti_mod(['/projects/b1081/member_directories/aporter/Ariana_Captures/sub-INET001_diff_spCorr_vs_WU120.dtseries.nii']);
mask = double(logical(read_file_pre.data==0));
read_file_pre.data(mask,1)= NaN;
ft_write_cifti_mod(['/projects/b1081/member_directories/aporter/Ariana_Captures/sub-INET001_diffmask_spCorr_vs_WU120.dtseries.nii'],read_file_pre);
