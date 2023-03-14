function make_sulc_mask_HCP()
addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master
addpath /projects/b1081/Scripts/CIFTI_RELATED/Resources
addpath /projects/b1081/Scripts/Create_Dconns
addpath /projects/b1081/Scripts/Scripts_general/FCprocessing 
%addpath /projects/b1081/Scripts/Scripts_general/Workbench_Scripts
%function called gifti will need to addpath

%directory
%Alexis: I think the location for our project is here 
% Need to ask Ally about 32K for post and pre sessions
fsDir = '/projects/b1081/Longitudinal_iNetworks/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';
%topDir we may want to have within your Arianas_BIDS folder
topDir = '/projects/b1081/Longitudinal_iNetworks/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/Ariana_BIDS/';
%making cifti template
dataFile = strcat('/projects/b1081/member_directories/aporter/threshold_CG/130730_RW2111_binarySpatialCorrMap_thresh10%.dtseries.nii');
cifti= ft_read_cifti_mod(dataFile);


%topDir = '/home/cgv5452/for_hubsR01/';
%fsDir = [topDir 'fsaverage_LR32k/']; %'/projects/b1081/MSC/MSCdata_v1/FREESURFER_fs_LR/';
parcellation_label_file.L = '/projects/b1081/Atlases/32k_ConteAtlas_v2/parcellations_VGD11b.L.32k_fs_LR.label.gii';
parcellation_label_file.R = '/projects/b1081/Atlases/32k_ConteAtlas_v2/parcellations_VGD11b.R.32k_fs_LR.label.gii';
outDir = [topDir 'HCP_sulc_masks/']; % '/projects/b1081/MSC/Atlases/IndParcels/PCcalcs/sulc_masks/';
if ~exist(outDir)
    mkdir(outDir);
end

% subjects
% get subject list
%Alexis: I fixed this one so it will make a list of your subjects 
fileinfo = dir('/projects/b1081/Longitudinal_iNetworks/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub*');

%fileinfo = dir([topDir 'HCP_Hub_IDs_forCaterina/*Hub_cluster_IDs*']);
for i = 1:length(fileinfo)
    subjects{i} = fileinfo(i).name();
end


% some constants to set
sulc_thresh = 0; % look for things above the mean in depth
medial_index = 3; % third map in parcellations file

for s = 1:length(subjects)
    sub = subjects{s};
    %hemis = {'lh','rh'};
    hemis = {'L','R'};
    for h = 1:2
        hemi = hemis{h};
        
        fssubDir = [fsDir sub '/NativeVol/fsaverage_LR32k/']; %same as NativeVol_fs_LR version
        %Alexis: adjusting this to inetworks filename
        %sulc_fname = [fssubDir hemi '.sulc'];
        %first step change this
        %'/projects/b1081/Longitudinal_iNetwroks/Nifti/derivatives/freesurfer-6.0.1/sub-INET001post/surf'
        %to 32k space
        sulc_fname = [fssubDir sub '.' hemi '.sulc.32k_fs_LR.shape.gii'];
        
        % The sulc measures used here are normalized... use toolbox below for
        % more exact mm based values
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6560124/
        % at a glance, tho, setting thresh to 0 does a good job of getting
        % crests of gyrim
        
        try
            sulc = gifti(sulc_fname);
        catch
            disp(sub);
        end
        sulc_data.(hemi) = sulc.cdata;
        
        % create a masked version, thresholded
        sulc_data_mask.(hemi) = sulc_data.(hemi) > sulc_thresh;
        
        % load medial wall mask
        parcellation_label = gifti(parcellation_label_file.(hemi));
        medial_wall.(hemi) = parcellation_label.cdata(:,medial_index);
    end
    
    % code based on attribute_gifti_data_to_cifti
    sulc_fname = [outDir sub '.sulc.32k_fs_LR.dtseries.nii'];
    cifti_order = [sulc_data.L(~medial_wall.L); sulc_data.R(~medial_wall.R)];
    %cifti_order(size(cifti_order,1):66697,1) = 0; % make everything else (subcortex?) zero
    cifti_write_wHDR(cifti_order,dataFile, sulc_fname, 'dtseries'); %to our future selves. Do this.
    
    sulc_mask_fname = [outDir sub '.sulc_MASKED.32k_fs_LR.dtseries.nii'];
    cifti_order = [sulc_data_mask.L(~medial_wall.L); sulc_data_mask.R(~medial_wall.R)];
    cifti_order(size(cifti_order,1):66697,1) = 0; % make everything else (subcortex?) zero
    cifti_write_wHDR(cifti_order,dataFile, sulc_mask_fname, 'dtseries');
end