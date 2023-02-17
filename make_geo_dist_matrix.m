%% create the distance matrix at the individual level
% necessary for running infomap

subject = {'LS02'};
%wb_path = '/projects/b1081/Scripts/workbench2/bin_linux64/';
wb_path = '/Applications/workbench/bin_macosx64/';
%right_surface = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-LS02/NativeVol/fsaverage_LR32k/sub-LS02.R.midthickness.32k_fs_LR.surf.gii';
%left_surface = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-LS02/NativeVol/fsaverage_LR32k/sub-LS02.L.midthickness.32k_fs_LR.surf.gii';
right_surface = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Pre-COVID/BIDS/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-LS02/NativeVol/fsaverage_LR32k/sub-LS02.L.midthickness.32k_fs_LR.surf.gii';
%left_surface = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-LS02/NativeVol/fsaverage_LR32k/sub-LS02.L.midthickness.32k_fs_LR.surf.gii';
num_verts = size(;
output = [];

for vertex = 1:32492
    system([wb_path ' wb_command -surface-geodesic-distance ' surface ' ' num2str(vertex) ' ' output_fname])
end
