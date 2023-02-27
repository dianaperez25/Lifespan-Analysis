%% Script to make distance matrix for individualized surfaces
% need the output of this script to make infomap network maps

%% PATHS
out_dir = '/Users/dianaperez/Desktop/distances/';
wb_path = '/Applications/workbench/bin_macosx64/';
surf_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Pre-COVID/BIDS/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';

if ~exist(out_dir)
    mkdir(out_dir);
end


%% VARIABLES
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};

for sub = 1:numel(subject)
    
    disp(['Creating distance matrix for subject ' subject{sub}]);
    
    % we don't need to load parcels, but do need to load the surfaces
    % geodesic distance matrices
    L_surf_fname = sprintf('%s/sub-%s/NativeVol/fsaverage_LR32k/sub-%s.L.midthickness.32k_fs_LR.surf.gii', surf_dir, subject{sub}, subject{sub});
    R_surf_fname = sprintf('%s/sub-%s/NativeVol/fsaverage_LR32k/sub-%s.R.midthickness.32k_fs_LR.surf.gii', surf_dir, subject{sub}, subject{sub});
    L_surf = gifti(L_surf_fname);
    R_surf = gifti(R_surf_fname);
    L_surf = L_surf.vertices;
    R_surf = R_surf.vertices;
    
    
    % euclidean center-to-center distances
    coords_all = [L_surf; R_surf];
    D = pdist(coords_all,'euclidean');
    Z = squareform(D);
    
    % and now replace with geodesic for all within hemisphere comparisons
    output_fname_L = sprintf('%s/sub-%s_distances_left.dconn.nii', out_dir, subject{sub});
    command = [wb_path 'wb_command -surface-geodesic-distance-all-to-all ' L_surf_name ' ' output_fname_L];
    system(command);
    output_fname_R = sprintf('%s/sub-%s_distances_right.dconn.nii', out_dir, subject{sub});
    command = [wb_path 'wb_command -surface-geodesic-distance-all-to-all ' R_surf_name ' ' output_fname_R];
    system(command);
    
    L_geo = ft_read_cifti_mod(output_fname_L);
    R_geo = ft_read_cifti_mod(output_fname_R);
    
%     Z_geo = Z;
%     Z_geo(1:length(indpos_L),1:length(indpos_L)) = distmat_surf_L(indpos_L,indpos_L); %need this function 
%     Z_geo(length(indpos_L)+1:length(indpos_R)+length(indpos_L),length(indpos_L)+1:length(indpos_R)+length(indpos_L)) = distmat_surf_R(indpos_R,indpos_R);
%     
%     save([out_dir subject '_vertices_distmat.mat'],'Z','Z_geo');
    
    %clear center_coords_* D Z Z_geo;
    
end

