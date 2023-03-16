%% Script to make distance matrix for individualized surfaces
% need the output of this script to make infomap network maps
% option to do vertex- or parcel-wise distance matrix
clear all

%% PATHS
out_dir = '/scratch/dcr8536/distances/'; %where the output files will be stored
wb_path = '/projects/b1081/Scripts/workbench2/bin_linux64/'; % path to wb_command, need this to make geo distance matrices
surf_root = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/'; % path to the individualized surfaces
parcels = 1; % set to 1 to make distance matrix for parcels
parcel_dir = '/scratch/dcr8536/parcellations/';
cifti_tmp = '/scratch/dcr8536/template.dtseries.nii';
gifti_tmp = '/projects/b1081/Atlases/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii';
%make output directory if it doesnt already exist
if ~exist(out_dir)
    mkdir(out_dir);
end

%% VARIABLES
subject = {'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};

for sub = 1:numel(subject)
    sub_ID = ['sub-' subject{sub}];
    disp(['Creating distance matrix for ' sub_ID]);
    
    % some subject specific paths
    surf_dir = sprintf('%s/%s/NativeVol/fsaverage_LR32k/', surf_root, sub_ID);
    L_surf_fname = sprintf('%s/%s.L.midthickness.32k_fs_LR.surf.gii', surf_dir, sub_ID);
    R_surf_fname = sprintf('%s/%s.R.midthickness.32k_fs_LR.surf.gii', surf_dir, sub_ID);
    
    % these are the geo distance matrices 
    output_fname_L = sprintf('%s/%s_distances_left.dconn.nii', out_dir, sub_ID);
    output_fname_R = sprintf('%s/%s_distances_right.dconn.nii', out_dir, sub_ID);
    % if they don't exist, we'll make them now
    if ~exist(output_fname_L)
        command = [wb_path 'wb_command -surface-geodesic-distance-all-to-all ' L_surf_fname ' ' output_fname_L];
        system(command);
    end
    if ~exist(output_fname_R)    
        command = [wb_path 'wb_command -surface-geodesic-distance-all-to-all ' R_surf_fname ' ' output_fname_R];
        system(command);
    end
    
    if parcels    
        % these are gifti files with the individual parcel indices
        left_gii_file = [parcel_dir sub_ID '/' sub_ID '.L.individual_parcels.func.gii']; %need to also make a .gii %scripts/cifti_related/resources/read_write_cifti_gifti_@gifti
        right_gii_file = [parcel_dir sub_ID '/' sub_ID '.R.individual_parcels.func.gii'];
        % if they don't exist, we'll make them now
        if ~exist(left_gii_file) 
            fname_stem = '_individual_parcels_edgethresh_0.5.dtseries.nii'; %the rest of the file name for the parcellation file after subject ID
            make_gifti_file(sub_ID, fname_stem, parcel_dir, cifti_tmp, gifti_tmp, parcel_dir)
        end
        %we need to find the center coordinates and indices for each parcel
        [center_coords_L indpos_L] = find_center_water_parcel_func_gw(sub_ID, left_gii_file,surf_dir,[out_dir sub_ID '_parcel_centers_L'],'LEFT'); %may need to find this script
        [center_coords_R indpos_R] = find_center_water_parcel_func_gw(sub_ID, right_gii_file,surf_dir,[out_dir sub_ID '_parcel_centers_R'],'RIGHT');
        % we put all the coordinates together
        coords_all = [center_coords_L; center_coords_R];
    else
        L_surf = gifti(L_surf_fname);
        R_surf = gifti(R_surf_fname);
        L_surf = L_surf.vertices;
        R_surf = R_surf.vertices;
        coords_all = [L_surf; R_surf];
    end
    
    % now we calculate the euclidean center-to-center distances
    % we'll use these for cross-hemisphere distances
    D = pdist(coords_all,'euclidean');
    dist_mat = squareform(D);
    
    % let's extract the geodesic distances
    L_geo = ft_read_cifti_mod(output_fname_L);
    R_geo = ft_read_cifti_mod(output_fname_R);

    % just a little bit of data manipulation to make sure we are indexing
    % the right values
    if parcels 
        % for parcels, we'll extract the distances for the indices of each
        % parcel center
        L_geo = L_geo.data(indpos_L', indpos_L);
        R_geo = R_geo.data(indpos_R', indpos_R);
        % and we'll specify where the left hemisphere ends, since each
        % subject has a different number of parcels
        Lhem_end_ind = length(indpos_L); % where left hemisphere distances end; note we don't specify where they start because they always start at 1
    else
        % for vertices, we'll use all the distance values, and we'll set
        % the index for where left hemisphere ends
        L_geo = L_geo.data;
        R_geo = R_geo.data;
        Lhem_end_ind = 32492;
    end
        
    % now we replace the within-hemisphere distances in the euclidean
    % distance matrix with geodesic distances
    % top-left quadrant = left-to-left distances
    % bottom-right quadrant = right-to-right distances
    dist_mat(1:Lhem_end_ind,1:Lhem_end_ind)=L_geo;
    dist_mat(Lhem_end_ind+1:end, Lhem_end_ind+1:end) = R_geo;
    
    if parcels 
        save([out_dir subject{sub} '_parcels_distmat.mat'], 'dist_mat', '-v7.3')
    else
        save([out_dir subject{sub} '_vertices_distmat.mat'], 'dist_mat', '-v7.3')
    end
    
    clear dist_mat D coords_all
   
end
function  make_gifti_file(subject, file_stem, file_path, cifti_fname, gifti_fname, output_dir)
    addpath('/projects/p31161/lateralization_code/PerezEtAl_HemAsymmetries')
    cifti_tmp = ft_read_cifti_mod(cifti_fname);
    gifti_tmp = gifti(gifti_fname);

    parcels = ft_read_cifti_mod([file_path subject '/' subject file_stem]);
    new_cifti = insert_nonbrain(parcels.data, 'both', cifti_tmp);
    left_hem = new_cifti(1:32492);
    right_hem = new_cifti(32493:end);
    gifti_tmp.cdata = left_hem;
    save(gifti_tmp, [output_dir subject '/' subject '.L.individual_parcels.func.gii']);
    gifti_tmp.cdata = right_hem;
    save(gifti_tmp, [output_dir subject '/' subject '.R.individual_parcels.func.gii']);

end

