outdir = '/scratch/dcr8536/distances';
surfdir = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';
wb_path = '/projects/b1081/Scripts/workbench2/bin_linux64/';


% geodesic distance matrices
load([surfdir '/geodesic_distance.L.32k_fs_LR.mat']);
load([surfdir '/geodesic_distance.R.32k_fs_LR.mat']);

subjects = 1;

for subs = subjects
    
    %subject = sprintf('MSC%02d',subjects(s));
    subject = 'MSC01';
    disp(subject);
    
    orig_parceldir = ['/scratch/dcr8536/parcellations/'];
    left_gii_file = [orig_parceldir 'sub-' subject '_400_L.func.gii']; %need to also make a .gii %scripts/cifti_related/resources/read_write_cifti_gifti_@gifti
    right_gii_file = [orig_parceldir subject '_400_R.func.gii'];
    
    
    [center_coords_L indpos_L] = find_center_water_parcel_func_gw(left_gii_file,surfdir,[outdir subject '_parcel_centers_L'],'LEFT'); %may need to find this script
    [center_coords_R indpos_R] = find_center_water_parcel_func_gw(right_gii_file,surfdir,[outdir subject '_parcel_centers_R'],'RIGHT');
    
    % euclidean center-to-center distances
    center_coords_all = [center_coords_L; center_coords_R];
    D = pdist(center_coords_all,'euclidean');
    Z = squareform(D);
    
    % and now replace with geodesic for all within hemisphere comparisons
    Z_geo = Z;
    Z_geo(1:length(indpos_L),1:length(indpos_L)) = distmat_surf_L(indpos_L,indpos_L); %need this function 
    Z_geo(length(indpos_L)+1:length(indpos_R)+length(indpos_L),length(indpos_L)+1:length(indpos_R)+length(indpos_L)) = distmat_surf_R(indpos_R,indpos_R);
    
    save([outdir subject '_indparcel_dists_mask.mat'],'Z','Z_geo');
    
    clear center_coords_* D Z Z_geo;
    
end