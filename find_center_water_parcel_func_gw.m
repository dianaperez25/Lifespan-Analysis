function [center_coords indpos] = find_center_water_parcel_func(waterfile,surfdir,outputfile,hems)
%TOL

HEMS = {'L';'R'};

switch hems
    case 'LEFT'
        h = 1;
    case 'RIGHT'
        h = 2;
end

for hem = h
    clear indpos
    sphere = gifti([surfdir '/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);

    % CG - was in Tim's original script, btu seems unused?
    %[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    
    water = gifti(waterfile);
    water = water.cdata;

    %Mask watersheds  %% commented this out of the distance file because
    %the gifti I used does not have the medial wall
    maskname = ['/gpfs/research/grattonlab/member_directories/gwulfekuhle/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    water(logical(mask)) = 0;
    
    
    waternum = unique(water);
    waternum(waternum==0) = [];
    
    for w = 1:length(waternum)
        
        ind = find(water==waternum(w));
           
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        
        
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        indpos(w) = ind(indval);
        
    end
    
    metric = zeros(32492,1);    
    metric(indpos) = 1;
    
    %Save out midthickness coordinates of centers
    midthickfile = dir([surfdir '/*' HEMS{hem} '*midthickness.32k_fs_LR.surf.gii']);
    midthickname = midthickfile(1).name;
    midthick = gifti([surfdir '/' midthickname]);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    save(gifti(single(metric)),[outputfile '_' HEMS{hem} '_parcel_center.func.gii'])
    dlmwrite([outputfile '_' HEMS{hem} '_parcel_center.txt'],indpos,'delimiter','\n')


end


