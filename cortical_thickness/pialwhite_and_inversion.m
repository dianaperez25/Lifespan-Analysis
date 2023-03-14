addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti-matlab-master/'))
L_white_coords=gifti('sub-114.L.white.32k_fs_LR.coord.gii');
L_white_coords=L_white_coords.vertices;
L_pial_coords=gifti('sub-114.L.pial.32k_fs_LR.coord.gii');
L_pial_coords=L_pial_coords.vertices;
thickness_matrix=zeros(32492,1);
for vertex=1: length(thickness_matrix)
    thickness_matrix(vertex,1)=pdist2(L_pial_coords(vertex,:), L_white_coords(vertex,:));
end



L_thickness_old=gifti('sub-114.L.thickness.32k_fs_LR.shape.gii');
L_thickness_new=L_thickness_old;
L_thickness_new.cdata=L_thickness_old.cdata .*-1;

save(L_thickness_new,'sub-114.L.thickness_inverted.32k_fs_LR.shape.gii');

scatter(thickness_matrix, L_thickness_new.cdata);

figure;hold on;scatter(thickness_matrix, L_thickness_new.cdata);plot(0:6,0:6);
axis square
xlim([0 7])
ylim([0 7])
plot(0:7,0:7);
ylabel('inverted thickness map after freesurf-to-32k');
xlabel('pial-white after fsLR_32k','interpreter','none');
corr(thickness_matrix, L_thickness_new.cdata);