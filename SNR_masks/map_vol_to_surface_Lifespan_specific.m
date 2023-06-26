function map_vol_to_surface_MSCspecific(volfile,subject)
% volfile should be without extension

%volfile = '1_L_Hand+2_R_Hand-3_L_Leg-4_R_Leg_zstat_333_t88';

% path to goodvoxels - use the one from the first day. For this purpose
% (just mapping), just use the first one in the folder. Could make a union
% mask eventually to be more accurate
goodvoxdir = ['/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/sub-' subject '/ses-1/goodvoxels/'];
goodvoxmask = [goodvoxdir 'sub-' subject '_ses-1_task-rest_goodvoxels.nii.gz'];
if ~exist(goodvoxmask)
    % calculate goodvoxels for session 1 only
    goodvox_fname = goodvoxels_wrapper(subject,s,tmask_names,preprocdata_names,surffuncdir,allstart_fstring2,run_nums,fs_LR_surfdir,goodvoxfolder,...
        T1name_end,space_short,force_remake_concat,force_ribbon,force_goodvoxels);
end

surfdir = '/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';
%timdir = '/data/nil-bluearc/GMT/Laumann/MSC'; % DP: not sure what this is
maskdir = '/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti_masks/subcortical_mask_LR_333_MNI.nii.gz';
medial_mask_L = '/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/projects/b1081/Scripts/CIFTI_RELATED/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';
workbenchdir = '/projects/b1081/Scripts/workbench2/bin_linux64/';

% resample to the surface using the subject's own pial/white/etc. data
HEMS = {'L','R'};
for hem = 1:2

    midsurf = [surfdir '/sub-' subject '/MNI/Native/sub-' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
    midsurf_LR32k = [surfdir '/sub-' subject '/MNI/fsaverage_LR32k/sub-' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    whitesurf = [surfdir '/sub-' subject '/MNI/Native/sub-' subject '.' HEMS{hem} '.white.native.surf.gii'];
    pialsurf = [surfdir '/sub-' subject '/MNI/Native/sub-' subject '.' HEMS{hem} '.pial.native.surf.gii'];
    nativedefsphere = [surfdir '/sub-' subject '/MNI/Native/sub-' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
    outsphere = [surfdir '/sub-' subject '/MNI/fsaverage_LR32k/sub-' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];

    system([workbenchdir '/wb_command -volume-to-surface-mapping ' volfile '.nii.gz ' midsurf ' ' volfile '_' HEMS{hem} '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' goodvoxmask]);
    system([workbenchdir '/wb_command -metric-dilate ' volfile '_' HEMS{hem} '.func.gii ' midsurf ' 10 ' volfile '_' HEMS{hem} '_dil10.func.gii']);
    system([workbenchdir '/wb_command -metric-resample ' volfile '_' HEMS{hem} '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' volfile '_' HEMS{hem} '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
    
    % remove intermediate files
    delete([volfile '_' HEMS{hem} '.func.gii']);
    delete([volfile '_' HEMS{hem} '_dil10.func.gii']);

end

end

