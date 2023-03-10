


subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
edgemap = 0; %edgemap exists
for s = 1:numel(subs)
    addpath(genpath('/scratch/dcr8536/parcellations'))
    dconnfile = ['/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/Lifespan_Dconns/sub-' subs{s} '_allsess_tmasked.dconn.nii'];
    surfdir = ['/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-' subs{s} '/NativeVol/fsaverage_LR32k/'];
    outputdir = ['/scratch/dcr8536/parcellations/subsample_10/sub-' subs{s}];
    if ~exist(outputdir)
        mkdir(outputdir);
    end
    if edgemap == 0
        
        surface_parcellation_singlesub(['sub-' subs{s}], dconnfile, surfdir, 10, 0, outputdir);
    end
    edgeciftiname = [outputdir '/corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii'];
    filestem = ['sub-' subs{s} '_individual_parcels'];
    threshperc = .50;
    template = '/scratch/dcr8536/template.dtseries.nii';
    parcel_creator_cifti(edgeciftiname,filestem,threshperc,template)
    
end

