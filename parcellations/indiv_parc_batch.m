


%subs = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
subs = {'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};%;LS02'
edgemap = 0; %edgemap exists
for s = 1:numel(subs)
    %addpath(genpath('/scratch/dcr8536/parcellations'))
    dconnfile = ['/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/sub-' subs{s} '_allsess_tmasked.dconn.nii'];
    surfdir = ['/projects/b1081/Lifespan/Nifti/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/sub-' subs{s} '/NativeVol/fsaverage_LR32k/'];
    outputdir = ['/scratch/dcr8536/parcellations/sub-' subs{s}];
    if ~exist(outputdir)
        mkdir(outputdir);
    end
    if edgemap == 0        
        surface_parcellation_singlesub(['sub-' subs{s}], dconnfile, surfdir, 10, 0, outputdir);
    end
    edgeciftiname = [outputdir '/corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii'];
    filestem = ['sub-' subs{s} '_individual_parcels'];
    threshperc = .50;
    template = '/projects/b1081/member_directories/dperez/template.dtseries.nii';
    parcel_creator_cifti(edgeciftiname,filestem,threshperc,template)    
end

