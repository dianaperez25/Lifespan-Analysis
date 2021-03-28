datadir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_FCProc';
outdir = '/Users/diana/Desktop/Research';
if ~exist(outdir)
    mkdir(outdir)
end
sub = 'LS05';
ses = [1:5];
%runs = [9,9,11,8,9]; %LS03
runs = [8, 8, 8, 9, 9]; %LS05
%runs = [7, 8]; %LS07
load('/Users/diana/Desktop/Research/template.mat')

for s = 1:length(ses)    
    for r = 1:runs(s)
        fname = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%d_fmriprep_zmdt_resid_ntrpl_bpss_zmdt.nii.gz',datadir,sub,ses(s),sub,ses(s),r);
        %filename = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz',datadir,sub,ses(s),sub,ses(s),runs(r));
        nifti = load_nii_wrapper(fname);
        vox = [];
        new_nii = [];
        for v = 1:902629
          vox = nifti(v,:)';
          if max(vox) > 75 || max(vox) < -75
              new_nii(v,1) = 1;
          else new_nii(v,1) = 0;
          end
        end
        template = tmp;
        img_out = reshape(new_nii, [91 109 91]);
        template.img = img_out;        
        outname = [outdir '/sub-' sub '_ses-' num2str(ses(s)) '_run-' num2str(r) '_outlier_voxels.nii.gz'];
        template.fileprefix = outname;
        template.hdr.dime.dim(5) = 1;
        save_nii(template, outname)
    end
end