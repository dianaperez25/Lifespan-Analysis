%% Script to check if fmriprep is correct variable type

%datadir = '/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.2.0/fmriprep';
datadir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/29Nov2020_transfer/preproc_fmriprep-20.2.0/fmriprep';
sub = 'LS03';
ses = [1:5];
runs = [9,9,11,8,9];

for s = 1:length(ses)
    for r = 1:runs(s)
        filename = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz', datadir, sub, ses(s), sub, ses(s), r);
        load_nii(filename);
    end
end
        