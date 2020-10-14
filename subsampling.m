clear all

dataDir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_fmriprep-20.0.6/fmriprep';
subject = 'LS03';
sessions = 5;
runs = [9,9,11,8,9];
catData = [];
catTmask = [];
pts2sample = 8181; %roughly equivalent to 150 minutes

for i = 1:sessions
    for j = 1:runs(i)
        
        rest_run = load_untouch_nii_wrapper(sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz',dataDir, subject, i, subject, i, j));
        tmask = table2array(readtable(sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%02d_desc-tmask_fFD.txt',dataDir, subject, i, subject, i, j)));
        tmask = logical(tmask);
        masked_data = rest_run(:,tmask);
        
%        catData = [catData masked_data];
        catData = [catData rest_run];
        catTmask = [catTmask tmask];
    end
end

disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject, size(catData,1), size(catData,2)))

restmask = catTmask;
samplepts = find(restmask == 1);
restmask(samplepts) = 0;
sampleselect = datasample(catData, samplepts, 'Replace', false);
restmask(sampleselect) = 1;
data_subsample = catData(1:voxnum, logical(restmask));



