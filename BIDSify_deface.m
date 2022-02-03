
subjects = {'LS05', 'LS08', 'LS11'};
sessions = [5, 3, 5];
root_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/';
nifti_path = [root_dir 'Nifti/'];
config_path = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/config.json';
% cd(nifti_path)

for sub = 1:numel(subjects)
    for ses = 1:sessions(sub)
        output_file = ['sub-' subjects{sub} '_ses-' num2str(ses) '_BIDSify_output.txt'];
        data_path = [root_dir 'DICOMS/sub-' subjects{sub} '/' subjects{sub} '_' num2str(ses) '/'];
        [status, output] = system(['dcm2bids -d ' data_path ' -p ' subjects{sub} ' -s ' num2str(ses) ' -c ' config_path ' -o ' nifti_path]);
        disp(['Data for subject ' subjects{sub} ' session ' num2str(ses) ' has been BIDSified!']);
        save_msg = output;
        if ses == 1
            [status, output] = system(['pydeface sub-' subjects{sub} '_ses-' num2str(ses) '_acq-ADNI3_T1w.nii.gz --outfile sub-' subjects{sub} '_ses-' num2str(ses) '_acq-ADNI3_T1w.nii.gz --force']);
        elseif ses == 2
            [status, output] = system(['pydeface sub-' subjects{sub} '_ses-' num2str(ses) '_acq-iso_run-01_T2w.nii.gz --outfile sub-' subjects{sub} '_ses-' num2str(ses) '_acq-iso_run-01_T2w.nii.gz --force']);
            [status, output] = system(['pydeface sub-' subjects{sub} '_ses-' num2str(ses) '_acq-iso_run-02_T2w.nii.gz --outfile sub-' subjects{sub} '_ses-' num2str(ses) '_acq-iso_run-02_T2w.nii.gz --force']);            
        elseif ses == 3
            [status, output] = system(['pydeface sub-' subjects{sub} '_ses-' num2str(ses) '_acq-RMS_T1w.nii.gz --outfile sub-' subjects{sub} '_ses-' num2str(ses) '_acq-RMS_T1w.nii.gz --force']);
        end 
        save_msg = [save_msg; output];
        fid = fopen(output_file, 'wt');
        fprint(fid, save_msg);
        fclose(fid);
    end
end

        
