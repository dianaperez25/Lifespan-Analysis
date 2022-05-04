%% get fd calc averages and amount of data

subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_fmriprep-20.2.0/fmriprep/';

for sub = 1:numel(subs)
    for ses = 1:5
        FD_nums = readtable(sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_desc-framepers_fFD.txt', data_dir, subs{sub}, ses, subs{sub}, ses));
        FD_nums = FD_nums{:,1};
        mean_FD(sub,ses) = mean(FD_nums);
        num_frames = readtable(sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_desc-framenums_fFD.txt', data_dir, subs{sub}, ses, subs{sub}, ses));
        num_frames = num_frames{:,1};
        mean_numframes(sub,ses) = sum(num_frames)*1.1;
    end
end
