%% script to check data for completeness
clear all

subject = 'LS11';
session = 3;

dcm_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/DICOM';

cd dcm_dir

subject_dir = ['sub-' subject '/' subject '_' num2str(session) '/'];

d = dir(subject_dir);

for x = 1:length(d)
    if contains(d(x).name,'FMRI_')
        cd([subject_dir d(x).name])
        continue
    end
end

d = dir;
rest_dirs = [];
count = 1;
for x = 1:length(d)
    if contains(d(x).name,'MB4_REST') && ~contains(d(x).name,'PHYSIOLOG')
        rest_dirs{count,1} = d(x).name;
        rd = dir(d(x).name);
        rest_dirs{count,2} = length(rd);
        clear rd
        count = count + 1;
    end
end

disp(['Number of rest runs is ' num2str(length(rest_dirs))])




