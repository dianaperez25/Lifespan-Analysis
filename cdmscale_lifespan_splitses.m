% Running cmdscale_mat with Lifespan subs, split all sess data into two
clear all
% make the input matrix
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
colors = [1, 0, 0; 0, 1, 0; 0, 0, 1; 0, 1, 1; 1, 0, 1; 0, 0, 0; 0.4940, 0.1840, 0.5560; 0.9290, 0.6940, 0.1250];
cmat_loc = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
atlas = 'Seitzman300';
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
input = [];
count = 1;
groups = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8];
concat_data = [];
%atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);
%roi_data = load_nii_wrapper(atlas_params.MNI_nii_file); %vox by 1


%% add code to separate by session, get minimum amt of data per session
for sub = 1:numel(subs)  
    for ses = 1:5
        load([cmat_loc 'sub-' subs{sub} '/sub-' subs{sub} '_sess-' num2str(ses) '_task-rest_corrmat_Seitzman300.mat'])
        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));
        allSubs_amtData(ses,sub) = size(masked_data,2);
        concat_data = [concat_data masked_data];
    end
    all_subs_data{sub}=concat_data;
    concat_data = [];
    clear sess_roi_timeseries_concat
end

total_data = sum(allSubs_amtData,1);
smallest_amt = min(total_data);
amt_data = floor(smallest_amt/2);

    
for sub = 1:numel(subs) 
    
    masked_data = all_subs_data{sub};        
    matched_data_1 = masked_data(:,1:amt_data);
    matched_data_2 = masked_data(:,(size(masked_data,2)-amt_data+1):end);
    corrmat_1 = paircorr_mod(matched_data_1');
        %figure_corrmat_GrattonLab(corrmat,atlas_params,-1,1);
    input(:,:,count) = corrmat_1;
    count = count + 1;
    corrmat_2 = paircorr_mod(matched_data_2');
    input(:,:,count) = corrmat_2;
    count = count + 1;
end

figure;
[Y E] = cmdscale_mat_MSC(input,groups,'euclidean',colors,subs)
figure;
[Y E] = cmdscale_mat_MSC(input,groups,'correlation',colors,subs)