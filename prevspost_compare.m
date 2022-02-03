%% compare pre and post covid corrmats

sub = 'LS05';
pre_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Pre-COVID/BIDS/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
post_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
pre_sessions = 5;
post_sessions = 5;

for pre_ses = 1:pre_sessions
    load([pre_dir 'sub-' sub '/sub-' sub '_sess-' num2str(pre_ses) '_task-rest_corrmat_Seitzman300.mat'])
    corrmats_pre(:,:,pre_ses) = corrmat; 
end

corr_values_pre = [];
combos = combntns(1:pre_sessions,2);
for pair = 1:length(combos)
    similarity = corrcoef(corrmats_pre(:,:,combos(pair,1)), corrmats_pre(:,:,combos(pair,2)));
    corr_values_pre(pair) = similarity(2);
end
avg_corr_pre = mean(corr_values_pre)

for post_ses = 1:post_sessions
    load([post_dir 'sub-' sub '/sub-' sub '_sess-' num2str(post_ses) '_task-rest_corrmat_Seitzman300.mat'])
    corrmats_post(:,:,post_ses) = corrmat;
end

corr_values_post = [];
combos = combntns(1:post_sessions,2);
for pair = 1:length(combos)
    similarity = corrcoef(corrmats_post(:,:,combos(pair,1)), corrmats_post(:,:,combos(pair,2)));
    corr_values_post(pair) = similarity(2);
end
avg_corr_post = mean(corr_values_post) 

%calculate similarity between timepoints (pre vs post)
corr_values_pre_post = [];
count = 1;
for pre_ses = 1:pre_sessions
    for post_ses = 1:post_sessions
        similarity = corrcoef(corrmats_pre(:,:,pre_ses),corrmats_post(:,:,post_ses));
        corr_values_pre_post(count) = similarity(2);
        count = count + 1;
    end
end
avg_corr_pre_post = mean(corr_values_pre_post)