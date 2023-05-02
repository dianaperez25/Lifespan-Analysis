% create a group average correlation matrix
clear all

dataset = 'iNetworks';
if strcmpi(dataset, 'Lifespan')
    subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
    sessions = 5;
    corrmat_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/preproc_FCProc/corrmats_Seitzman300/';
elseif strcmpi(dataset, 'iNetworks')   
    subs = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};%, 'INET034', 'INET035', 'INET036', 'INET039', 'INET040', 'INET041'};
    sessions = 4;
    corrmat_dir = '/Users/dianaperez/Desktop/iNet_corrmats/';
    %corrmat_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Segregation_analyses/iNetworks/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300';
end
out_dir = '/Users/dianaperez/Desktop/';
fout_str = sprintf('%s/%s_average_correlation_matrix', out_dir, dataset);

count = 1;
for sub = 1:numel(subs)
    for ses = 1:sessions
    % load the matrix
    corrmat_fname = sprintf('%s/sub-%s/sub-%s_sess-%d_task-rest_corrmat_Seitzman300.mat', corrmat_dir, subs{sub}, subs{sub}, ses);
    load(corrmat_fname)
    all_corrmats(count,:,:,:) = corrmat;
    count = count + 1;
    end
end

%calculate mean
mean_corrmat = squeeze(mean(all_corrmats));

%clear some memory
clear all_corrmats count corrmat corrmat_dir sess_roi* tmask*

%load atlas parameters for plotting
atlas_params = atlas_parameters_GrattonLab('Seitzman300','/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/');

%plot mean correlation matrix
figure_corrmat_GrattonLab(mean_corrmat,atlas_params,-1,1);

%save mean correlation matrix image
saveas(gcf,[fout_str '.tiff'],'tiff');
close(gcf);

%save the corrmat data
save([fout_str '.mat'], 'mean_corrmat');

%start by thresholding the group matrix at a set threshold (2% edge density)
thr = 0.02;
adj_mat,adj_mat_sym = gfns.threshold_matrix_density(groupmat,thr)
