function template_matching_parcels(subject)

%subject = 'LS02';
sessions = 5;
runs = 13; % I think this is the max number of runs
% some paths
cifti_data_dir = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/';
tmask_dir = '/projects/b1081/Lifespan/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
%ind_parcels_dir = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/indiv_parcels/';
ind_parcels_dir = '/scratch/dcr8536/Lifespan_KongParcellation/';
output_dir = ['/scratch/dcr8536/parcellations/sub-' subject '/'];
% load network templates
load('/projects/b1081/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates

templates = templates(1:59412,:)';
template_values_sorted = sort(templates(:), 'descend');
threshval= template_values_sorted(round(numel(template_values_sorted) .* 0.05));
threshtemplates= templates >= threshval;
clear allTvals allTvals_sorted templates threshval

% Set variables
if ~exist('output_dir')
    output_dir = pwd;
end

% load individual parcellation file for subject
%ind_parcels_fname = sprintf('%s/sub-%s/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', ind_parcels_dir, subject, subject);
ind_parcels_fname = sprintf('%s/%s_cMSHBM.dtseries.nii', ind_parcels_dir, subject);
ind_parcels = ft_read_cifti_mod(ind_parcels_fname);
unique_parcels = unique(ind_parcels.data);
if unique_parcels(1) == 0
    unique_parcels(1) = [];
end

parcel_corrmap_fname = [output_dir 'sub-' subject '_corrmap_by_parcel_KONG.mat'];
output_fname = sprintf('%s/sub-%s_corr_parcel_by_template_binarized_KONG.mat', output_dir, subject);

if exist(parcel_corrmap_fname)
    load(parcel_corrmap_fname)
else
%initialize some variables
cifti_ts_concat = []; tmask_concat = [];
for ses = 1:sessions
    for run = 1:runs
        %load cifti timeseries
        cifti_ts_fname = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii', cifti_data_dir, subject, ses, subject, ses, run);  
        tmask_fname = sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%d_desc-tmask_fFD.txt', tmask_dir, subject, ses, subject, ses, run);
        if exist(cifti_ts_fname)
            cifti_ts = ft_read_cifti_mod(cifti_ts_fname);
            tmask = table2array(readtable(tmask_fname));
            cifti_ts_concat = [cifti_ts_concat cifti_ts.data];
            tmask_concat = [tmask_concat tmask'];
        end 
    end
end

masked_data = cifti_ts_concat(1:59412,logical(tmask_concat));

    % Loop through each parcel and calculate its connectivity map
    for ind=1:length(unique_parcels)
        parcel_verts = find(ind_parcels.data==unique_parcels(ind));
        corr_map = paircorr_mod(mean(masked_data(parcel_verts,:))',masked_data'); 
        corr_map(isnan(corr_map)) = 0;
        corr_map_parcel(:,ind) = corr_map';
        clear corr_map
    end
    clear masked_data
    save(parcel_corrmap_fname, 'corr_map_parcel', '-v7.3');
end
    
    % Match each variant's connectivity map to each template connectivity map
    disp('Matching to templates and computing goodness of fit')
    corr_coeff = zeros(length(unique_parcels),size(threshtemplates,1));
    
    corrmap_linear = corr_map_parcel(:);
    sorted_vals = sort(corrmap_linear, 'descend');
    subject_thresh = sorted_vals(round(0.05 * numel(sorted_vals)));
    corr_mat_thresh = corr_map_parcel >= subject_thresh;
    for parcel=1:length(unique_parcels)
        disp(['Subject ' subject ', parcel #' num2str(parcel) ' out of ' num2str(length(unique_parcels))]);
        for templatenum = 1:size(threshtemplates,1) %% this isn't doing dice correlation, is it? It's just straight up correlation....
            corr_coeff(parcel,templatenum) = paircorr_mod(corr_mat_thresh(:,parcel),threshtemplates(templatenum,:)');
        end
    end
    
%     for parcel=1:length(unique_parcels)
%         disp(['Subject ' subject ', parcel #' num2str(parcel) ' out of ' num2str(length(unique_parcels))]);
%         sorted_vals = sort(corr_map_parcel(:,parcel), 'descend');
%         subject_thresh = sorted_vals(round(0.05 * numel(sorted_vals)));
%         corr_mat_thresh = corr_map_parcel(:,parcel) >= subject_thresh;
%         for templatenum = 1:size(threshtemplates,1) %% this isn't doing dice correlation, is it? It's just straight up correlation....
%             corr_coeff(parcel,templatenum) = paircorr_mod(corr_mat_thresh,threshtemplates(templatenum,:)');
%         end
%     end
    
    % Determine the 1st and 2nd place template network (max and next max match)
    corr_coeff(isnan(corr_coeff)) = 0;
    [~,maxi] = max(corr_coeff,[],2);
    tempCorrCoeff = corr_coeff;
    for qq=1:length(unique_parcels)
        tempCorrCoeff(qq,maxi(qq))=0;
    end
    [~,nextMax] = max(tempCorrCoeff,[],2);
    clear tempCorrCoeff
    
    
% Save the correlation between each parcel and each template network
% (between their connectivity maps)
save(output_fname, 'corr_coeff', 'maxi', 'nextMax', '-v7.3');

%now make a cifti network map
%load a template
template_fname = '/projects/b1081/member_directories/dperez/template.dtseries.nii';
template = ft_read_cifti_mod(template_fname);

%loop through each parcel and label its vertices with the right network
% but first, adjust indices to add the "nothing" networks
color_change = [1:14; 1:3, 5, 7:16];
parcel_net_assign = raw2colors_mat(maxi,color_change');
network_map = zeros(size(template.data));
for ind=1:length(unique_parcels)
    parcel_verts = find(ind_parcels.data==unique_parcels(ind));
    network_map(parcel_verts) = parcel_net_assign(ind);
end
template.data = network_map;

netmap_fname = sprintf('%s/sub-%s_indiv_parcels_net_assigned_binarized_KONG.dtseries.nii', output_dir, subject);
ft_write_cifti_mod(netmap_fname, template)
clear maxi corr_coeff nextMax 
end
