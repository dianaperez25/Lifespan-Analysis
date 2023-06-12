addpath(genpath('/Users/dianaperez/Documents/Dependencies'))
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/parcellations/';
subjects = {'LS05', 'LS08', 'LS11','LS14', 'LS16', 'LS17'}; %'LS02', 'LS03',  
parcellation = 'Kong';
fname_stem = 'individual_parcels_edgethresh_0.5';
fname_stem = 'cMSHBM';
for sub = 1:numel(subjects)
    %load necessary files
    parcel_id = ft_read_cifti_mod(sprintf('%s/%s_Parcellation/sub-%s/sub-%s_%s.dtseries.nii', data_dir, parcellation, subjects{sub}, subjects{sub}, fname_stem));
    parcel_nets = ft_read_cifti_mod(sprintf('%s/%s_Parcellation/sub-%s/sub-%s_indiv_parcels_net_assigned_binarized_%s.dtseries.nii', data_dir, parcellation, subjects{sub}, subjects{sub}, parcellation));
    avg_timecourse = load(sprintf('%s/%s_Parcellation/avg_timecourses/sub-%s_individual_parcels_average_timecourses_%s.mat', data_dir, parcellation, subjects{sub}, parcellation));
    %make a matrix with parcel ID and network assignment 
    unique_parcels = unique((parcel_id.data(~isnan(parcel_id.data))));
    if unique_parcels(1) == 0
        unique_parcels(1) = [];
    end
    parcel_netid = [];
    for parcel = 1:length(unique_parcels)
        parcel_netid(parcel,1) = unique_parcels(parcel);
        vertices = find(parcel_id.data==unique_parcels(parcel));
        nets = unique(parcel_nets.data(vertices));
        if length(nets)>1
            disp('Alright dude, something is wrong. This parcel has more than one network asigned to its vertices!')
            nets = parcel_nets.data(vertices);
            unique_nets = unique(nets);
            net_size = [];
            for net = 1:length(unique_nets)
                net_size(net,1) = unique_nets(net);
                net_size(net,2) = length(find(nets==unique_nets(net)));
            end
            [big_net, ind] = max(net_size(:,2));
            nets = net_size(ind);        
        end
        parcel_netid(parcel,2) = nets;
    end
    
    % now I should probably rearrange the corrmats by network
    unique_nets = unique(parcel_netid(:,2));
    sorted_parcels = [];
    for net = 1:length(unique_nets)
        parcels_by_net = find(parcel_netid(:,2)==unique_nets(net));
        sorted_parcels = [sorted_parcels; parcels_by_net];       
    end
    sorted_ts = avg_timecourse.avg_parc_ts(sorted_parcels,:);
    sorted_corrmat = paircorr_mod(sorted_ts');
    imagesc(sorted_corrmat, [-1 1]) % corrmat after sorting, SUCCESS!!
    fout_str = sprintf('%s/%s_Parcellation/avg_timecourses/corrmats/sub-%s_indiv_parcels_corrmats_sorted',data_dir, parcellation, subjects{sub});
    save([fout_str '.mat'],'sorted_corrmat', 'sorted_ts', 'sorted_parcels', 'parcel_netid');
end
        