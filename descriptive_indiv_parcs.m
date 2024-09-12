%% Comparison of properties of individualized networks in older vs younger adults

clear all

%% PATHS
% path to individualized parcellations
indivparc_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/indiv_parcs/';
% path to group average parcellation
groupparc_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
% path to surface area file
surf_areas_path = '/Users/dianaperez/Documents/Dependencies/surf_areas_verts.dtseries.nii';
output_dir = '/Users/dianaperez/Desktop/';
if ~exist(output_dir)
    mkdir(output_dir)
end

% datasets
datasets = {'Lifespan-NU', 'iNet-NU', 'Lifespan-FSU', 'iNet-FSU'};% 
exclude_subs = {'LS46', 'INET108', 'LS08'};
atlas = 'Parcels333';


% load atlas params
% atlas_params = atlas_parameters_GrattonLab(atlas, groupparc_dir);
% load surface area file
surf_areas = ft_read_cifti_mod(surf_areas_path);

for d = 1:numel(datasets)
    num_parcs = []; % create a variable to store number of parcels
    avg_parc_size = []; % create a variable to store average parcel size
    by_subject = [];
    % get subject ID's
    [subjects, sessions, N] = get_subjects(datasets{d}, exclude_subs);
    for sub = 1:numel(subjects)
        sub_parcs = [];

        % load individual parcels
        indiv_parc_file = sprintf('%s/sub-%s/sub-/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', indivparc_dir, subjects{sub}, subjects{sub});
        indiv_parcs = ft_read_cifti_mod(indiv_parc_file);
        
        unique_parcs = unique(indiv_parcs.data);        
        if unique_parcs(1) == 0
            unique_parcs(1) = [];
        end
        
        num_parcs(sub) = length(unique_parcs);% calculate number of parcels (add it to variable)

        
        for parc = 1:num_parcs(sub)
            sub_parcs(parc,1) = unique_parcs(parc);
            parc_verts = find(indiv_parcs.data==unique_parcs(parc));
            sub_parcs(parc,2) = length(parc_verts);
            sub_parcs(parc,3) = sum(surf_areas.data(parc_verts));
        end
        
        avg_parc_size(sub,1) = mean(sub_parcs(:,2));
        avg_parc_size(sub,2) = mean(sub_parcs(:,3));
        by_subject{sub} = sub_parcs; clear sub_parcs parc_verts unique_parcs
    end
    save(sprintf('%s/%s_descriptiveIndivParcelInfo.mat', output_dir, datasets{d}), 'avg_parc_size', 'num_parcs', 'by_subject', 'subjects')
end

% get number of parcels, and average parcel size for group average
group_parcs = ft_read_cifti_mod(sprintf('%s/Evan_parcellation/Published_parcels/Parcels_LR.dtseries.nii', groupparc_dir));
unique_parcs = unique(group_parcs.data);

if unique_parcs(1) == 0
    unique_parcs(1) = [];
end

num_parcs = length(unique_parcs);
for parc = 1:num_parcs
    atlas_parcs(parc,1) = unique_parcs(parc);
    parc_verts = find(group_parcs.data==unique_parcs(parc));
    atlas_parcs(parc,2) = length(parc_verts);
    atlas_parcs(parc,3) = sum(surf_areas.data(parc_verts));
end
        
avg_parc_size(1) = mean(atlas_parcs(:,2));
avg_parc_size(2) = mean(atlas_parcs(:,3));

function [subject, sessions, N] = get_subjects(dataset, exclude_subs)

if strcmpi(dataset, 'Lifespan-NU')
    subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'}; %
    sessions = 5;
elseif strcmpi(dataset, 'iNet-NU')
    subject = {'INET002', 'INET003', 'INET005', 'INET006','INET010',...
    'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
    'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET062', 'INET063',...
    'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %'INET061', 'INET001', 
    sessions = 4;
elseif strcmpi(dataset, 'iNet-FSU')
    subject = {'INET074', 'INET075', 'INET077', 'INET078', 'INET083', 'INET084',...
    'INET085', 'INET086', 'INET087', 'INET088',...
    'INET091', 'INET093', 'INET094', 'INET096', 'INET098', 'INET099', 'INET101',...
    'INET103', 'INET104', 'INET105', 'INET106', 'INET107', 'INET108', 'INET109',...
    'INET110', 'INET111', 'INET112', 'INET114', 'INET115', 'INET116', 'INET120',...
    'INET123', 'INET133', 'INET136', 'INET137', 'INET140', 'INET141', 'INET143',...
    'INET156', 'INET158', 'INET159', 'INET160', 'INET165', 'INET168'};
    sessions = 4;
elseif strcmpi(dataset, 'Lifespan-FSU')
    subject = {'LS31', 'LS32', 'LS33', 'LS39', 'LS43', 'LS44', 'LS45', 'LS46',...
    'LS47', 'LS54', 'LS61', 'LS62', 'LS63', 'LS68', 'LS70', 'LS71', 'LS72',...
    'LS76', 'LS77', 'LS79', 'LS85', 'LS89', 'LS94', 'LS108'};
    sessions = 5;
end
subject(:,find(contains(subject, exclude_subs))) = [];
N = numel(subject);

end


