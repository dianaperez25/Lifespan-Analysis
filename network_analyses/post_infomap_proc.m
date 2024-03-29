%% post infomap proc
% regularizes and combines the output of infomap across different
% thresholds into one map assigning different weights to smaller vs. larger
% thresholds
clear all
%addpath(genpath('/projects/b1081/Scripts/graphtools/'))
addpath(genpath('/Users/dianaperez/Documents/Dependencies/graphtools/'))

%% CHANGE THIS
data_type = 'parcels'; % parcels or vertices
minNetSize = 20; %400 for vertices, 20 for parcels?
sub = 'LS05';
str = ['sub-' sub '_infomap'];
%cd(['/scratch/dcr8536/infomap/' sub '/parcels/'])
cd(['/Users/dianaperez/Desktop/Lifespan_infomap/infomap/' sub '/parcels/'])
%template_fname = '/scratch/dcr8536/template.dtseries.nii';
template_fname = '/Users/dianaperez/Desktop/template.dtseries.nii';

%% REGULARIZATION
simple = modify_clrfile('simplify','rawassn.txt',minNetSize); %makes a file called rawassn_minsizeX.txt
regularized = rawoutput2clr(simple);
regularized(regularized < 2) = 0;
regularized = regularized - 1;
dlmwrite(['rawassn_minsize' num2str(minNetSize) '_regularized.txt'],regularized,'delimiter',' ') %makes a file called rawassn_minizeX_regularized.txt

%% COMBINING ACROSS THRESHOLDS
thresholdarray = [0.0005:0.0005:0.004 0.005:0.005:0.05];%[0.003 0.004 0.005:0.005:0.05];
infoassn = dlmread(['rawassn_minsize' num2str(minNetSize) '_regularized.txt']);
%run conBensus
conBensus(infoassn, str, [], thresholdarray*100, 'parcel', minNetSize)
load([str '_conBensus_weighted_minsize' num2str(minNetSize) '.mat'])
%load template and save output of conBensus
tmp = ft_read_cifti_mod(template_fname);
switch data_type
    case 'parcels'
        % load parcellation file:
        parcels = ft_read_cifti_mod(['/Users/dianaperez/Desktop/Lifespan_parcellations/GordonParcellation/sub-' sub '/sub-' sub '_individual_parcels_edgethresh_0.5.dtseries.nii']);
        %parcels = ft_read_cifti_mod(['/scratch/dcr8536/parcellations/sub-' sub '/sub-' sub '_individual_parcels_edgethresh_0.5.dtseries.nii']);
        unique_parcels = unique(parcels.data);
        if unique_parcels(1) == 0 
            unique_parcels(1) = [];
        end
        
        for parc = 1:length(unique_parcels)
            inds = find(parcels.data==unique_parcels(parc));
            tmp.data(inds) = consen(parc);
        end
        ft_write_cifti_mod([str '_conBensus_weighted_minsize' num2str(minNetSize) '_parcels.dtseries.nii'], tmp);
    case 'vertices'
        tmp.data = consen;
        ft_write_cifti_mod([str '_conBensus_weighted_minsize' num2str(minNetSize) '_vertices.dtseries.nii'], tmp);
end

%re-label the networks so that they're all consecutive numbers
colorChange = [];
inds = unique(consen);
colorChange(:,1) = inds;
colorChange(:,2) = 0:length(inds)-1;
outmat = raw2colors_mat(consen,colorChange);
tmp.data = zeros(size(tmp.data));

switch data_type
    case 'parcels'        
        for parc = 1:length(unique_parcels)
            inds = find(parcels.data==unique_parcels(parc));
            tmp.data(inds) = outmat(parc);
        end
        ft_write_cifti_mod([str '_conBensus_weighted_minsize' num2str(minNetSize) '_consecutive_individual_parcels.dtseries.nii'], tmp);
    case 'vertices'
        tmp.data = outmat;
        ft_write_cifti_mod([str '_conBensus_weighted_minsize' num2str(minNetSize) '_consecutive_individual_vertices.dtseries.nii'], tmp);
end



