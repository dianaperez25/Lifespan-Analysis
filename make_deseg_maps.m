%% Segregation Index Brain Maps
%This script creates a brain map of segregation index values (and also of 
%average within- and between-system correlations. 
%The inputs are a network assignment map (so that we know which network
%each vertex is assigned to) and a dconn(or vertex-wise correlation matrix;
%so that we know the FC between each vertex and every other vertex). The
%script is written to handle three ways of dealing with negative
%correlations, which is to be specified under the variable neg_corrs:
%1. 'nan' - this will delete all negative correlations when calculating
%   average within- and between-system correlations, and therefore when
%   calculating the segregation index
%2. 'zero' - this will set all negative correlations to zero 
%3. 'asis' - this will leave and consider all the negative correlations in the
%   calculations to be performed
%The end result is a cifti with a segregation index value per vertex that
%in theory should tell us whether there are regional differences in
%desegregation.

clear all
%% PATHS
root_dir = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/';
netmap_dir = [root_dir '/template_matching/'];
dconn_dir = [root_dir '/dconn_cifti_normalwall/Lifespan_Dconns/'];
out_dir = '/scratch/dcr8536/deseg_maps/';
if ~exist(out_dir)
    mkdir(out_dir);
end

%% VARIABLES
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
neg_corrs = 'zero'; %nan, zero, or asis

for sub = 1:numel(subject)
    disp(['working on subject ' subject{sub} '...'])
% STEP 1: load network map
    netmap_fname = [netmap_dir '/sub-' subject{sub} '_dice_WTA_map_kden0.05.dtseries.nii'];
    netmap = ft_read_cifti_mod(netmap_fname);
    clear netmap_fname
% STEP 2: load dconn (aka our correlation matrix)
    dconn_fname = [dconn_dir '/sub-' subject{sub} '_allsess_tmasked.dconn.nii'];
    dconn = ft_read_cifti_mod(dconn_fname);
    clear dconn_fname
% STEP 3: create a desegregation map to be filled in with vertex-wise values   
    deseg_map= zeros(size(netmap.data));
    within_corr_map = zeros(size(netmap.data));
    between_corr_map = zeros(size(netmap.data));
% STEP 3: calculate segregation index at each vertex
    for vertex = 1:size(netmap.data,1)
        % first, we have to get the within-system correlations
        network_assigned = netmap.data(vertex); %this is the network assigned to this vertex
        if network_assigned > 0 % if this vertex is assigned a network...
            % we are finding all the vertices that are also assigned to
            % the same network
            within_nodes = find(netmap.data==network_assigned);
            % but also let's make sure that we are not counting this
            % vertex's correlation with itself
            delete_vertex = find(within_nodes==vertex);
            within_nodes(delete_vertex) = [];     
            clear delete_vertex
            within_corrs = dconn.data(vertex,within_nodes);
            if strcmpi(neg_corrs, 'nan')
                within_corrs(within_corrs<0) = [];
            elseif strcmpi(neg_corrs, 'zero')
                within_corrs(within_corrs<0) = 0;
            end
            average_within = mean(within_corrs);
            
            % now we find between-network correlations
            % find the vertices NOT assigned to the same network
            between_nodes = find(netmap.data~=network_assigned);
            clear network_assigned
            between_corrs = dconn.data(vertex,between_nodes);
            if strcmpi(neg_corrs, 'nan')
                between_corrs(between_corrs<0) = [];
            elseif strcmpi(neg_corrs, 'zero')
                between_corrs(between_corrs<0) = 0;
            end
            average_between = mean(between_corrs);
            
            % now we calculate the segregation index
            deseg_map(vertex) = (average_within - average_between)/average_within;
            within_corr_map(vertex) = average_within;
            between_corr_map(vertex) = average_between;
        end           
    end
template = netmap;
template.data = deseg_map;
ft_write_cifti_mod([out_dir '/sub-' subject{sub} '_desegregation_brain_map_negcorrs_' neg_corrs '.dtseries.nii'], template);
template.data = within_corr_map;
ft_write_cifti_mod([out_dir '/sub-' subject{sub} '_within_correlations_brain_map_negcorrs_' neg_corrs '.dtseries.nii'], template);
template.data = between_corr_map;
ft_write_cifti_mod([out_dir '/sub-' subject{sub} '_between_correlations_brain_map.dtseries_' neg_corrs '.nii'], template);
end