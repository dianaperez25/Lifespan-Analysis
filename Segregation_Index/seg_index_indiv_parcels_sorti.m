%% SEGREGATION INDEX CALCULATION
% This script calculates a measure of system segregation for individualized parcels.

% Inputs: 
    %   individualized parcellation, 
    %   individualized parcels timeseries 
    %   parcel network assignment
% Outputs: 
    %   mean within-system correlation 
    %   mean between-system correlation
    %   segregation index per network
    %   segregation index per subject
    %   (detailed explanation of outputs at the end of the script)
    
% Based on the analysis describe in: 
% Chan, Micaela Y., et al. "Decreased segregation of brain systems across 
% the healthy adult lifespan." Proceedings of the National Academy of 
% Sciences 111.46 (2014): E4997-E5006.
%
% NOTES: in original analysis negative z-values were set to zero.
% Within-system connectivity was calculated as the mean node-to-node 
% z-value of all nodes of that system to each other. 
% Between-system connectivity was calculated as the mean node-to-node 
% z-value between each node of a system and all nodes of all other systems. 

% 11/28/23 - DP: adding some code to sort ses-specific timecourses for inet
% subs
% ------------------------------------------------------------------------

clear all

% ------------------------------------------------------------------------
%% PATHS
% ------------------------------------------------------------------------
root_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/parcellations/';
parcel_timecourse_dir = [root_dir 'Gordon_Parcellation/avg_timecourses/'];
network_map_dir = [root_dir 'Gordon_Parcellation/'];
atlas_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Atlases/';
out_dir = '/Users/dianaperez/Desktop/ohbm_poster/';

if ~exist(out_dir)
    mkdir(out_dir);
end

% ------------------------------------------------------------------------
%% OPTIONS
% ------------------------------------------------------------------------
match_data = 1; % if 1, will calculate the minimum possible amount of data available and will force all subs to have that amount of data
amt_data = 1021; % if this is commented out or set to 0, then the script will calculate it
atlas = 'Seitzman300';
datasets = {'lifespan', 'inetworks'}; %'lifespan', 'inetworks'

% ------------------------------------------------------------------------
%% VARIABLES
% ------------------------------------------------------------------------
% load subject id's
subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17',...
'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %'INET001', 'INET002', 

% load atlas information
atlas_params = atlas_parameters_GrattonLab(atlas, atlas_dir);%[root_dir 'Atlases']);%
network_z.networks = atlas_params.networks(2:end); %copy names of each network

% ------------------------------------------------------------------------
%% DATA MATCHING
% ------------------------------------------------------------------------
if match_data
    if amt_data == 0 || ~exist('amt_data')
        allSubs_amtData = [];
        % get minimum amt of data per session
        allSubs_amtData = [];
        for sub = 1:numel(subjects)            
            for ses = 1:10
                %if contains(subjects{sub}, 'LS')
                    fname = [parcel_timecourse_dir '/sub-' subjects{sub} '_ses-' num2str(ses) '_indiv_parcels_timseries_sorted.mat'];
                    if exist(fname)
                        load(fname)
                        allSubs_amtData(sub,ses) = size(ses_tseries,2); clear ses_tseries
                    else continue;
                    end
                %elseif contains(subjects{sub}, 'INET')
                    %fname = [parcel_timecourse_dir '/sub-' subjects{sub} '_ses-' num2str(ses) '_individual_parcels_timecourse.mat'];
                    %if exist(fname)
                        %load(fname)
                        %allSubs_amtData(sub,ses) = size(avg_parc_ts,2); clear avg_parc_ts
                    %end
                %end                
            end
        end
        amt_data = min(min(allSubs_amtData));
    end
end

% ------------------------------------------------------------------------
%% MAIN FOR-LOOP
% ------------------------------------------------------------------------

for group = 1:numel(datasets)
    
    % initializing variables
    corrmats = []; group_struct = {};
    
    if strcmpi(datasets{group}, 'lifespan')
        inds = find(contains(subjects, 'LS'));
        subs = subjects(inds);
        sort = 0;
    elseif strcmpi(datasets{group}, 'inetworks')
        inds = find(contains(subjects, 'INET'));
        subs = subjects(inds);
        sort = 0;
    end
    
    for sub = 1:numel(subs)
    
        % initialize a couple variables
        net_size = []; parcels_sorted = []; sub_struct = {}; sub_within_corrs = []; sub_between_corrs = [];

        % load cifti with individual parcels
        parcel_fname = sprintf('%s/sub-%s/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', network_map_dir, subs{sub}, subs{sub});
        parcel_map = ft_read_cifti_mod(parcel_fname);

        % get unique parcel IDs
        unique_IDs = unique(parcel_map.data);
        if unique_IDs(1)==0
            unique_IDs(1)=[]; %delete parcels equal to 0
        end

        % parcel IDs are not consecutive; making consecutive here 
        consec_parc_map = zeros(size(parcel_map.data));    
        for parc = 1:length(unique_IDs)
            consec_parc_map(find(parcel_map.data==unique_IDs(parc))) = parc;
        end        %load network map to obtain network assignments
        netmap_fname = sprintf('%s/sub-%s/sub-%s_indiv_parcels_net_assigned_binarized.dtseries.nii', network_map_dir, subs{sub}, subs{sub});
        netmap = ft_read_cifti_mod(netmap_fname);

        % get unique network IDs
        unique_net_IDs = unique(netmap.data);
        if unique_net_IDs(1) == 0
            unique_net_IDs(1) = []; %but delete the 0
        end
        
            % now iterate across networks...
        for net = 1:length(unique_net_IDs)
            % ...to find the parcels that are assigned to each network...
            net_size(net,1) = unique_net_IDs(net); 
            net_verts = find(netmap.data==unique_net_IDs(net)); %...now get indices for the vertices assigned to each network...
            net_parcels = unique(consec_parc_map(net_verts)); %...and get the parcel ID's that correspond to those vertices...
            net_size(net,2) = length(net_parcels); %...get the size of the network in number of parcels...
            net_size(net,3) = length(net_verts); % ...and also get the size of each network in number of vertices...
            parcels_sorted = [parcels_sorted; net_parcels]; %...now add the parcel IDs for this network to a structure...
        end
            
        total_num_parcels = length(parcels_sorted); %...the total number of parcels for this subject...
        
        data_concat = [];
        for ses = 1:5
            %load parcels average timecourse; these have already been masked 
            parcel_tc_fname = sprintf('%s/sub-%s_ses-%d_indiv_parcels_timseries_sorted.mat', parcel_timecourse_dir, subs{sub},ses);
            if exist(parcel_tc_fname)
                mat_struct = load(parcel_tc_fname);
            else continue;
            end
            if sort == 1
                parc_tc_sorted = mat_struct.ses_tseries(parcels_sorted,:); %...now sorting the timeseries information for the parcels by network...
            else
                parc_tc_sorted = mat_struct.ses_tseries;
            end
            clear mat_struct
            % ... match the amount of data ...
            if match_data == 0 %if we don't care about matching data, then use the max amount of data available per subject/session
                amt_data = size(parc_tc_sorted,2);
            end
            if size(parc_tc_sorted,2)<amt_data
                continue;
            else
                matched_data = parc_tc_sorted(:,1:amt_data);
                data_concat = [data_concat matched_data];
            end
        end

        % ... calculate person correlation ...
        matrix_sorted = single(FisherTransform(paircorr_mod(data_concat')));% fisher transform r values

        count = 1; % start a count

        for net = 1:size(net_size,1) % iterate through each network
            % ...extract the parcels belonging to system net...
            net_parcels_all{net} = count:(count+net_size(net,2)-1); % put them all in a structure so that we can use them later             
            net_parcels = net_parcels_all{net}; % putting the cell array in a matrix so that we can use them now
            weights(sub,net) = net_size(net,2)/length(parcels_sorted); % calculate weights based on the size of the network

            %% ...GET WITHIN-SYSTEM CORRELATIONS...
            % a temporary matrix with correlation values within system net
            tmp = matrix_sorted(net_parcels(1):net_parcels(end),net_parcels(1):net_parcels(end)); 

            %make a mask to only keep correlation values in upper triangle
            maskmat = ones(net_size(net,2)); 
            maskmat = logical(triu(maskmat,1));
            within = tmp(maskmat); % mask out values in lower triangle
            within(within<0) = 0; % set negative z values to 0  

            % put those z values in a structure to use later
            withinFC_by_net(sub,net) = mean(within); 
            sub_within_corrs = [sub_within_corrs; within];
            clear maskmat tmp % clear up some variables

            %% ...GET BETWEEN SYSTEM CORRELATIONS...
            % a temporary matrix with z values for network net
            tmp = matrix_sorted(net_parcels(1):net_parcels(end), 1:total_num_parcels);
            % make a mask the same size of tmp
            maskmat = ones(size(tmp)); 
            maskmat(:,net_parcels(1):net_parcels(end)) = 0; % mask out z values within system net
            between = tmp(maskmat==1); % put only between system z values in a variable
            between(between<0) = 0;% set negative z values to 0

            % put those z values in a structure to use later
            sub_between_corrs = [sub_between_corrs; between];
            betweenFC_by_net(sub,net) = mean(between); 

            clear maskmat between within

            %% ...CALCULATE SEGREGATION FOR THIS NETWORK...
            segregation_by_net(sub,net) = (withinFC_by_net(sub,net) - betweenFC_by_net(sub,net))/withinFC_by_net(sub,net);
            count = count + net_size(net,2);
        end
        
        for net = 1:size(net_size,1) % iterate through each network again
            % get correlation with network net and every other network
            net_parcels = net_parcels_all{net}; 
            tmp = matrix_sorted(net_parcels(1):net_parcels(end), 1:total_num_parcels);
            for net2 = 1:size(net_size,1)
                net2_parcels = net_parcels_all{net2};            
                maskmat = zeros(size(tmp)); %create another mask the size of tmp
                maskmat(:,net2_parcels(1):net2_parcels(end)) = 1; % extract z values between net and net2
                FC = tmp(maskmat==1); 
                sub_corrmat(net,net2) = mean(FC); % put these in a matrix
            end
        end
        corrmats(sub,:,:,:) = sub_corrmat; % put the sub corrmat in a structure with all other subjects
        %% create an image with those correlations
        %% STILL NEED TO TEST THIS FUNCTION
        % add line to save the file
        outfile_fig = sprintf('%s/sub-%s_corrmat_indiv_parcels', out_dir, subs{sub});
        make_indparc_corrmat(sub_corrmat)
        saveas(gcf,[outfile_fig '.tiff'],'tiff');
        close(gcf);
        %% calculate segregation index across all networks 
        %segregation index = (mean within system Z - mean between system Z)/mean within system Z
        segregation_by_sub(sub,1) = (mean(sub_within_corrs) - mean(sub_between_corrs))/mean(sub_within_corrs);
        segregation_by_sub(sub,2) = mean(sub_within_corrs);
        segregation_by_sub(sub,3) = mean(sub_between_corrs);

        sub_struct.averages = segregation_by_sub(sub,:);
        sub_struct.corrmat = sub_corrmat;
        sub_struct.networks(:,1) = segregation_by_net(sub,:)';
        sub_struct.networks(:,2) = withinFC_by_net(sub,:)';
        sub_struct.networks(:,3) = betweenFC_by_net(sub,:)';
        sub_struct.networks(:,4) = weights(sub,:)';
        sub_struct.net_size = net_size;

        clear sub_within_corrs sub_between_corrs sub_corrmat

        outfile_sub = sprintf('%s/sub-%s_segregation_index_indiv_parcels.mat', out_dir, subs{sub});
        save(outfile_sub, 'sub_struct');
        count = count + net_size(net,2);
    end
    group_struct.averages = segregation_by_sub;
    group_struct.networks.withinFC = withinFC_by_net;
    group_struct.networks.betweenFC = betweenFC_by_net;
    group_struct.networks.segregation = segregation_by_net;
    group_struct.corrmats.subjects = corrmats;
    group_struct.corrmats.average = squeeze(mean(corrmats(:,:,:,:),1));
    outfile_group = sprintf('%s/group_%s_segregation_index_indiv_parcels.mat', out_dir, datasets{group});
    save(outfile_group, 'group_struct');
    
    outfile_fig = sprintf('%s/group_%s_corrmat_indiv_parcels', out_dir, datasets{group});
    make_indparc_corrmat(squeeze(mean(corrmats(:,:,:,:),1)))
    saveas(gcf,[outfile_fig '.tiff'],'tiff');
    close(gcf);
end

function make_indparc_corrmat(matrix)

% information about networks and colors
networks = {'DMN','Vis','FPN', 'DAN','VAN','Sal','CON','SMd','SMl','Aud','Tpole','MTL','PMN','PON'};

% open figure
h = figure('Color',[0.8275 0.8275 0.8275],'Position',[56 143 1295 807]); %[56 143 1095 807]

% plot out matrix
climlow = -1;
climhigh = 1;
imagesc(matrix,[climlow climhigh]);

% My favorite colormap - edit if you want a different one
load better_jet_colormap.mat; % assume this is in the same folder
colormap(better_jet_colormap_diff);

% put lines between the networks
vline_new([1:14]+.5,'k',3);
hline_new([1:14]+.5,'k',3);
tickpos = 1:14;
ax = axis;

% put ticks in the right places
set(gca,'XTick',tickpos,'Xlim',[ax(1) ax(2)]);

% take out the current tick labels
set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

% label the networks on the two axes
tx= text(tickpos,ones(1,length(tickpos))*(14+1),networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
set(tx,'Color',[0,0,0],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
set(gca,'FontWeight','bold','FontSize',10);

ty= text(-.05*ones(1,length(tickpos)),tickpos,networks);
set(ty,'Color',[0,0,0],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

colorbar;
set(gca,'FontWeight','bold','FontSize',10);

% make square
axis square;
end
