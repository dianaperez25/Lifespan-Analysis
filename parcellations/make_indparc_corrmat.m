clear all

root_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/parcellations/';
parcel_timecourse_dir = [root_dir 'Gordon_Parcellation/avg_timecourses/'];
network_map_dir = [root_dir 'Gordon_Parcellation/'];
out_dir = '/Users/dianaperez/Desktop/ohbm_poster/';

subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17',...
    'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
    'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
    'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
    'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
    'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
    'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
    'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; 
match_data = 1;
amt_data = 1021;

for sub = 1:numel(subjects)
        % load cifti with individual parcels
        parcel_fname = sprintf('%s/sub-%s/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', network_map_dir, subjects{sub}, subjects{sub});
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
        end       

        %load network map to obtain network assignments
        netmap_fname = sprintf('%s/sub-%s/sub-%s_indiv_parcels_net_assigned_binarized.dtseries.nii', network_map_dir, subjects{sub}, subjects{sub});
        netmap = ft_read_cifti_mod(netmap_fname);

         % get unique network IDs
        unique_net_IDs = unique(netmap.data);
        if unique_net_IDs(1) == 0
            unique_net_IDs(1) = []; %but delete the 0
        end
        
        % now iterate across networks...
        net_size = [];
        for net = 1:length(unique_net_IDs)
            % ...to find the parcels that are assigned to each network...
            net_size(net,1) = unique_net_IDs(net); 
            net_verts = find(netmap.data==unique_net_IDs(net)); %...now get indices for the vertices assigned to each network...
            net_parcels = unique(consec_parc_map(net_verts)); %...and get the parcel ID's that correspond to those vertices...
            net_size(net,2) = length(net_parcels); %...get the size of the network in number of parcels...
        end
        
        data_concat = [];
        for ses = 1:5
            %load parcels average timecourse; these have already been masked 
            parcel_tc_fname = sprintf('%s/sub-%s_ses-%d_indiv_parcels_timseries_sorted.mat', parcel_timecourse_dir, subjects{sub},ses);
            if exist(parcel_tc_fname)
                mat_struct = load(parcel_tc_fname);
            else continue;
            end

            parc_tc_sorted = mat_struct.ses_tseries;
            
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

        make_corrmat(matrix_sorted, net_size)
        fout_str = sprintf('%s/corrmats/sub_%s_indparcs_corrmat.tiff',out_dir,subjects{sub});
        saveas(gcf,[fout_str '.tiff'],'tiff');
        close all

end


function make_corrmat(matrix, net_size)
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

lines = [];
count = 0;
for net = 1:size(net_size,1)
    tickpos(net) = count + (.5*net_size(net,2));
    count = count+net_size(net,2);
    lines = [lines; count];
end
% put lines between the networks
vline_new([lines]+.5,'k',3);
hline_new([lines]+.5,'k',3);
ax = axis;

% put ticks in the right places
set(gca,'XTick',tickpos,'Xlim',[ax(1) ax(2)]);

% take out the current tick labels
set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

% label the networks on the two axes
tx= text(tickpos,ones(1,length(tickpos))*(lines(end)+1),networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
set(tx,'Color',[0,0,0],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
set(gca,'FontWeight','bold','FontSize',10);

ty= text(-10*ones(1,length(tickpos)),tickpos,networks);
set(ty,'Color',[0,0,0],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

set(gca,'FontWeight','bold','FontSize',10);

% make square
axis square;
colorbar;
end
