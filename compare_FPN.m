subjects = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17',...
    'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'};

% load surface area file
surf_area = ft_read_cifti_mod('/scratch/dcr8536/surf_areas_verts.dtseries.nii');
% load group average
group_avg = ft_read_cifti_mod('/scratch/dcr8536/WashU120_groupNetworks.dtseries.nii');
group_FP = zeros(size(group_avg.data)); group_FP(group_avg.data==3) = 1;
group_FP = [group_FP(group_avg.brainstructure==1); group_FP(group_avg.brainstructure==2)];
group_DMN = zeros(size(group_avg.data)); group_DMN(group_avg.data==1) = 1;
group_DMN = [group_DMN(group_avg.brainstructure==1); group_DMN(group_avg.brainstructure==2)];
clear group_avg
%initialize variables to store information
FPN = [];
DMN = [];
for sub = 1:numel(subjects)
    %load network map
    netmap_fname = sprintf('/scratch/dcr8536/postFCproc_CIFTI/template_matching/sub-%s_dice_WTA_map_kden0.05.dtseries.nii', subjects{sub});
    net_map = ft_read_cifti_mod(netmap_fname);
    %index FPN vertices in individual
    FPN_verts = find(net_map.data==3);
    indiv_FP = zeros(59412,1);
    indiv_FP(FPN_verts) = 1;
    %index DMN vertices
    DMN_verts = find(net_map.data==1);
    indiv_DMN = zeros(59412,1);
    indiv_DMN(DMN_verts) = 1;
    %calculate surface area of FPN
    FPN(sub,1) = sum(surf_area.data(FPN_verts));
    %calculate surface area of DMN
    DMN(sub,1) = sum(surf_area.data(DMN_verts));
    %calculate dice overlap bewteen individual FPN and group average FPN
    FPN(sub,2) = dice_coefficient_mod(group_FP, indiv_FP);
    %calculate dice overlap bewteen individual DMN and group average DMN
    DMN(sub,2) = dice_coefficient_mod(group_DMN, indiv_DMN);
    dconn_fname = sprintf('/scratch/dcr8536/iNetworks/Nifti/postFCproc_CIFTI/dconn_cifti_normalwall/sub-%s_allsess_tmasked.dconn.nii', subjects{sub});
    if exist(dconn_fname)
        %load subject's dconn
        dconn = ft_read_cifti_mod(dconn_fname);
        %index FPN functional connectivity
        iFPN_fc = dconn.data(FPN_verts,:); %individual FPN
        gFPN_fc = dconn.data(group_FP==1,:); %group FPN
        %index DMN functional connectivity
        iDMN_fc = dconn.data(DMN_verts,:); %individual DMN
        gDMN_fc = dconn.data(group_DMN==1,:); %group DMN
        clear dconn
        %calculate FPN within-network connectivity
        FPN(sub,3) = mean(iFPN_fc(:,FPN_verts)>0, 'all', 'omitnan');
        %calculate FPN within-network connectivity
        DMN(sub,3) = mean(iDMN_fc(:,DMN_verts)>0, 'all', 'omitnan');
        %calculate FPN between-network connectivity
        FPN_mask = ones(size(iFPN_fc));
        FPN_mask(:,FPN_verts) = 0;
        FPN(sub,4) = mean(iFPN_fc(logical(FPN_mask)));
        clear FPN_mask
        %calculate DMN between-network connectivity
        DMN_mask = ones(size(iDMN_fc));
        DMN_mask(:,DMN_verts) = 0;
        DMN(sub,4) = mean(iDMN_fc(logical(DMN_mask)));
        clear DMN_mask iDMN_fc 
        %calculate FP-DMN connectivity
        FPN(sub,7) = mean(iFPN_fc(:,DMN_verts)>0,'all','omitnan');
        clear iFPN_fc 
        %calculate FPN within-network connectivity
        FPN(sub,5) = mean(gFPN_fc(:,group_FP==1)>0, 'all', 'omitnan');
        %calculate FPN within-network connectivity
        DMN(sub,5) = mean(gDMN_fc(:,group_DMN==1)>0, 'all', 'omitnan');
         %calculate FPN between-network connectivity
        FPN_mask = ones(size(gFPN_fc));
        FPN_mask(:,group_FP==1) = 0;
        FPN(sub,6) = mean(gFPN_fc(logical(FPN_mask)));
        clear FPN_mask
        %calculate DMN between-network connectivity
        DMN_mask = ones(size(gDMN_fc));
        DMN_mask(:,group_DMN==1) = 0;
        DMN(sub,6) = mean(gDMN_fc(logical(DMN_mask)));
        clear DMN_mask gDMN_fc 
        %calculate FP-DMN connectivity
        FPN(sub,8) = mean(gFPN_fc(:,group_DMN==1)>0,'all','omitnan');
        clear gFPN_fc 
    end
    results = [FPN DMN];
end