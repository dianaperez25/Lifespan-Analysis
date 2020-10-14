?addpath(genpath('/projects/b1081/Scripts'))

cortexOnly = 0;
if cortexOnly == 1
    voxnum = 59412; % Cortex + subcortex/cerebellum: 65625 %Cortex only: voxnum = 59412;
else
    voxnum = 65625;
end

data_folder = '/projects/b1081/Lifespan/derivatives/postFCproc_CIFTI/';
tmask_folder = '/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.0.6/fmriprep/';
template_fname = '/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-INET003_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii';
subject = 'LS03';
sessions = [2:5];
runs = [9,9,11,8,9];
output_file = sprintf('%s/dconn_cifti_normalwall/sub-%s_tmasked.dconn.nii', data_folder, subject);
catData = [];
catTmask = [];

%LOAD DATA
        
for i = 1:numel(sessions)
    for j = 1:runs(i)
        disp(sprintf('Loading data for session %d run %02d...', sessions(i), j))
        data = ft_read_cifti_mod(sprintf('%s/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%02d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',data_folder,subject,sessions(i),j));
        input_data = data.data;
        clear data;

        disp(sprintf('Loading tmask for session %d run %02d...', sessions(i), j))
        tmask = sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%02d_desc-tmask_fFD.txt',tmask_folder,subject,sessions(i),subject,sessions(i),j);
        tmask_data = table2array(readtable(tmask)); 
        masked_data = input_data(:,logical(tmask_data));
        clear tmask_data;
        
        disp('concatenating data')
        catData = [catData masked_data];        
        disp(sprintf('data size is now %d by %d', size(catData,1), size(catData,2)))
        clear masked_data;
    end
%end
%% Take correlations among vertices
% a single run takes about 1 min
disp(sprintf('Starting dconn correlations on data size %d by %d...', size(catData,1), size(catData,2)))
tic;
dconn_dat = paircorr_mod(catData');
disp('Finished correlations')
toc

% save out correlations (if desired)    
% load template:
template = ft_read_cifti_mod(template_fname);
%     
% % modify as needed
template.data = [];
template.dimord = 'pos_pos';
template.hdr.dim(7) = voxnum;
template.hdr.dim(6) = template.hdr.dim(7);
template.hdr.intent_code = 3001;
template.hdr.intent_name = 'ConnDense';
    
% substitute data
template.data = dconn_dat;
output_file = sprintf('/projects/p31161/sub-%s_ses-%d_tmasked.dconn.nii', subject, sessions(i));
% write dconn
% a single run takes about 30 seconds
tic;
disp('writing dconn...');
ft_write_cifti_mod(output_file,template);
disp('done');
toc
end
%dconn = ft_read_cifti_mod('/projects/b1081/Lifespan/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/sub-LS03_tmasked.dconn.nii');
%template = ft_read_cifti_mod('/projects/b1081/Lifespan/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/sub-LS03_tmasked.dconn.nii')
tic;
disp('making spatial correlations...')
%createSptlcorr_MSCdconns_timeseries('/projects/b1081/Atlases', '120_allsubs_corr',1, '/projects/p31161/', dconn_dat, 'spatialCorrMapLS03')
groupAvgLoc = '/projects/b1081/Atlases';
groupAvgName = '120_allsubs_corr';
outputdir = '/projects/p31161/';
outputname = 'spatialCorrMapLS03';

cifti_corrmap = dconn_dat(1:voxnum,1:voxnum);
cifti_corrmap(isnan(cifti_corrmap)) = 0;
cifti_corrmap1 = cifti_corrmap(1:voxnum,1:(voxnum/4));
disp(['First quarter of Cifti corrmap is ' num2str(size(cifti_corrmap1,1)) ' by ' num2str(size(cifti_corrmap1,2)) ': ' datestr(now)])
% Remove NaNs (produced if vertices have no data)


for dconnrow = 1:size(cifti_corrmap1,1)            
    cifti_corrmap1(dconnrow,:) = single(FisherTransform(cifti_corrmap1(dconnrow,:)));
end

disp(sprintf('First quarter of Fisher Transform Finished: %s', datestr(now)));
                    
disp(sprintf('Loading Template: %s', datestr(now)));

% Load group-average corrmat (assumes it's a dconn)
group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
group = group.data(1:voxnum,1:voxnum);
group1 = group(1:voxnum,1:(voxnum/4));
sizedata = size(group1,2);

if cortexOnly==1
    % Load template variants file
    template = ft_read_cifti_mod('/projects/b1081/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
    template.data = [];
else
    template = ft_read_cifti_mod('/projects/b1081/MSC/TaskFC/FCProc_MSC05_motor_pass2/cifti_timeseries_normalwall_native_freesurf/vc39006_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');
    template.data = [];
    template.hdr(6).dim = 1;
end

disp(['First quarter of Cifti group data is ' num2str(size(group1,1)) ' by ' num2str(size(group1,2)) ': ' datestr(now)])

disp(sprintf('Template Loaded: %s', datestr(now)));
            
            % Compare to group average
for i=1:sizedata
    template.data(i,1) = paircorr_mod(group1(:,i),cifti_corrmap1(:,i));
end
            
clear cifti_corrmap1
clear group1

disp(sprintf('First quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap2 = cifti_corrmap(1:voxnum,(voxnum/4)+1:(voxnum/2));

disp(['Second quarter of Cifti corrmap is ' num2str(size(cifti_corrmap2,1)) ' by ' num2str(size(cifti_corrmap2,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap2,1)
    cifti_corrmap2(dconnrow,:) = single(FisherTransform(cifti_corrmap2(dconnrow,:)));
end

disp(sprintf('Second quarter of Fisher Transform Finished: %s', datestr(now)));

group2 = group(1:voxnum,(voxnum/4)+1:(voxnum/2));
sizedata = size(group2,2);

disp(['Second quarter of Cifti group data is ' num2str(size(group2,1)) ' by ' num2str(size(group2,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data((i+(voxnum/4)),1) = paircorr_mod(group2(1:voxnum,i),cifti_corrmap2(1:voxnum,i));
end

clear cifti_corrmap2
clear group2

disp(sprintf('Second quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap3 = cifti_corrmap(:,(voxnum/2)+1:(voxnum*.75));

disp(['Third quarter of Cifti corrmap is ' num2str(size(cifti_corrmap3,1)) ' by ' num2str(size(cifti_corrmap3,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap3,1)
    cifti_corrmap3(dconnrow,:) = single(FisherTransform(cifti_corrmap3(dconnrow,:)));
end

disp(sprintf('Third quarter of Fisher Transform Finished: %s', datestr(now)));

group3 = group(:,(voxnum/2)+1:(voxnum*.75));
sizedata = size(group3,2);

disp(['Third quarter of Cifti group data is ' num2str(size(group3,1)) ' by ' num2str(size(group3,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum/2),1) = paircorr_mod(group3(:,i),cifti_corrmap3(:,i));
end

clear cifti_corrmap3
clear group3

disp(sprintf('Third quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap4 = cifti_corrmap(:,(voxnum*.75)+1:end);

disp(['Fourth quarter of Cifti corrmap is ' num2str(size(cifti_corrmap4,1)) ' by ' num2str(size(cifti_corrmap4,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap4,1)
    cifti_corrmap4(dconnrow,:) = single(FisherTransform(cifti_corrmap4(dconnrow,:)));
end

disp(sprintf('Fourth quarter of Fisher Transform Finished: %s', datestr(now)));

group4 = group(:,(voxnum*.75)+1:end);
sizedata = size(group4,2);

disp(['Fourth quarter of Cifti group data is ' num2str(size(group4,1)) ' by ' num2str(size(group4,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum*.75),1) = paircorr_mod(group4(:,i),cifti_corrmap4(:,i));
end

clear cifti_corrmap4
clear group4

disp(sprintf('Fourth quarter of Correlation Finished: %s', datestr(now)));


disp(['Final size of spatial correlation map is ' num2str(size(template.data,1)) ' by ' num2str(size(template.data,2)) ': ' datestr(now)])


% Write out the variants

if cortexOnly == 1

    ft_write_cifti_mod([outputdir '/' outputname '_vs_' groupAvgName '_cortex_corr'],template)
    template.data = [];

elseif cortexOnly == 0

    ft_write_cifti_mod([outputdir '/' outputname '_vs_' groupAvgName '_subcort_corr'],template)
    template.data = [];

end
 


disp('done')
toc