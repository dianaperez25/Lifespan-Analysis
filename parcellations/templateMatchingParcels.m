function templateMatchingVariants(paramsFile,outputdir)
subject = 'LS02';
sessions = 5;
runs = 13; % I think this is the max number of runs
% some paths
cifti_data_dir = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/';
tmask_dir = '/projects/b1081/Lifespan/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
ind_parcels_dir = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/indiv_parcels/';

% load network templates
load('/projects/b1081/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates
count = 0;
for i=1:length(IDNames)
    if strcmp(IDNames{i},'skip')
        continue;
    else
        count = count + 1;
        tempIDNames{count} = IDNames{i};
    end
end
IDNames = tempIDNames;

templates = templates(1:59412,:)';
template_values_sorted = sort(templates(:), 'descend');
threshval= template_values_sorted(round(numel(template_values_sorted) .* 0.05));
threshtemplates= templates >= threshval;
clear allTvals allTvals_sorted templates threshval

% Set variables
if ~exist('outputdir')
    outputdir = pwd;
end

% load individual parcellation file for subject
ind_parcels_fname = sprintf('%s/sub-%s/sub-%s_individual_parcels_edgethresh_0.5.dtseries.nii', ind_parcels_dir, subject, subject);
ind_parcels = ft_read_cifti_mod(ind_parcels_fname);
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
            cifti_ts_concat = [cifti_ts_concat cifti_ts];
            tmask_concat = [tmask_concat tmask];
        end 
    end
end

masked_data = cifti_ts_concat(:,logical(tmask_concat'));
unique_parcels = unique(ind_parcels.data);
if unique_parcels(1) == 0;
    unique_parcels(1) = [];
end
    % Loop through each variant and calculate its connectivity map
    for ind=1:length(unique_parcels)
        parcel_verts = logical(ind_parcels.data==unique_parcels(ind));
        corr_map = paircorr_mod(mean(masked_data(parcel_verts,:))',masked_data'); 
        corr_map(isnan(corr_map)) = 0;
        corr_map_parcel(:,ind) = corr_map';
        clear corr_map
    end
    clear masked_data
    
    
    % Match each variant's connectivity map to each template connectivity map
    disp('Matching to templates and computing goodness of fit')
    corr_coeff = zeros(length(unique_parcels),size(threshtemplates,2));
    %% threshtemplates here is supposed to have the number of templates
    for parcel=1:length(unique_parcels)
        disp(['Subject ' subject ', parcel #' num2str(parcel) ' out of ' num2str(length(unique_parcels))]);
        for templatenum = 1:size(threshtemplates,2); %% this isn't doing dice correlation, is it? It's just straight up correlation....
            corr_coeff(parcel,templatenum) = paircorr_mod(corr_map_parcel(:,parcel),threshtemplates(:,templatenum));
        end
    end
    
    
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
    for parcel=1:length(unique_parcels)
        parcel_verts = logical(ind_parcels.data==unique_parcels(parcel));
        corrVals = corr_coeff(parcel,:);
        networkIDs(parcel_verts,count) = IDs(maxi(parcel));
        networkRatio(parcel_verts,count) = corr_coeff(parcel,maxi(parcel))./corr_coeff(parcel,nextMax(parcel));
        for kk=1:size(ThreshTemplates,2)
            networkPercent(parcel_verts,kk,count)=corrVals(kk);
        end
    end
    clear maxi corr_coeff nextMax corrVals
end


% Write out the results
out_template.data = zeros(ncortverts,count);

out_template.data(1:size(networkIDs,1),:) = networkIDs;
ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_networkIDs'],out_template);
clear networkIDs

out_template.data(1:size(networkRatio,1),:) = networkRatio;
ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_ratioOfTopTwoTemplates'],out_template);
clear networkRatio

for kk = 1:size(ThreshTemplates,2)
    out_template.data(1:size(networkPercent,1),:) = squeeze(networkPercent(:,kk,:));
    ft_write_cifti_mod([outputdir '/variantTemplatematch_allSubjects_corrCoeffFor' IDNames{kk}],out_template);
end

end

function [pathstr, name, ext] = fileparts(file)
%FILEPARTS Filename parts.
%   [PATHSTR,NAME,EXT] = FILEPARTS(FILE) returns the path, file name, and
%   file name extension for the specified FILE. The FILE input is a string
%   containing the name of a file or folder, and can include a path and
%   file name extension. The function interprets all characters following
%   the right-most path delimiter as a file name plus extension.
%
%   If the FILE input consists of a folder name only, be sure that the
%   right-most character is a path delimiter (/ or \). Othewise, FILEPARTS
%   parses the trailing portion of FILE as the name of a file and returns
%   it in NAME instead of in PATHSTR.
%
%   FILEPARTS only parses file names. It does not verify that the file or
%   folder exists. You can reconstruct the file from the parts using
%      fullfile(pathstr,[name ext])
%
%   FILEPARTS is platform dependent.
%
%   On Microsoft Windows systems, you can use either forward (/) or back
%   (\) slashes as path delimiters, even within the same string. On Unix
%   and Macintosh systems, use only / as a delimiter.
%
%   See also FULLFILE, PATHSEP, FILESEP.

%   Copyright 1984-2012 The MathWorks, Inc.
%   $Revision: 1.18.4.18 $ $Date: 2012/04/14 04:15:41 $

pathstr = '';
name = '';
ext = '';

if ~ischar(file)
    error(message('MATLAB:fileparts:MustBeChar'));
elseif isempty(file) % isrow('') returns false, do this check first
    return;
elseif ~isrow(file)
    error(message('MATLAB:fileparts:MustBeChar'));
end

if ispc
    ind = find(file == '/'|file == '\', 1, 'last');
    if isempty(ind)
        ind = find(file == ':', 1, 'last');
        if ~isempty(ind)       
            pathstr = file(1:ind);
        end
    else
        if ind == 2 && (file(1) == '\' || file(1) == '/')
            %special case for UNC server
            pathstr =  file;
            ind = length(file);
        else 
            pathstr = file(1:ind-1);
        end
    end
    if isempty(ind)       
        name = file;
    else
        if ~isempty(pathstr) && pathstr(end)==':' && ...
                (length(pathstr)>2 || (length(file) >=3 && file(3) == '\'))
                %don't append to D: like which is volume path on windows
            pathstr = [pathstr '\'];
        elseif isempty(deblank(pathstr))
            pathstr = '\';
        end
        name = file(ind+1:end);
    end
else    % UNIX
    ind = find(file == '/', 1, 'last');
    if isempty(ind)
        name = file;
    else
        pathstr = file(1:ind-1); 

        % Do not forget to add filesep when in the root filesystem
        if isempty(deblank(pathstr))
            pathstr = '/';
        end
        name = file(ind+1:end);
    end
end

if isempty(name)
    return;
end

% Look for EXTENSION part
ind = find(name == '.', 1, 'last');

if isempty(ind)
    return;
else
    ext = name(ind:end);
    name(ind:end) = [];
end

end