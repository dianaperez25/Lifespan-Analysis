clear all
addpath(genpath('/scratch/dcr8536/Lifespan-Analysis/'))
subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};%, 'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
iterations = 1;
% surface
data_dir = '/scratch/dcr8536/TimeB/Nifti/postFCproc_CIFTI/FC_Parcels_333/';
atlas = 'Parcels333';
% surface nan
neg_corrs = 'nan';
output_str = 'neg_corrs_nan_Parcels333'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all
% surface zero
neg_corrs = 'zero';
output_str = 'neg_corrs_zero_Parcels333'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all
%surface as is
neg_corrs = 'asis';
output_str = 'neg_corrs_asis_Parcels333'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all
% volume
data_dir = '/scratch/dcr8536/TimeB/Nifti/preproc_FCProc/corrmats_Seitzman300/';
atlas = 'Seitzman300';
% volume nan
neg_corrs = 'nan';
output_str = 'neg_corrs_nan_Seitzman300'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all
% volume zero
neg_corrs = 'zero';
output_str = 'neg_corrs_zero_Seitzman300'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all
% volume as is
neg_corrs = 'asis';
output_str = 'neg_corrs_asis_Seitzman300'; %something to add to the filename for the output figures to differentiate it from others?
seg_index_reliability(subject, neg_corrs, data_dir, atlas, output_str, iterations);
close all





