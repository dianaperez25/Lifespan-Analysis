%% compare methods for dealing with negative correlations
% I need to make a plot that will show the difference in segregation
% indices across three different methods of dealing with negative
% correlations and across volume and surface parcellations

clear all
%% some paths
data_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';
output_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';

%% some variables
% I don't think I'll need these, but just in case
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17','INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
rgb_colors = {[1 0 0], [0, 0, 1]};% red and blue
%% first load the first batch
neg_corrs = 'asis'; % choose: 'nan', 'zero', 'asis'
atlas = 'Seitzman300'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_vol_asis = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_vol_asis = OA_vol_asis.ses_SI; 
% these are the segregation indices for young adults
YA_vol_asis = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_vol_asis = YA_vol_asis.ses_SI;

%% load the next batch
neg_corrs = 'zero'; % choose: 'nan', 'zero', 'asis'
atlas = 'Seitzman300'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_vol_zero = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_vol_zero = OA_vol_zero.ses_SI; 
% these are the segregation indices for young adults
YA_vol_zero = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_vol_zero = YA_vol_zero.ses_SI;

%% load the next batch
neg_corrs = 'nan'; % choose: 'nan', 'zero', 'asis'
atlas = 'Seitzman300'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_vol_nan = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_vol_nan = OA_vol_nan.ses_SI; 
% these are the segregation indices for young adults
YA_vol_nan = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_vol_nan = YA_vol_nan.ses_SI;

%% load the first batch
neg_corrs = 'asis'; % choose: 'nan', 'zero', 'asis'
atlas = 'Parcels333'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_surf_asis = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_surf_asis = OA_surf_asis.ses_SI; 
% these are the segregation indices for young adults
YA_surf_asis = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_surf_asis = YA_surf_asis.ses_SI;

%% load the next batch
neg_corrs = 'zero'; % choose: 'nan', 'zero', 'asis'
atlas = 'Parcels333'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_surf_zero = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_surf_zero = OA_surf_zero.ses_SI; 
% these are the segregation indices for young adults
YA_surf_zero = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_surf_zero = YA_surf_zero.ses_SI;
YA_surf_zero = YA_surf_zero(:);
%% load the next batch
neg_corrs = 'nan'; % choose: 'nan', 'zero', 'asis'
atlas = 'Parcels333'; % 'Seitzman300' or 'Parcels333'
%% load all the files
% these are the TimeB segregation indices for older adults
OA_surf_nan = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
OA_surf_nan = OA_surf_nan.ses_SI; 
OA_surf_nan = OA_surf_nan(:);
% these are the segregation indices for young adults
YA_surf_nan = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
YA_surf_nan = YA_surf_nan.ses_SI(:);
YA_surf_nan = YA_surf_nan(:);
%% Now how do I plot these?






