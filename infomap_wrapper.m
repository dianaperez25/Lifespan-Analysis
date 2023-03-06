% infomap_wrapper
clear all
addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/Resources/'))
addpath(genpath('/scratch/dcr8536/Lifespan-Analysis/'))
subject = {'LS03', 'LS05', 'LS11'};
thresholdarray = [0.003 0.004 0.005:0.005:0.05];
for sub = 1:numel(subject)
outdir = ['/scratch/dcr8536/infomap/sub-' subject{sub} '/'];

% dmat
dmat_fname = ['/scratch/dcr8536/distances/' subject{sub} '_vertices_distmat.mat'];
dmat = smartload(dmat_fname);
% rmat
dconn_fname = ['/scratch/dcr8536/TimeB/Nifti/postFCproc_CIFTI/dconn_cifti_normalwall/Lifespan_Dconns/sub-' subject{sub} '_allsess_tmasked.dconn.nii'];
dconn = ft_read_cifti_mod(dconn_fname);
rmat = dconn.data(1:size(dmat,1), 1:size(dmat,1));
clear dconn



Run_Infomap_adaptive(rmat, dmat, 30, thresholdarray, 0, outdir, [],[],6,1);
clear dmat dconn
end