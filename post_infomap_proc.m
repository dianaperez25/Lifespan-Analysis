% post infomap proc
addpath(genpath('/projects/b1081/Scripts/graphtools/'))
minNetSize = 400;
sub = 'LS02';
str = ['sub-' sub '_infomap'];
cd(['/scratch/dcr8536/infomap/' sub '/'])
simple = modify_clrfile('simplify','rawassn.txt',minNetSize); %makes a file called rawassn_minsizeX.txt
regularized = rawoutput2clr(simple);
regularized(regularized < 2) = 0;
regularized = regularized - 1;
dlmwrite(['rawassn_minsize' num2str(minNetSize) '_regularized.txt'],regularized,'delimiter',' ') %makes a file called rawassn_minizeX_regularized.txt

thresholdarray = [0.003 0.004 0.005:0.005:0.05];
infoassn = dlmread(['rawassn_minsize' num2str(minNetSize) '_regularized.txt']);
conBensus(infoassn, str, [], thresholdarray*100, 'voxel', minNetSize)
load([str '_conBensus_weighted_minsize' num2str(minNetSize) '.mat'])
template_fname = '/scratch/dcr8536/template.dtseries.nii';
tmp = ft_read_cifti_mod(template_fname);
tmp.data = consen;%(tmp.brainstructure>0);
ft_write_cifti_mod([str '_conBensus_weighted_minsize' num2str(minNetSize) '.dtseries.nii'], tmp)

