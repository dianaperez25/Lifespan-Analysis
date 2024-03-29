%% This script plots segregation indices by session for Lifespan and iNet 

clear all
%% some paths
data_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';
output_dir = '/Users/dianaperez/Desktop/Research/Segregation_Analyses/';

%% some options
neg_corrs = 'asis'; % choose: 'nan', 'zero', 'asis'
atlas = 'Seitzman300'; % 'Seitzman300' or 'Parcels333'

%% some variables
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17','INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};
rgb_colors = [1 0 0;%LS02
            0, 1, 0;%LS03
            0, 0, 1;%LS05
            0, 1, 1;%LS08
            1, 0, 1;%LS11
            0.4660 0.6740 0.188;%LS14
            0.9290 0.6940 0.1250;%LS16
            0.4940 0.1840 0.5560;%LS17
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0];
        
        
%% load all the files
% these are the TimeB segregation indices for older adults
all_OA_subs = load([data_dir 'Lifespan_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
all_OA_subs = all_OA_subs.ses_SI; 
% these are the TimeA segregation indices for older adult longitudinal subjects
long_OA_subs = load([data_dir 'Lifespan_longitudinalsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
long_OA_subs = long_OA_subs.ses_SI;
% these are the segregation indices for young adults
all_YA_subs = load([data_dir 'iNetworks_allsubs_seg_index_ses_negcorrs' neg_corrs '_' atlas '.mat']);
all_YA_subs = all_YA_subs.ses_SI;
 % iNet subs do not have 5th sessions so setting to -1 so that they don't show up in the plot

%% some data curation
% setting some sessions to negative 1 so that they don't show up in the plot (subjects did not complete those sessions)
long_OA_subs(1,4:5) = -1;
all_YA_subs(:,5) = -1;
% delete first session because we will plot those separately; but FIRST
% save to new variable
first_ses_OA = all_OA_subs(:,1);
all_OA_subs(:,1) = [];
first_ses_YA = all_YA_subs(:,1);
all_YA_subs(:,1) = [];
first_ses_long_OA = long_OA_subs(:,1);
long_OA_subs(:,1) = [];
% concat OA and YA data to plot
all_subs = [all_OA_subs' all_YA_subs']; % all segregation indices together

for n = 1:numel(subs)
rgb_for_plot{n} = rgb_colors(n,:);
end

rgb_for_long = rgb_for_plot(:,1:3);
%% First, plot Time A data
% session 2+ data with unfilled circles
handles = plotSpread(long_OA_subs', 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_long', 'binWidth', .5);
hold on
% session 1 data with a cross
handles = plotSpread(first_ses_long_OA', 'distributionMarkers', {'+'}, 'distributionColors', rgb_for_long', 'xNames', subs);
clear rgb_for_long
hold on
%% Then, plot Time B data with filled circles 
handles = plotSpread_v2(all_subs, 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_plot', 'xNames', subs, 'binWidth', 1);
hold on
%% Then, plot session 1 data
handles = plotSpread([first_ses_OA' first_ses_YA'], 'distributionMarkers', {'*'}, 'distributionColors', rgb_for_plot', 'xNames', subs);
hold on
%% Lastly, we are going to increase visibility of the circles
% setting rgb colors to blue for light colors and white for dark colors
rgb_colors = [0 0 0;
            0 0 0;
            1 1 1;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1;
            1 1 1];
for n = 1:numel(subs)
    rgb_for_plot{n} = rgb_colors(n,:);
end
% then drawing unfilled circles around existing data points
handles = plotSpread(all_subs, 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_plot', 'xNames', subs, 'binWidth', 1);

% some parameters for the plot
axis([0 19 .25 .6]) % for the nan neg corrs seitzman300
axis([0 19 .8 1.05]) % for the as is neg corrs parcels333
axis([0 19 0.6 0.8]) % for the zero neg corrs parcels333
axis([0 19 0.4 0.6]) % for the nan neg corrs parcels333
axis([0 19 .5 .8]) % for the zero neg corrs seitzman300
axis([0 19 .85 1.15]) % for the as is neg corrs seitzman300
ax = gca;
ax.FontSize = 20;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 1, 0.7]);
print(gcf, [output_dir 'SegInd_stability_negcorrs' neg_corrs '_' atlas '.jpg'], '-dpng', '-r300')
close all                        
