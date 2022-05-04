load('/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/allsubs_seg_index_ses_longitudinal.mat')
longitudinal_SI = ses_SI;
longitudinal_SI(:,6:10) = [];
longitudinal_SI(1,4:5) = -1;
longitudinal_SI(:,1) = [];
clear rgb_for_plot
rgb_colors = [1 0 0;%LS02
            0, 1, 0;%LS03
            0, 0, 1];%LS05
for n = 1:3
    rgb_for_plot{n} = rgb_colors(n,:);
end
%% First plot longitudinal data -- unfilled
handles = plotSpread(longitudinal_SI', 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_plot', 'binWidth', .5)

axis([0 19 .25 .6])
 ax = gca;
 ax.FontSize = 20;
 
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
subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17','INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030'};

for n = 1:numel(subs)
rgb_for_plot{n} = rgb_colors(n,:);
end
load('/Users/dianaperez/Desktop/Segregation_Analyses/allsubs_seg_index_ses.mat');
lifespan_SI = ses_SI;
load('/Users/dianaperez/Desktop/Segregation_Analyses/INET_allsubs_seg_index_ses.mat');
iNet_SI =  ses_SI';
iNet_SI(5,:) = -1;

lifespan_SI(:,1) = [];
iNet_SI(1,:) = [];
SIs = [lifespan_SI' iNet_SI];

%% plot recent data -- add marker edge color: 'MarkerFaceColor', plotColors{iData,iCategory},...
handles = plotSpread([lifespan_SI' iNet_SI], 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_plot', 'xNames', subs, 'binWidth', 1)
hold on
%% plot session 1's
pre_first_ses = [0.40465483 0.47502443 0.42045888];
post_first_ses = [0.35866937 0.35903123 0.39794636 0.39622128 0.37635228 0.48386648 0.45387810 0.50429422 0.49092054 0.35866398 0.49722296 0.45307758 0.44976237 0.40819031 0.39566854];
handles = plotSpread([pre_first_ses post_first_ses], 'distributionMarkers', {'*'}, 'distributionColors', rgb_for_plot', 'xNames', subs)
pre_first_ses = [-1 -1 -1];

%% increase visibility -- delete marker edge color: 'MarkerFaceColor', plotColors{iData,iCategory},...
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
handles = plotSpread([lifespan_SI' iNet_SI], 'distributionMarkers', {'o'}, 'distributionColors', rgb_for_plot', 'xNames', subs, 'binWidth', 1)



                        
