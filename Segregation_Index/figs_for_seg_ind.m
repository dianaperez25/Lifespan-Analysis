%% By subject
% NaN
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/INET_allsubs_seg_index_sub_withNaNneg.mat')
iNet_SI =  sub_SI';
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/LS_allsubs_seg_index_sub_withNaNneg.mat')
lifespan_SI = sub_SI';
ind=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];
handles = plotSpread([iNet_SI; lifespan_SI], 'distributionIdx', ind,'distributionMarkers', {'x'},'xNames', {'YA', 'OA'}, 'distributionColors', {'b', 'r'},'binWidth', 1)
hold on
% zeros
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/INET_allsubs_seg_index_sub_with0neg.mat')
iNet_SI =  sub_SI';
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/LS_allsubs_seg_index_sub_with0neg.mat')
lifespan_SI = sub_SI';
handles = plotSpread([iNet_SI; lifespan_SI], 'distributionIdx', ind,'distributionMarkers', {'o'},'xNames', {'YA', 'OA'}, 'distributionColors', {'b', 'r'},'binWidth', 1)
ax = gca;
ax.FontSize = 24;
axis([0, 3, 0.3, .8])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.7]); %first and second control position on screen, third controls width, and fourth controls height
print(gcf, '/Users/dianaperez/Desktop/SI_by_sub_NanVs0.jpg', '-dpng', '-r300')

%% By session
% Nan
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/allsubs_seg_index_ses.mat')
lifespan_SI = ses_SI;
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/INET_allsubs_seg_index_ses.mat')
iNet_SI = ses_SI;
iNet_SI(:,5) = -1;
subs = {'INET001','INET002', 'INET003','INET005','INET006','INET010','INET016','INET018','INET019','INET030', 'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16','LS17'};
handles = plotSpread([iNet_SI; lifespan_SI]', 'categoryIdx', [ind;ind;ind;ind;ind],'categoryMarkers', {'x','x'},'xNames', subs, 'categoryColors', {'b', 'r'},'binWidth', 1)
hold on
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/LS_allsubs_seg_index_ses_with0negs.mat')
lifespan_SI = ses_SI;
load('/Users/dianaperez/Desktop/Research/Segregation_Analyses/INET_allsubs_seg_index_ses_with0negs.mat')
iNet_SI = ses_SI;
iNet_SI(:,5) = -1;
handles = plotSpread([iNet_SI; lifespan_SI]', 'categoryIdx', [ind;ind;ind;ind;ind],'categoryMarkers', {'o','o'},'xNames', subs, 'categoryColors', {'b', 'r'},'binWidth', 1)
ax = gca;
ax.FontSize = 24;
axis([0.5, 18.5, .3, .8])
print(gcf, '/Users/dianaperez/Desktop/SI_by_ses_NanVs0.jpg', '-dpng', '-r300')


%% by network
% Nan
data = [];
avg = [];
nets = {'unassign' 'SMd'	'SMl' 'CO' 'Aud' 'DMN'	'PMN' 'Vis'	'FP' 'Sal' 'Lang' 'DAN' 'MTL' 'PON'};

for sub = 1:numel(subs)
    load(['/Users/dianaperez/Desktop/Research/Segregation_Analyses/' subs{sub} '_seg_index_net_withNaNnegs.mat'])
    sub_struct.seg_ind(15,:) = [];
    avg = [avg mean(sub_struct.seg_ind, 2)];
    if contains(subs{sub}, 'INET')
        sub_struct.seg_ind(:,5) = -1
    end
    data = [data sub_struct.seg_ind];
end
data = data';
avg_nan = avg';
cx_ind = [];
cx_ind(1:50) = 1; 
cx_ind(51:90) = 2;
cx_ind = repmat(cx_ind, 14,1);
%ind = repmat([1:14]',18,1);

handles = plotSpread(data, 'categoryIdx', cx_ind','categoryMarkers', {'x','x'},'xNames', nets, 'categoryColors', {'b', 'r'},'binWidth', 1)
hold on

% 0
data = [];
avg = [];
for sub = 1:numel(subs)
    load(['/Users/dianaperez/Desktop/Research/Segregation_Analyses/' subs{sub} '_seg_index_net_with0negs.mat'])
    sub_struct.seg_ind(15,:) = [];
    avg = [avg mean(sub_struct.seg_ind, 2)];
    if contains(subs{sub}, 'INET')
        sub_struct.seg_ind(:,5) = -1
    end
    data = [data sub_struct.seg_ind];
end
avg_zero = avg';
data = data';
handles = plotSpread(data, 'categoryIdx', cx_ind','categoryMarkers', {'.','.'},'xNames', nets, 'categoryColors', {'b', 'r'},'binWidth', 1)
ax = gca;
ax.FontSize = 24;
axis([0, 15, -0.3, .9])
print(gcf, '/Users/dianaperez/Desktop/SI_by_net_allsess_NanVs0.jpg', '-dpng', '-r300')
close gcf
% only averages across sessions
cx_ind = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]';
cx_ind = repmat(cx_ind, 1,14);
handles = plotSpread(avg_nan, 'categoryIdx', cx_ind,'categoryMarkers', {'x','x'},'xNames', nets, 'categoryColors', {'b', 'r'},'binWidth', 1)
hold on
handles = plotSpread(avg_zero, 'categoryIdx', cx_ind,'categoryMarkers', {'o','o'},'xNames', nets, 'categoryColors', {'b', 'r'},'binWidth', 1)
print(gcf, '/Users/dianaperez/Desktop/SI_by_net_NanVs0.jpg', '-dpng', '-r300')
