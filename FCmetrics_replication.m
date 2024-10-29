
%ohbm analyses
datasets = {'Lifespan-NU', 'iNet-NU', 'Lifespan-FSU', 'iNet-FSU'};
exclude_subs = {'INET001', 'INET061', 'LS46', 'INET108'};

parcel_data_dir = '/Users/dianaperez/Desktop/FC_metrics_replication/new/';

for d = 1:numel(datasets)
    group_fname = sprintf('%s/group_%s_segregation_index_group_parcels.mat', parcel_data_dir, datasets{d});
    indiv_fname = sprintf('%s/group_%s_segregation_index_indiv_parcels.mat', parcel_data_dir, datasets{d});

    load(group_fname); 
    if strcmpi(datasets{d}, 'Lifespan-NU')
        OA_NU.group = group_struct;
    elseif strcmpi(datasets{d}, 'iNet-NU')
        YA_NU.group = group_struct;
    elseif strcmpi(datasets{d}, 'iNet-FSU')
        YA_FSU.group = group_struct;
    elseif strcmpi(datasets{d}, 'Lifespan-FSU')
        OA_FSU.group = group_struct;
    end

    load(indiv_fname);
    if strcmpi(datasets{d}, 'Lifespan-NU')
        OA_NU.indiv = group_struct;
    elseif strcmpi(datasets{d}, 'iNet-NU')
        YA_NU.indiv = group_struct;
    elseif strcmpi(datasets{d}, 'iNet-FSU')
        YA_FSU.indiv = group_struct;
    elseif strcmpi(datasets{d}, 'Lifespan-FSU')
        OA_FSU.indiv = group_struct;
    end
    clear group_struct
end

% set up data for anova

% values for individualized parcellations
seg_indiv = [OA_NU.indiv.averages(:,1); YA_NU.indiv.averages(:,1); OA_FSU.indiv.averages(:,1); YA_FSU.indiv.averages(:,1)];
wFC_indiv = [OA_NU.indiv.averages(:,2); YA_NU.indiv.averages(:,2); OA_FSU.indiv.averages(:,2); YA_FSU.indiv.averages(:,2)];
bFC_indiv = [OA_NU.indiv.averages(:,3); YA_NU.indiv.averages(:,3); OA_FSU.indiv.averages(:,3); YA_FSU.indiv.averages(:,3)];

% values for group parcellations
seg_group = [OA_NU.group.averages(:,1); YA_NU.group.averages(:,1); OA_FSU.group.averages(:,1); YA_FSU.group.averages(:,1)];
wFC_group = [OA_NU.group.averages(:,2); YA_NU.group.averages(:,2); OA_FSU.group.averages(:,2); YA_FSU.group.averages(:,2)];
bFC_group = [OA_NU.group.averages(:,3); YA_NU.group.averages(:,3); OA_FSU.group.averages(:,3); YA_FSU.group.averages(:,3)];

% put them together
seg = [seg_indiv; seg_group];
wFC = [wFC_indiv; wFC_group];
bFC = [bFC_indiv; bFC_group];

data = [seg wFC bFC];

% set up factor labels
age_group = [];
age_group(1:8) = 1; % older adults NU
age_group(9:54) = 2; % younger adults NU
age_group(55:77) = 1; % older adults FSU
age_group(78:120) = 2; % younger adults FSU
age_group = [age_group'; age_group'];

networks = [];
networks(1:120) = 3; % individualized parcels
networks(121:240) = 4; % group parcels

site = [];
site(1:54) = 5; % NU
site(55:120) = 6; % group parcels
site = [site'; site'];

factors = {age_group, networks', site};
p_vals = [];
% anova for segregation index
aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation" "Site"])
groupmeans(aov_seg,["Age Group" "Parcellation" "Site"])
stats_aov = stats(aov_seg)
p_vals = [p_vals; stats_aov.pValue];

aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation" "Site"])
groupmeans(aov_wFC,["Age Group" "Parcellation" ])
stats_aov = stats(aov_wFC)
p_vals = [p_vals; stats_aov.pValue];

aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation" "Site"])
groupmeans(aov_bFC,["Age Group" "Parcellation"])
stats_aov = stats(aov_bFC)
p_vals = [p_vals; stats_aov.pValue];



%% VERTEX LEVEL DATA

load('/Users/dianaperez/Desktop/ohbm_poster/vertex_level/all_subs_info_vertex.mat')

% set up data for anova

% values for individualized parcellations
seg_indiv = indiv_networks.FC_metrics(:,1);
wFC_indiv = indiv_networks.FC_metrics(:,2);
bFC_indiv = indiv_networks.FC_metrics(:,3);

% values for group parcellations
segind_group = group_networks.FC_metrics(:,1);
wFC_group = group_networks.FC_metrics(:,2);
bFC_group = group_networks.FC_metrics(:,3);

% put them together
seg = [seg_indiv; segind_group];
wFC = [wFC_indiv; wFC_group];
bFC = [bFC_indiv; bFC_group];

data = [seg wFC bFC];

% set up factor labels
age_group = [];
age_group(1:8) = 1; % older adults
age_group(9:26) = 2; % younger adults
age_group = [age_group'; age_group'];

networks = [];
networks(1:26) = 3; % individualized networks
networks(27:52) = 4; % group networks

factors = {age_group, networks'};

p_vals = []
% anova for segregation index
aov_seg = anova(factors, seg, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation"])
groupmeans(aov_seg,["Age Group" "Parcellation"])
stats_aov = stats(aov_seg)
p_vals = [p_vals; stats_aov.pValue];

aov_wFC = anova(factors, wFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation"])
groupmeans(aov_wFC,["Age Group" "Parcellation"])
stats_aov = stats(aov_wFC)
p_vals = [p_vals; stats_aov.pValue];

aov_bFC = anova(factors, bFC, 'ModelSpecification', 'interactions', FactorNames=["Age Group" "Parcellation"])
groupmeans(aov_bFC,["Age Group" "Parcellation"])
stats_aov = stats(aov_bFC)
p_vals = [p_vals; stats_aov.pValue];


[p_fdr, p_masked] = FDR(p_vals, 0.05);

% t-tests

%segregation
oa_seg_indiv = seg_indiv(1:8);
oa_seg_group = segind_group(1:8);
ya_seg_indiv = seg_indiv(9:end);
ya_seg_group = segind_group(9:end);

% comparison of young vs older adults using individual parcels: significant
[h,p,ci,stats] = ttest2(ya_seg_indiv,oa_seg_indiv, 'VarType', 'unequal')
p_vals = [p_vals; p];

% comparison of young vs older adults using group parcels: significant
[h,p,ci,stats] = ttest2(ya_seg_group,oa_seg_group, 'VarType', 'unequal')
p_vals = [p_vals; p];

% comparison of individualized and group parcels in young adults:
% significant
[h,p,ci,stats] = ttest2(ya_seg_indiv,ya_seg_group)
p_vals = [p_vals; p];

% comparison of individualized and group parcels in older adults:
% significant
[h,p,ci,stats] = ttest2(oa_seg_indiv,oa_seg_group)
p_vals = [p_vals; p];

% comparison of individualized and group parcels in all:
% significant
[h,p,ci,stats] = ttest2(seg_indiv,segind_group)
p_vals = [p_vals; p];

% wFC
oa_wFC_indiv = wFC_indiv(1:8);
oa_wFC_group = wFC_group(1:8);
ya_wFC_indiv = wFC_indiv(9:end);
ya_wFC_group = wFC_group(9:end);

% comparison of young vs older adults using individual parcels: significant
[h,p,ci,stats] = ttest2(ya_wFC_indiv,oa_wFC_indiv, 'VarType', 'unequal')
p_vals = [p_vals; p];

% comparison of young vs older adults using group parcels: significant
[h,p,ci,stats] = ttest2(ya_wFC_group,oa_wFC_group, 'VarType', 'unequal')
p_vals = [p_vals; p];

%significant
[h,p,ci,stats] = ttest2(ya_wFC_indiv,ya_wFC_group)
p_vals = [p_vals; p];

% non-significant
[h,p,ci,stats] = ttest2(oa_wFC_indiv,oa_wFC_group)
p_vals = [p_vals; p];

% comparison of individualized and group parcels in all:
% significant
[h,p,ci,stats] = ttest2(wFC_indiv,wFC_group)
p_vals = [p_vals; p];

% bFC
oa_bFC_indiv = bFC_indiv(1:8);
oa_bFC_group = bFC_group(1:8);
ya_bFC_indiv = bFC_indiv(9:end);
ya_bFC_group = bFC_group(9:end);

%non-significant
[h,p,ci,stats] = ttest2(ya_bFC_indiv,oa_bFC_indiv, 'VarType', 'unequal')
p_vals = [p_vals; p];

%non-significant
[h,p,ci,stats] = ttest2(ya_bFC_group,oa_bFC_group, 'VarType', 'unequal')
p_vals = [p_vals; p];

%significant
[h,p,ci,stats] = ttest2(ya_bFC_indiv,ya_bFC_group)
p_vals = [p_vals; p];

%not significant
[h,p,ci,stats] = ttest2(oa_bFC_indiv,oa_bFC_group);
p_vals = [p_vals; p];

[p_fdr, p_masked] = FDR(p_vals, 0.05);


%permutations
%let's first compare older adult indiv vs group
p_vals_perms = [];
%segregation index
oa_seg_indiv = seg_indiv(1:8);
oa_seg_group = segind_group(1:8);
oa_seg = [oa_seg_indiv; oa_seg_group];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = oa_seg(rand_labels==1);
    pseudo_group = oa_seg(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(oa_seg_indiv - oa_seg_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);
% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
%young adults
ya_seg_indiv = seg_indiv(9:end);
ya_seg_group = segind_group(9:end);
ya_seg = [ya_seg_indiv; ya_seg_group];
labels = ones(1,36);
labels(19:end) = 2;
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = ya_seg(rand_labels==1);
    pseudo_group = ya_seg(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(ya_seg_indiv - ya_seg_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);

% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
% diff between indiv - group compared across YA and OA
diff = seg_indiv- segind_group;
labels = ones(1,26);
labels(9:end) = 2;

num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_oa = diff(rand_labels==1);
    pseudo_ya = diff(rand_labels==2);
    perm_diffs(p) = mean(pseudo_oa) - mean(pseudo_ya);
end
true_diff = mean(oa_seg_indiv - oa_seg_group) - mean(ya_seg_indiv - ya_seg_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);

% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
%within net FC
oa_wFC_indiv = wFC_indiv(1:8);
oa_wFC_group = wFC_group(1:8);

oa_wFC = [oa_wFC_indiv; oa_wFC_group];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = oa_wFC(rand_labels==1);
    pseudo_group = oa_wFC(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(oa_wFC_indiv - oa_wFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);
% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
%young adults
ya_wFC_indiv = wFC_indiv(9:end);
ya_wFC_group = wFC_group(9:end);
ya_wFC = [ya_wFC_indiv; ya_wFC_group];
labels = ones(1,36);
labels(19:end) = 2;
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = ya_wFC(rand_labels==1);
    pseudo_group = ya_wFC(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(ya_wFC_indiv - ya_wFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);

% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
% diff between indiv - group compared across YA and OA
diff = wFC_indiv- wFC_group;
labels = ones(1,26);
labels(9:end) = 2;

num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_oa = diff(rand_labels==1);
    pseudo_ya = diff(rand_labels==2);
    perm_diffs(p) = mean(pseudo_oa) - mean(pseudo_ya);
end
true_diff = mean(oa_wFC_indiv - oa_wFC_group) - mean(ya_wFC_indiv - ya_wFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);

% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
%between net FC
oa_bFC_indiv = bFC_indiv(1:8);
oa_bFC_group = bFC_group(1:8);

oa_bFC = [oa_bFC_indiv; oa_bFC_group];

labels = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2];
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = oa_bFC(rand_labels==1);
    pseudo_group = oa_bFC(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(oa_bFC_indiv - oa_bFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);
% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
%young adults
ya_bFC_indiv = bFC_indiv(9:end);
ya_bFC_group = bFC_group(9:end);
ya_bFC = [ya_bFC_indiv; ya_bFC_group];
labels = ones(1,36);
labels(19:end) = 2;
num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_indiv = ya_bFC(rand_labels==1);
    pseudo_group = ya_bFC(rand_labels==2);
    perm_diffs(p) = mean(pseudo_indiv - pseudo_group);
end
true_diff = mean(ya_bFC_indiv - ya_bFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);

% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
p_vals_perms = [p_vals_perms; p];
% diff between indiv - group compared across YA and OA
diff = bFC_indiv- bFC_group;
labels = ones(1,26);
labels(9:end) = 2;

num_perms = 1000;
for p = 1:num_perms
    rng('shuffle');
    ind = randperm(length(labels))';
    rand_labels = labels(ind);
    pseudo_oa = diff(rand_labels==1);
    pseudo_ya = diff(rand_labels==2);
    perm_diffs(p) = mean(pseudo_oa) - mean(pseudo_ya);
end
true_diff = mean(oa_bFC_indiv - oa_bFC_group) - mean(ya_bFC_indiv - ya_bFC_group);
p = min([length(find(perm_diffs<true_diff))/num_perms, length(find(perm_diffs>true_diff))/num_perms]);
p_vals_perms = [p_vals_perms; p];
% ind = zeros(num_perms+1,1);
% ind(1,1) = 1;
% handles = plotSpread([true_diff; perm_diffs'], 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'})
% scatter(1, true_diff, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
% 

[p_fdr_perms, p_masked] = FDR(p_vals_perms, 0.05)



% plotting with individual lines
oa_pairs_wFC = [oa_wFC_group oa_wFC_indiv];
ya_pairs_wFC = [ya_wFC_group ya_wFC_indiv];
means_wFC = [mean(oa_wFC_group) mean(oa_wFC_indiv) mean(ya_wFC_group) mean(ya_wFC_indiv)] ;
% ya_means_wFC = [mean(ya_wFC_group) mean(ya_wFC_indiv)] ;
%standard error
err = [(std(oa_wFC_group)/sqrt(8)) (std(oa_wFC_indiv)/sqrt(8)) (std(ya_wFC_group)/sqrt(18)) (std(ya_wFC_indiv)/sqrt(18))];
figure;
colors = [1.00000, 0.82353, 0.30196; 1.00000, 0.90196, 0.60000; 1.00000, 0.82353, 0.30196; 1.00000, 0.90196, 0.60000];
for m = 1:length(means_wFC)
    if m < 3
        bg(m) = bar(m, means_wFC(m), 'FaceColor', colors(m,:))
        hold on
        errorbar(m, means_wFC(m), err(m), 'Color', 'k')
    else
        bg(m) = bar(m+1, means_wFC(m), 'FaceColor', colors(m,:))
        hold on
        errorbar(m+1, means_wFC(m), err(m), 'Color', 'k')
    end
end

hold on

for sub = 1:length(oa_pairs_wFC)
    plot(1:2,[oa_pairs_wFC])
    hold on
end

hold on

for sub = 1:length(ya_pairs_wFC)
    plot(4:5,[ya_pairs_wFC])
    hold on
end

xticks([1.5 4.5])
xticklabels({'Older Adults', 'Younger Adults'})
ax = gca;
ax.FontSize = 18;



oa_pairs_bFC = [oa_bFC_group oa_bFC_indiv];
ya_pairs_bFC = [ya_bFC_group ya_bFC_indiv];
means_bFC = [mean(oa_bFC_group) mean(oa_bFC_indiv) mean(ya_bFC_group) mean(ya_bFC_indiv)] ;

%standard error
err = [(std(oa_bFC_group)/sqrt(8)) (std(oa_bFC_indiv)/sqrt(8)) (std(ya_bFC_group)/sqrt(18)) (std(ya_bFC_indiv)/sqrt(18))];
figure;
colors = [0.20000, 0.60000, 1.00000; 0.50196, 0.74902, 1.00000; 0.20000, 0.60000, 1.00000; 0.50196, 0.74902, 1.00000];
for m = 1:length(means_bFC)
    if m < 3
        bg(m) = bar(m, means_bFC(m), 'FaceColor', colors(m,:))
        hold on
        errorbar(m, means_bFC(m), err(m), 'Color', 'k')
    else
        bg(m) = bar(m+1, means_bFC(m), 'FaceColor', colors(m,:))
        hold on
        errorbar(m+1, means_bFC(m), err(m), 'Color', 'k')
    end
end

hold on

for sub = 1:length(oa_pairs_bFC)
    plot(1:2,[oa_pairs_bFC])
    hold on
end

hold on

for sub = 1:length(ya_pairs_bFC)
    plot(4:5,[ya_pairs_bFC])
    hold on
end

xticks([1.5 4.5])
xticklabels({'Older Adults', 'Younger Adults'})
ax = gca;
ax.FontSize = 18;
xlabel('Age Group')
ylabel('Average Functional Connectivity')






oa_pairs_seg = [oa_seg_group oa_seg_indiv];
ya_pairs_seg = [ya_seg_group ya_seg_indiv];
means_seg = [mean(oa_seg_group) mean(oa_seg_indiv) mean(ya_seg_group) mean(ya_seg_indiv)] ;

%standard error
err = [(std(oa_seg_group)/sqrt(8)) (std(oa_seg_indiv)/sqrt(8)) (std(ya_seg_group)/sqrt(18)) (std(ya_seg_indiv)/sqrt(18))];
figure;
colors = [0.32549, 0.77647, 0.32549; 0.54902, 0.85098, 0.54902; 0.32549, 0.77647, 0.32549; 0.54902, 0.85098, 0.54902];
for m = 1:length(means_seg)
    if m < 3
        bg(m) = bar(m, means_seg(m), 'FaceColor', colors(m,:), 'FaceAlpha', .5)
        hold on
        errorbar(m, means_seg(m), err(m), 'Color', 'k')
    else
        bg(m) = bar(m+1, means_seg(m), 'FaceColor', colors(m,:), 'FaceAlpha', .5)
        hold on
        errorbar(m+1, means_seg(m), err(m), 'Color', 'k')
    end
end

hold on

for sub = 1:length(oa_pairs_seg)
    plot(1:2,[oa_pairs_seg])
    hold on
end

hold on

for sub = 1:length(ya_pairs_seg)
    plot(4:5,[ya_pairs_seg])
    hold on
end

xticks([1.5 4.5])
xticklabels({'Older Adults', 'Younger Adults'})
ax = gca;
ax.FontSize = 18;
xlabel('Age Group')
ylabel('Segregation Index')