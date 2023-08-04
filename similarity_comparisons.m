%similarity comparisons

% all networks

LS_path_to_data = '/Users/diana/Desktop/Lifespan_surf_similarityMat_matchedData.mat';
load(LS_path_to_data)
LS_subject = iNet_subject; clear iNet_subject
sessions = 5;
count = 1;
LS_within = [];
LS_between = [];
for s = 1:numel(LS_subject)
    lines = [count:(count+sessions-1)];
    sub_vals = simmat(lines,:);
    maskmat = ones(sessions,sessions);
    maskmat = logical(triu(maskmat, 1));
    within_sub = sub_vals(:,lines);
    LS_within = [LS_within; mean(within_sub(maskmat))];
    maskmat = ones(size(sub_vals));
    maskmat(:,lines) = 0;
    LS_between = [LS_between; mean(sub_vals(maskmat==1))];
    count = count+sessions;
end
%mean_between = mean(between);
%mean_within = mean(within);

iNet_path_to_data = '/Users/diana/Desktop/iNet_simmat_info.mat';
load(iNet_path_to_data)
sessions = 4;
count = 1;
iNet_within = [];
iNet_between = [];
for s = 1:numel(iNet_subject)
    lines = [count:(count+sessions-1)];
    sub_vals = simmat(lines,:);
    maskmat = ones(sessions,sessions);
    maskmat = logical(triu(maskmat, 1));
    within_sub = sub_vals(:,lines);
    iNet_within = [iNet_within; mean(within_sub(maskmat))];
    maskmat = ones(size(sub_vals));
    maskmat(:,lines) = 0;
    iNet_between = [iNet_between; mean(sub_vals(maskmat==1))];
    count = count+sessions;
end

[h_within,p_within] = ttest2(LS_within,iNet_within,'Vartype','unequal')
[h_between,p_between] = ttest2(LS_between,iNet_between,'Vartype','unequal')

