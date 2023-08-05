%similarity comparisons
clear all
% all networks

LS_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/Lifespan_Parcels333_memory-default_similarityMat_matchedData.mat';
load(LS_path_to_data)
LS_subject = LS_subject; clear iNet_subject
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
LS_avg_between = mean(LS_between);
LS_avg_within = mean(LS_within);

iNet_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/iNetworks_Parcels333_memory-default_similarityMat_matchedData.mat';
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

iNet_avg_between = mean(iNet_between);
iNet_avg_within = mean(iNet_within);

[h_within,p_within] = ttest2(LS_within,iNet_within,'Vartype','unequal')
[h_between,p_between] = ttest2(LS_between,iNet_between,'Vartype','unequal')


%% COMPARE SYSTEM CLASSES

iNet_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/iNetworks_Parcels333_memory-default_similarityMat_matchedData.mat';
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

iNet_avg_between = mean(iNet_between);
iNet_avg_within = mean(iNet_within);
mem_wit = iNet_within;
mem_bet = iNet_between;


iNet_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/iNetworks_Parcels333_control-related_similarityMat_matchedData.mat';
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

iNet_avg_between = mean(iNet_between);
iNet_avg_within = mean(iNet_within);
con_wit = iNet_within;
con_bet = iNet_between;
iNet_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/iNetworks_Parcels333_sensorimotor_similarityMat_matchedData.mat';
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

iNet_avg_between = mean(iNet_between);
iNet_avg_within = mean(iNet_within);
SM_wit = iNet_within;
SM_bet = iNet_between;
[h_within,p_within] = ttest2(LS_within,iNet_within,'Vartype','unequal')

% control vs sensorimotor 
[hw1,pw1] = ttest(con_wit,SM_wit)
[hb1,pb1] = ttest(con_bet,SM_bet)
%control vs memory/default
[hw2,pw2] = ttest(con_wit,mem_wit)
[hb2,pb2] = ttest(con_bet,mem_bet)
%memory/default vs sensorimotor
[hw3,pw3] = ttest(mem_wit,SM_wit)
[hb3,pb3] = ttest(mem_bet,SM_bet)

LS_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/Lifespan_Parcels333_memory-default_similarityMat_matchedData.mat';
load(LS_path_to_data)
LS_subject = LS_subject; clear iNet_subject
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
LS_avg_between = mean(LS_between);
LS_avg_within = mean(LS_within);
mem_wit = LS_within;
mem_bet = LS_between;

LS_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/Lifespan_Parcels333_control-related_similarityMat_matchedData.mat';
load(LS_path_to_data)
LS_subject = LS_subject; clear iNet_subject
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
LS_avg_between = mean(LS_between);
LS_avg_within = mean(LS_within);
con_wit = LS_within;
con_bet = LS_between;

LS_path_to_data = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/similarity/Lifespan_Parcels333_sensorimotor_similarityMat_matchedData.mat';
load(LS_path_to_data)
LS_subject = LS_subject; clear iNet_subject
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
LS_avg_between = mean(LS_between);
LS_avg_within = mean(LS_within);
SM_wit = LS_within;
SM_bet = LS_between;

% control vs sensorimotor 
[hw1,pw1] = ttest(con_wit,SM_wit)
[hb1,pb1] = ttest(con_bet,SM_bet)
%control vs memory/default
[hw2,pw2] = ttest(con_wit,mem_wit)
[hb2,pb2] = ttest(con_bet,mem_bet)
%memory/default vs sensorimotor
[hw3,pw3] = ttest(mem_wit,SM_wit)
[hb3,pb3] = ttest(mem_bet,SM_bet)
