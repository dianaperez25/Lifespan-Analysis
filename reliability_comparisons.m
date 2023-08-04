%reliability comparisons


LS_subject = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
iNet_subject = {'INET001', 'INET002', 'INET003', 'INET005', 'INET006','INET010',...
'INET018','INET019', 'INET026', 'INET030',  'INET032', 'INET033',...
'INET034', 'INET035', 'INET036', 'INET038', 'INET039', 'INET040', 'INET041',...
'INET042', 'INET043', 'INET044', 'INET045', 'INET046', 'INET047', 'INET048',...
'INET049', 'INET050', 'INET051', 'INET052', 'INET053', 'INET055', 'INET056',...
'INET057', 'INET058', 'INET059', 'INET060', 'INET061', 'INET062', 'INET063',...
'INET065', 'INET067', 'INET068', 'INET069', 'INET070', 'INET071', 'INET072', 'INET073'}; %

path_to_data = '/Users/diana/Desktop/Lifespan_reliability_v1.mat';
load(path_to_data);
LS_times = times_all;
LS_means = means;
LS_mean_of_means = mean_of_means;
LS_allsubs_corrs = allsubs_corrs;
path_to_data = '/Users/diana/Desktop/iNetworks_reliability_v1.mat';
load(path_to_data);
iNet_times = times_all;
iNet_means = means;
iNet_mean_of_means = mean_of_means;
iNet_allsubs_corrs = allsubs_corrs;
clear times_all means mean_of_means allsubs_corrs subject 

% compare peak reliability
LS_peak_reliability = [];
for sub = 1:numel(LS_subject);
    subs_reliability = LS_means{1,sub};
    LS_peak_reliability(sub,1) = max(subs_reliability);
end
    
iNet_peak_reliability = [];
for sub = 1:numel(iNet_subject);
    subs_reliability = iNet_means{1,sub};
    iNet_peak_reliability(sub,1) = max(subs_reliability);
end

% h = 0 does not reject the null hypothesis, 1 does
% p = p-value
[h1,p1] = ttest2(LS_peak_reliability,iNet_peak_reliability,'Vartype','unequal')

% comparison of the time it takes to reach .85
LS_time_to_reach = [];
for sub = 1:numel(LS_subject);
    subs_reliability = round(LS_means{1,sub},2);
    for t = 1:numel(subs_reliability)
        if subs_reliability(t) >= .85
            LS_time_to_reach(sub,1) = t*2.5;
            break;
        end
    end
end
    
iNet_time_to_reach = [];
for sub = 1:numel(iNet_subject);
    subs_reliability = round(iNet_means{1,sub},2);
    for t = 1:numel(subs_reliability)
        if subs_reliability(t) >= .85
            iNet_time_to_reach(sub,1) = t*2.5;
            break;
        end
    end
end

[h2,p2] = ttest2(LS_time_to_reach,iNet_time_to_reach,'Vartype','unequal')

