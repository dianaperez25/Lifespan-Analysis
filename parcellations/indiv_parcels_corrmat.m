% make corrmats for individual parcels

subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};

timecourse_dir = '/Volumes/fsmresfiles-1/PBS/Gratton_Lab/Lifespan/parcellations/avg_timecourses/';

%make corrmat

for sub = 1:numel(subs)
    load([timecourse_dir 'sub-' subs{sub} '_individual_parcels_average_timecourses.mat']);
    fout_str = sprintf('%s/corrmats/sub-%s_indiv_parcels_corrmats',timecourse_dir, subs{sub});
    corrmat = paircorr_mod(avg_parc_ts');
    save([fout_str '.mat'],'corrmat');
end