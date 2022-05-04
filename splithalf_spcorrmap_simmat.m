%% spatial corr map similarity split-half

subs = {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14', 'LS16', 'LS17'};
data_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/spCor_Maps/split-half/';
out_dir = '/Volumes/RESEARCH_HD/Lifespan/CNS_analyses/';
corrs = [];
for sub1 = 1:numel(subs)
    for sub2 = 1:numel(subs)
    map1 = ft_read_cifti_mod([data_dir 'sub-' subs{sub1} '_first-half_vs_120_avg_corr_LR_cortex_corr.dtseries.nii']);
    map2 = ft_read_cifti_mod([data_dir 'sub-' subs{sub2} '_second-half_vs_120_avg_corr_LR_cortex_corr.dtseries.nii']);
    corrs(sub1,sub2) = corr(map1.data, map2.data);
    end
end


figure('Position',[1 1 1000 800]);
imagesc(corrs,[0 1]); colormap('jet');
hline_new([0:8]+0.5,'k',2);
vline_new([0:8]+0.5,'k',2);
set(gca,'XTick',[1:8], 'YTick', [1:8], 'XTickLabel',...
     {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'}, 'YTickLabel', {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'});
axis square;
colorbar;
%title('Correlation Matrix Similarity');
saveas(gcf,[out_dir 'SimilarityMat_SpCorrMap.tiff'],'tiff');
