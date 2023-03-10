%% Stability of Segregation Index

output_dir = '/Users/dianaperez/Desktop/Segregation_Analyses/';
if ~exist(output_dir)
    mkdir(output_dir)
end
load('/Users/dianaperez/Desktop/Segregation_Analyses/INET_allsubs_seg_index_ses.mat')
diff_score = [];
count = 1;
for sub1 = 1:10
    for ses1 = 1:4
        for sub2 = 1:10
            for ses2 = 1:4
                diff_score(count) = ses_SI(sub1,ses1) - ses_SI(sub2,ses2);
                count = count + 1;
            end
        end
    end
end
matrix = abs(diff_score);
matrix = reshape(matrix, [40,40]);
figure('Position',[1 1 1000 800]);
imagesc(matrix,[0 .14]); colormap('jet');
hline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
vline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
set(gca,'XTick',[3,8,13,18,23,28,33,38], 'YTick', [3,8,13,18,23,28,33,38], 'XTickLabel',...
     {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'}, 'YTickLabel', {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'});
axis square;
colorbar;
%title('Correlation Matrix Similarity');
saveas(gcf,[output_dir 'SimilarityMat_SegregationIndex.tiff'],'tiff');

% FT_matrix = FisherTransform(matrix);
% figure('Position',[1 1 1000 800]);
% imagesc(matrix,[0 .14]); colormap('jet');
% hline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
% vline_new([0,5,10,15,20,25,30,35,40]+0.5,'k',2);
% set(gca,'XTick',[3,8,13,18,23,28,33,38], 'YTick', [3,8,13,18,23,28,33,38], 'XTickLabel',...
%      {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'}, 'YTickLabel', {'LS02', 'LS03', 'LS05', 'LS08', 'LS11', 'LS14','LS16','LS17'});
% axis square;
% colorbar;
