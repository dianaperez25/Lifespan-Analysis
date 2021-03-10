clear all
%272 = 4.99 min, 
dataDir = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_FCProc/corrmats_Seitzman300/';
subject = {'LS03', 'LS05'};
sessions = [5,5];
runs = [9,9,11,8,9;8,8,8,9,9];%LS03= 7,7,8,0,0;
%pts2sample = 8181; %8181 roughly equivalent to 150 minutes
pts2sample = 5440; %roughly equals 100 minutes
%pts2sample = 3780; %~70 minutes
sampstep=272; %5 minutes
atlas = 'Seitzman300';
atlas_dir = '/Users/dianaperez/Box/Quest_Backup/Atlases/';
output_dir = '/Users/dianaperez/Desktop/';


%% ROI info
atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);
roi_data = load_nii_wrapper(atlas_params.MNI_nii_file); %vox by 1

for s = 1:numel(subject)
    catData = [];
    catTmask = [];
    
    for i = 1:sessions
        
        %load mat file
        load([dataDir 'sub-' subject{s} '/sub-' subject{s} '_sess-' num2str(i) '_task-rest_corrmat_Seitzman300.mat'])

        masked_data = sess_roi_timeseries_concat(:,logical(tmask_concat'));

    %        catData = [catData masked_data];
        catData = [catData masked_data];
        catTmask = [catTmask tmask_concat'];
    end
    
    %I think this should be 10,816
    disp(sprintf('Total number of sample points for subject %s is %d by %d...', subject{s}, size(catData,1), size(catData,2)))


    %true-half is 150min=8181 samp points (TR = 1.1; (8181*1.1)/60=149.99)
    truehalf = catData(:,1:pts2sample);
    corrmat_truehalf = paircorr_mod(truehalf');
    
    fout_truehalf = sprintf('%s/sub-%s_truehalf_100min_corrmat_reliabilityanalysis',output_dir, subject{s});
    figure_corrmat_GrattonLab(corrmat_truehalf,atlas_params,-1,1);
    saveas(gcf,[fout_truehalf '.tiff'],'tiff');
    close(gcf);
    
    maskmat = ones(300);
    maskmat = logical(triu(maskmat, 1));
    truehalf_corrlin(1,:) = corrmat_truehalf(maskmat);

    times = floor((size(catData,2)-pts2sample)/sampstep);
    times = [5:5:(times*5)];
    %sample data
    for t = 1:numel(times)
        sampledDatas{t} = catData(:, (pts2sample+1):(pts2sample+(sampstep*t)));
        corrmat = paircorr_mod(sampledDatas{t}');
        fout = sprintf('%s/sub-%s_sampled_%dmin_corrmat_reliabilityanalysis', output_dir, subject{s}, times(t));
        figure_corrmat_GrattonLab(corrmat,atlas_params,-1,1);
        saveas(gcf,[fout '.tiff'],'tiff');
        close(gcf);
        corrs{t} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
        corrlins(t,:) = corrmat(maskmat);
    end

    sampledDatas{t+1} = catData(:, (pts2sample+1):end);
    corrmat = paircorr_mod(sampledDatas{t+1}');
    corrs{t+1} = paircorr_mod(triu(corrmat_truehalf), triu(corrmat));
    corrlins((t+1),:) = corrmat(maskmat);
    times = [5:5:((numel(times)+1)*5)];
    
    for j = 1:size(corrlins,1)
        tmpcorr = corrcoef(truehalf_corrlin', corrlins(j,:)', 'rows', 'complete');
        corr(s,j) = tmpcorr(2,1);
        clear tmpcorr
    end
    
    times_all(s,1:size(times,2)) = times;
    
    clear times
    clear catData
    clear catTmask
    clear corrlins
    clear corrmat
    clear corrmat_truehalf
    clear corrs
    clear masked_data
    clear sampledDatas
    clear tmask
    clear tmask_concat
    clear truehalf
    clear truehalf_corrlin
    clear sess_roi_timeseries
    clear sess_roi_timeseries_concat
    
end

times =[5:5:100];
figure;
% plot(times(1:5),corr(1,1:5),'Color',[1, 0.5, 0],'LineWidth', 3)
% ylim([0 1]);
% hold on
plot(times, corr(1,:), 'Color', [1,0.5,0], 'LineWidth', 3)
ylim([0 1]);
hold on
plot(times(1:18),corr(2,1:18),'Color',[0, 0, 1],'LineWidth', 3)
hold on

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');

m = findobj(gca,'Type','line');

hleg1 = legend(m(1:2), 'LS05', 'LS03', 'Location', 'SouthEast');
hleg1.FontSize = 14;
ax = gca;
ax.FontSize = 17;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);

print(gcf,[output_dir '/ReliabilityLifespanRestDatatruhalf' num2str(pts2sample) '.jpg'],'-dpng','-r300');





