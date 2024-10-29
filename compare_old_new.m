% Ok, I need to check which subjects have changed so much after the post freesurfer pipeline was fixed...
%     I have the old FCParcels333 data, and I'm going to measure the correlation between the old and new matrices

datasets = {'Lifespan-NU', 'iNet-NU'};
exclude_subs = {};
old_data_dir = '/Users/dianaperez/Downloads/old_FC_Parcels_333/';
new_data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Diana/Diss/Nifti/FC_Parcels_333/';

maskmat = ones(333);
maskmat = logical(triu(maskmat, 1));

for d = 1:numel(datasets)
    [subjects, sessions, N] = get_subjects(datasets{d}, exclude_subs);
    count = 1; ses_lims = []; old_corrlin = []; new_corrlin = [];
    for s = 1:numel(subjects)
        ses_count = 0;
        for t = 1:sessions            
            old_data_fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', old_data_dir, subjects{s}, t);
            new_data_fname = sprintf('%s/sub-%s_rest_ses-%d_parcel_corrmat.mat', new_data_dir, subjects{s}, t);
            old_data = load(old_data_fname); old_data = old_data.parcel_corrmat;
            new_data = load(new_data_fname); new_data = new_data.parcel_corrmat;
            old_corrlin(count, :) = old_data(maskmat);
            new_corrlin(count, :) = new_data(maskmat);
            count = count + 1;
        end
        ses_lims = [ses_lims; ses_count];
    end
    simmat = corr(old_corrlin', new_corrlin');

    figure('Position',[1 1 1000 800]);
load better_jet_colormap.mat
imagesc(simmat,[0 1]); colormap(better_jet_colormap_diff);
tick_marks = [0:sessions:(5*numel(subjects))]+0.5;
hline_new(tick_marks,'k',2);
vline_new(tick_marks,'k',2);
set(gca,'XTick',tick_marks(1:numel(subjects))+(sessions/2), 'YTick', tick_marks(1:numel(subjects))+(sessions/2), 'XTickLabel',...
    subjects, 'YTickLabel', subjects);
axis square;
colorbar;
saveas(gcf,['ComparisonFCParcels333_beforeandafterPostFreeSurferUpdate_' datasets{d} '.jpg'],'jpg');
close('all');
end