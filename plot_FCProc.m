datadir = '/projects/b1081/Lifespan/Nifti/derivatives/preproc_FCProc';
outdir = '/projects/p31161/Lifespan';
mkdir(outdir)
sub = 'LS05';
ses = [1:5];
%runs = [9,9,11,8,9]; %LS03
runs = [8, 8, 8, 9, 9]; %LS05
%runs = [7, 8]; %LS07

for s = 1:length(ses)    
    for r = 1:runs(s)
        filename = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%d_fmriprep_zmdt_resid_ntrpl_bpss_zmdt.nii.gz',datadir,sub,ses(s),sub,ses(s),r);
        %filename = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz',datadir,sub,ses(s),sub,ses(s),runs(r));
        nifti = load_nii_wrapper(filename);
        run =[];
        for v = [1:1000:902629]
          run = [run nifti(v,:)'];
            %pause(2);
        end
        %addptitle = sprintf('%s/sub_%s_ses-%d_run%02d.png',outdir, sub, ses(s), r);
        %plot(run)
        plot(run(:,45))
        hold on
        %title(gcf, ['sub' sub ', ses ' num2str(ses(s)) ', run ' num2str(r)]);
        %saveas(gcf,title)
    end
end
        