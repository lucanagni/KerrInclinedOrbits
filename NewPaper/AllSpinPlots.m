function AllSpinPlots()
    warning('off','MATLAB:print:ContentTypeImageSuggested')

    % ==================== PLOTTING PARAMETERS ====================
    my_linewidth = 2.0;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    spins = [ 02 05 09];
    iotas = [0 30 45 60];
    outdir = '/home/luca/repos/teobiresumsprecessing/Latex/NewPaper/Notes/Figs/catalogue';

for j=1:length(spins)
    a = spins(j);
    for k=1:length(iotas)
        i = iotas(k);
        for m=1:2
            DB_waveforms_test(a,i,2,m,'rhoresum','hybrid','deltaresum','hybrid','4p5')
            exportgraphics(gcf, sprintf('%s/a0%di%d_h2%d.pdf',outdir,a,i,m), 'ContentType', 'vector')
            close

            if a==0
                DB_waveforms_test(a,i,2,m,'rhoresum','hybrid','deltaresum','hybrid','4p5','geod')
                exportgraphics(gcf, sprintf('%s/a0%di%d_h2%d_geod.pdf',outdir,a,i,m), 'ContentType', 'vector')
                close
            end
        end
    end
end
return
