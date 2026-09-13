function RecapPlots(spins)
    addpath Plot_scripts
    addpath Paper

    incl = [00,30,45,60];
    %spins = string(spinvec);

    for j=1:length(spins)
        a = spins(j);
        for i=1:length(incl)
            iota = incl(i);
            
            if iota==0
                label = sprintf('a0%si0%d',string(a),iota);
            else
                label = sprintf('a0%si%d',string(a),iota);
            end

            DB_waveforms_test(a,iota,2,2,'recap');
            exportgraphics(gcf, sprintf('~/repos/teobiresumsprecessing/Latex/NewPaper/Figs/IotaExp/%s_h22.pdf',label), 'ContentType', 'vector')
            %exportgraphics(gcf, sprintf('~/Desktop/IotaExp/%s_h22.pdf',label), 'ContentType', 'vector')
            close 

            DB_waveforms_test(a,iota,2,1,'recap');
            exportgraphics(gcf, sprintf('~/repos/teobiresumsprecessing/Latex/NewPaper/Figs/IotaExp/%s_h21.pdf',label), 'ContentType', 'vector')
            %exportgraphics(gcf, sprintf('~/Desktop/IotaExp/%s_h21.pdf',label), 'ContentType', 'vector')
            close 
        end
    end
return