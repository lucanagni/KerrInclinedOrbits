function catalog_figs
    % ==========================================================================================================
    % Automatically generate and save figures for the appendix of PaperI
    % ==========================================================================================================

    warning('off','MATLAB:print:ContentTypeImageSuggested')
    iotas = [0 30 45 60 120 135 150 180];
    spins = [2 5 9];
    for i=1:length(spins)
        for j=1:length(iotas)
            mydir = finddir(spins(i),iotas(j));
            mylabel = dir2label(spins(i),iotas(j));
            fprintf('dir:%s, label: %s\n',mydir,mylabel)
            
            dir = sprintf('~/waveforms/K/plunge/%s/wf.mat',mydir);
            load(dir)

            [af,dyn,hp] = hlm_plots(s,2,2,'emminus','hpc');
            af_label = sprintf('catalog/%s_af.pdf',mylabel);
            dyn_label = sprintf('catalog/%s_dyn.pdf',mylabel);
            hp_label = sprintf('catalog/%s_hp.pdf',mylabel);
            exportgraphics(af, af_label, 'ContentType', 'vector')
            exportgraphics(dyn, dyn_label, 'ContentType', 'vector')
            exportgraphics(hp, hp_label, 'ContentType', 'vector')

            af_fig = sprintf('catalog/figs/%s_af.fig',mylabel);
            dyn_fig = sprintf('catalog/figs/%s_dyn.fig',mylabel);
            hp_fig = sprintf('catalog/figs/%s_hp.fig',mylabel);
            saveas(af, af_fig)
            saveas(dyn, dyn_fig)
            saveas(hp, hp_fig)
            close all
            
        end
    end


    

return

function dir = finddir(a,iota)
    if iota>90
        sign = '-';
    else
        sign = '';
    end
    th = abs(90-iota);
    th0 = num2str(th);

    dir = sprintf('th0_%s/a%s0%d',th0,sign,a);
return

function label = dir2label(a,iota)


    label = sprintf('a0%di%d',a,iota);

return

%{
function label = dir2label(dir)
    th0 = dir(strfind(dir,'th0_')+4:strfind(dir,'th0_')+5);
    if contains(dir,'-')
        sign = 'm';
        spin = dir(strfind(dir,'-')+1:strfind(dir,'-')+2);
        iota = num2str(90+str2double(th0));
    else
        sign = '';
        spin = dir(strfind(dir,'a')+1:strfind(dir,'a')+2);
        iota = num2str(90-str2double(th0));
    end

    label = sprintf('a%s%si%s',sign,spin,iota);

return
%}