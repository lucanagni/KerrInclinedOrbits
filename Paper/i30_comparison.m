function i30_comparison(save)
    % ==========================================================================================================
    % Produce Fig. 15 of Paper I
    % Assumes Schwarzschild waveforms is located in ~/waveforms/S/fluxes/th0_60/wf.mat
    % Assumes Kerr waveforms are located in ~/waveforms/K/plunge/th0_60/a0<value>/wf.mat
    % wf.mat is a struct containing the dynamics (element of DB_class) as wf.dyn and modes as wf.ell(l).emm(m+1).hlm (emminus for negative m)
    % ==========================================================================================================

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    load('~/waveforms/S/fluxes/th0_60/wf.mat');
    a0 = s;
    load('~/waveforms/K/plunge/th0_60/a02/wf.mat');
    a02 = s;
    load('~/waveforms/K/plunge/th0_60/a05/wf.mat');
    a05 = s;
    load('~/waveforms/K/plunge/th0_60/a09/wf.mat');
    a09 = s;
    clear s

    addpath '/home/luca/repos/kerrorbitsolver/KerrDynamics' 
    addpath Kerr_stuff

    dynamics = [a0, a02, a05, a09];
    l=2;
    m=2;
    norm = sqrt((l+2).*(l+1).*l.*(l-1));
        
    for i=1:length(dynamics)
        s = dynamics(i);
        config_label = KerrLabel(s.dyn);

        % =======================
        % QNM frequency
        % =======================
        a = s.dyn.chi1(3);
        if a>0
            sigma = kerr_FitKerrQNMs(l,m,a,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
        else
            sigma = kerr_FitKerrQNMs(l,-m,abs(a),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
        end
        omg_QNM = imag(sigma);

        % =======================
        % Dynamical quantities
        % =======================
        dyn = s.dyn;
        T = dyn.t;
        tLR = tLR_splined(dyn);
        tLSSO = tLSSO_splined(dyn);

        t_in = tLR - 200;
        t_end = tLR + 100;

        % =======================
        % Waveform
        % =======================
        t = s.ell(l).emm(m+1).t;
        hlm = s.ell(l).emm(m+1).hlm;
        hlmminus = s.ell(l).emminus(m+1).hlm;
    
        % =======================
        % Amplitude-Frequency plots
        % =======================
        figure
        tiledlayout(2, 1,'TileSpacing', 'compact', 'Padding', 'compact');
        %f1.Position(1) = f1.Position(1) + 300;
        nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,abs(hlm)./norm,'LineWidth',1,'Color','k','DisplayName',sprintf('$|\\Psi_{%d%d}|/\\nu$',l,m))
        plot(t,abs(hlmminus)./norm,'LineWidth',1.5,'LineStyle',':','Color','r','DisplayName',sprintf('$|\\Psi_{%d-%d}|/\\nu$',l,m))

        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(tLSSO,'Color',[.7 .7 .7],'LineStyle','-.','HandleVisibility','off')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
        xlim([t_in,t_end]) 
        lims = ylim;
        ylim(lims*1.1)
        xlength = t_end-t_in;
        ylength = abs(lims(2)-lims(1))*1.1;
        text(t_end-xlength/4,lims(2)*1.1 - ylength/8,config_label,'FontSize',legend_fontsize,'Interpreter','latex')

        nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        lnk = plot(t,freq(hlm,t),'LineWidth',1,'Color','k','DisplayName',sprintf('$\\omega_{%d%d}$',l,m));
        yline(omg_QNM,'LineStyle','--','HandleVisibility','off')
        plot(t,freq((-1).^l.*conj(hlmminus),t),'LineWidth',1.5,'LineStyle',':','Color','r','DisplayName',sprintf('$-\\omega_{%d-%d}$',l,m))

        lnh = plot(T,2*dyn.Omg_orb,'Color',[MyColors('g1'),0.5],'LineWidth',2,'DisplayName','$2\Omega_{\rm orb}$');
        %lnh1 = plot(T,2*(sin(dyn.th)).*DB_D1(dyn.th,dyn.t,4),'Color',[.8 .8 .8],'DisplayName','$2\Omega$');
        uistack(lnh,"bottom");
        uistack(lnk,'top');
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(tLSSO,'Color',[.7 .7 .7],'LineStyle','-.','HandleVisibility','off')
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        yline(0,'Color',[.9 .9 .9],'LineWidth',.25,'HandleVisibility','off')
        L = legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7);
        L.Direction = 'reverse';
        xlim([t_in,t_end])
        xlength = t_end-t_in;
        lims = ylim;
        if abs(lims(2)-omg_QNM)<0.07
            lims = ylim;
            ylim([lims(1),lims(2)*1.1])
        end
        ylength = abs(lims(2)-lims(1));


        if save==1
            pdf_label = config_label(strfind(config_label,'a'):strfind(config_label,'}')-1);
            dir = sprintf('/home/luca/repos/teobiresumsprecessing/Latex/Paper/Figs/%s_comparison.pdf',pdf_label);
            %exportgraphics(gcf, dir, 'ContentType', 'vector')
            saveas(gcf,dir)
            sprintf(['pdfcrop ', dir, ' ', dir])
            system(['pdfcrop ', dir, ' ', dir]);
        end




    end