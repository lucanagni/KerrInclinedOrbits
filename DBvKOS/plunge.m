function plunge
    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    load('/home/luca/waveforms/DB_KOS/a09/KOS.mat')
    wf_KOS = s;
    load('/home/luca/waveforms/DB_KOS/a09/DB.mat')
    wf_DB = s;

    h1 = wf_DB.ell(2).emm(3).hlm;
    t1 = wf_DB.ell(2).emm(3).t;
    DB = wf_DB.dyn;

    h2 = wf_KOS.ell(2).emm(3).hlm;
    t2 = wf_KOS.ell(2).emm(3).t;
    KOS = wf_KOS.dyn;

    tLR = tLR_splined(DB);
    tLR2 = KOS_tLR_splined(KOS);
    dt = tLR2-tLR;

    tLSSO = tLSSO_splined(DB);
    tLSSO2 = KOS_tLSSO_splined(KOS);
    %dt = tLSSO2 - tLSSO;

    t_start = tLR-2200;
    t_end = tLR+100;

    figure
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact')

    nexttile(1,[2,1])
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontName','Times','LineWidth',0.5,'FontSize',axes_fontsize)
    ylabel('$|h_{22}|/\nu$','FontSize',labels_fontsize,'Interpreter','latex')
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    hold on
    plot(t1,abs(h1),'LineWidth',my_linewidth,'LineStyle','-','Color',MyColors('b1'),'DisplayName','RR of Ref. [41]')
    plot(t2-dt,abs(h2),'LineWidth',my_linewidth,'LineStyle','--','Color',MyColors('r1'),'DisplayName','EOB-resummed RR')
    xline(tLR,'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')
    xline(tLR2-dt,'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')
    xline(tLSSO,'LineStyle',':','LineWidth',1,'Color',[.7 .7 .7],'HandleVisibility','off')
    xlim([t_start,t_end])
    lims = ylim;
    ylim(1.2*lims)



    nexttile(3)
    grid on
    ylabel('$|\Delta A_{22}|$','FontSize',labels_fontsize,'Interpreter','latex')
    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')
    hold on
    plot(t1,abs(abs(h1)-spline(t2-dt,abs(h2),t1))./abs(h1),'LineWidth',my_linewidth,'Color','k')
    xlim([t_start,t_end])