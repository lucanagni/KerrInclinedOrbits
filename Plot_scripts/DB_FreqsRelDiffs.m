function DB_FreqsRelDiffs(l,m)
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    load('~/waveforms/K/plunge/th0_90/a05/wf.mat')
    se = s;
    eq_label = KerrLabel(s.dyn);
    load('~/waveforms/K/plunge/th0_45/a05/wf.mat')
    sp = s;
    prec_label = KerrLabel(s.dyn);

    %f1 = freq(se.ell(l).emm(m+1).hlm,se.ell(l).emm(m+1).t);
    f1 = freq(sp.ell(l).emm(m+1).hlm,sp.ell(l).emm(m+1).t);
    f2 = freq((-1).^l.*conj(sp.ell(l).emminus(m+1).hlm),sp.ell(l).emminus(m+1).t);
    t1 = sp.ell(l).emm(m+1).t;
    t2 = sp.ell(l).emminus(m+1).t;

    f2 = spline(t2,f2,t1);

    tLR = tLR_splined(se.dyn);
    %tLR2 = tLR_splined(sp.dyn);

    %dt = tLR2 - tLR1;
    
    t_start = tLR - 500;
    t_end = tLR + 210;
    figure
    ax1 = subplot(3,1,[1,2]);
    hold on
    plot(t1,f1,'LineWidth',my_linewidth,'DisplayName','$h_{\l m}$','LineStyle','-','Color',MyColors('b1'))
    plot(t1,f2,'LineWidth',my_linewidth,'DisplayName','$(-)^\ell h_{\ell -m}^*$','LineStyle','-','Color',MyColors('r1'))
    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    ylabel(sprintf('$\\omega_{%d%d}$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)
    xlims = xlim;
    xlength = xlims(2) - xlims(1);
    ylims = ylim;
    ylength = ylims(2) - ylims(1);
    text(gca,xlims(1) + 4*xlength/9,ylims(2) - ylength/10,prec_label,'FontSize',14,'Interpreter','latex')
    
    subplot(3,1,3)
    plot(t1,abs(f2-f1),'LineWidth',my_linewidth,'Color','k')
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    xlim([t_start,t_end])
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel(sprintf('$\\Delta \\omega_{%d%d}$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)
    ylims = ylim;
    ylim(1.1*[ylims(1),ylims(2)]);     

    grid on
