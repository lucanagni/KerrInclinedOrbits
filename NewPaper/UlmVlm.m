function UlmVlm(m)
    l=2;
    load("~/waveforms/K/plunge/th0_30/a05/wf.mat")
    DB = s.dyn;
    tLR = tLR_splined(DB);
    flags.r=1;
    flags.Newt_switch = 1;
    flags.source = 0;
    flags.analytical_label = 'N';
    [Ulm,Vlm] = DB_Multipoles(DB,l,m,flags);

    w_num = s.ell(l).emm(m+1).hlm;
    T = s.ell(l).emm(m+1).t;
    tau = DB.t;

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    f = figure;
    f.Position(3:4) = [f.Position(3)*2,f.Position(4)*.75];
    tl = tiledlayout(3,1,'Padding','tight','TileSpacing','tight');
    t_in = tLR-350;
    t_end = tLR+15;
    label = KerrLabel(DB);
    wf_label = sprintf('h_{%d%d}',l,m);

    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize-2,'BackgroundAlpha',0.95,'NumColumns',2)%,'Orientation','horizontal')
    xlim([t_in,t_end])
    hold on

    %plot(T,real(w_num),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s]/\\nu$',wf_label),'Color',[0 0 0 0.2],'LineStyle',':')
    plot(tau,(1/sqrt(2)).*abs(Ulm),'LineWidth',my_linewidth,'DisplayName',sprintf('$(N,0)$'),'LineStyle','-','Color',[0,0,1,0.15])
    plot(tau,(1/sqrt(2)).*abs(Vlm),'LineWidth',my_linewidth,'DisplayName',sprintf('$(N,1)$'),'LineStyle','-','Color',[1,0,0,0.15])
    plot(T,abs(w_num),'LineWidth',my_linewidth+0.5,'DisplayName',sprintf('Numerical'),'Color',[0 0 0 .8])
    plot(tau,(1/sqrt(2)).*abs(Ulm - 1i*Vlm),'LineWidth',my_linewidth+1,'DisplayName',sprintf('Analytical (Newt)'),'LineStyle','-.','Color',[1 0 0],'LineWidth',1)
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    ylabel(sprintf('$|%s|/\\nu$',wf_label),'Interpreter','latex','Fontsize',labels_fontsize)

    ys = ylim;
    xs = xlim;

    text(xs(2)-abs(xs(2)-xs(1))./4,ys(2) - abs(ys(2)-ys(1))./10,label,'Interpreter','latex','FontSize',16)

    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')

return
