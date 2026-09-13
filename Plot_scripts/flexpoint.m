function [t_flex,tLR,dt] = flexpoint(s,l,m)
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    if abs(s.dyn.th0-pi/2)>1e-5
        wf = struct;
        wf = DB_mode_rotate_timedep(s,s.dyn,l,m,wf);
        tilde_flag = '\tilde';
    else
        wf = s;
        tilde_flag = '';
    end
    dyn = wf.dyn;
    t = wf.ell(l).emm(m+1).t;
    h = wf.ell(l).emm(m+1).hlm;
    f = freq(h,t);
    dynlabel = KerrLabel(dyn);

    tLR = tLR_splined(dyn);

    tnew = t(1):2:t(end);
    tnew = tnew';
    %tnew = t;
    sf = spline(t,f,tnew);

    df = DB_D1(sf,tnew,4);
    d2f = DB_D1(DB_D1(sf,tnew,4),tnew,4);

    t_cut = (tLR-10):1e-5:(tLR+10);
    sd2f = spline(tnew,d2f,t_cut);
    idx_flex  = find(sd2f<0,1);
    t_flex = t_cut(idx_flex);
    dt = t_flex-tLR;

    midpoint = min([t_flex,tLR]) + abs(dt/2)-tLR;
    x1_inset = min([midpoint - dt,midpoint + dt]);
    x2_inset = max([midpoint - dt,midpoint + dt]);

    t_start = tLR - 200 - tLR;
    t_end = s.dyn.t(end)-5 - tLR;

    %{
    figure
    ax0 = gca;
    set(ax0,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','box','on');
    hold on
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel('$d\omega/dt $','Interpreter','latex','FontSize',labels_fontsize)
    plot(tnew-tLR,df,'LineWidth',my_linewidth,'Color','b')
    yline(0,'LineStyle',':','Color','k')
    xline(0,'Color',[.7,.7,.7],'LineStyle','--')
    xlim([t_start,t_end])

    fig = figure;

    ax1 = gca;
    set(ax1,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','box','on');
    hold on
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel('$d^2\omega/dt^2 $','Interpreter','latex','FontSize',labels_fontsize)
    plot(tnew-tLR,d2f,'LineWidth',my_linewidth,'Color','b')
    plot(ax1,t_cut-tLR,sd2f,'Color','r')
    yline(0,'LineStyle',':','Color','k')
    xline(0,'Color',[.7,.7,.7],'LineStyle','--')
    xlim([t_start,t_end])

    x1_inset = min([midpoint - dt,midpoint + dt]);
    x2_inset = max([midpoint - dt,midpoint + dt]);
    
    x_length = x2_inset - x1_inset;
    ax2 = axes('Position',[.17 .57 .3 .3]);
    box on
    hold on
    xlim([x1_inset,x2_inset])
    set(ax2,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
    plot(ax2,t_cut-tLR,sd2f,'Color','r');
    plot(tnew-tLR,d2f,'LineWidth',my_linewidth,'Color','b')
    xline(0,'Color',[.7,.7,.7],'LineStyle','--')
    yline(0,'LineStyle',':','Color','k')

    %y_rect = ylim;
    %rectangle('Position', [x1_inset,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax1);
    %}
    %close     
    f2 = figure;
    f2.Position(1) = 200;
    ax3 = gca;
    set(ax3,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','box','on');
    hold on
    xlabel('$u-u_{LR}$','Interpreter','latex','FontSize',labels_fontsize)
    %ylabel('$\omega_{\ell m}$','Interpreter','latex','FontSize',labels_fontsize)    
    legend('BackgroundAlpha',.7,'Location','northwest','Interpreter','latex','Fontsize',legend_fontsize)
    plot(t-tLR,f,'LineWidth',my_linewidth,'Color','b','DisplayName',sprintf('$%s \\omega_{22}$',tilde_flag))
    plot(s.dyn.t-tLR,2*s.dyn.Omg_orb,'LineWidth',my_linewidth,'DisplayName','$2\Omega^{\rm orb}$','Color',MyColors('g1'))
    xline(tLR-tLR,'Color',[.7,.7,.7],'LineStyle','--','DisplayName','$t_{\rm LR}$','HandleVisibility','off')
    xline(t_flex-tLR,'Color',[.7,.7,.7],'LineStyle',':','LineWidth',1.5,'DisplayName',sprintf('$u^{\\rm infl}_{%s\\omega_{22}}$',tilde_flag))
    xlim([t_start,t_end])
    vlims = ylim;
    text(-30,vlims(2)-.02,dynlabel,'FontSize',14,'Interpreter','latex')
    
    %ax4 = axes('Position',[.25 .35 .3 .3]);
    ax4 = axes('Position',[.42 .59 .28 .28]);
    box on
    hold on
    xlim([x1_inset,x2_inset])
    set(ax4,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
    plot(t-tLR,f,'LineWidth',my_linewidth,'Color','b')
    %plot(s.dyn.t-tLR,2*s.dyn.Omg_orb,'LineWidth',my_linewidth,'DisplayName','$2\Omega^{\rm orb}$','Color',MyColors('g1'))
    xline(tLR-tLR,'Color',[.7,.7,.7],'LineStyle','--')
    xline(t_flex-tLR,'Color',[.7,.7,.7],'LineStyle',':','LineWidth',1.5)

    %y_rect = ylim;
    %rectangle('Position', [x1_inset,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax3);


end