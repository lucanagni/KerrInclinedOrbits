function S_af_reldiffs(l,m)
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize;
    labels_fontsize = 14;

    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat');
    wf_eq = s;
    load('~/waveforms/S/plunge/th0_60/wf.mat');
    wf_prec = s;
    wf_rot = DB_mode_rotate_timedep(wf_prec,wf_prec.dyn,l,m);

    heq = wf_eq.ell(l).emm(m+1).hlm./norm;
    teq = wf_eq.ell(l).emm(m+1).t;
    
    hr = wf_rot.ell(l).emm(m+1).hlm./norm;
    tr = wf_rot.ell(l).emm(m+1).t;

    dt = find_deltat(heq,teq,hr,tr);
    %dt = 0;
    tr = tr - dt;

    heq = spline(teq,heq,tr);
    t = tr;

    [peak,t_peak] = Apeak(hr,t);

    t_start = t_peak-250;
    t_end = t_peak+100;

    f = figure;
    tl = tiledlayout(3,2,'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl,1,[2,1]);  
    set(ax,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(ax,t,abs(heq),'LineWidth',my_linewidth,'DisplayName','$\iota = 0^\circ$','LineStyle','-','Color','r')
    plot(ax,t,abs(hr),'LineWidth',my_linewidth,'DisplayName','$\iota = 30^\circ$ (rotated)','LineStyle','--','Color','b')
    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    if m==1
        legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','southeast')
    end
    ylabel(sprintf('$|\\Psi_{%d%d}|/\\nu$',l,m), 'Interpreter', 'latex')

    x1_rect = t_start + 100;
    x2_rect = t_start + 150;
    x_length = x2_rect - x1_rect;

    axPos = get(ax, 'Position');  % [x y w h]
    insetWidth = 0.4 * axPos(3);  
    insetHeight = 0.4 * axPos(4); 
    insetX = axPos(1) + axPos(3)./10; 
    insetY = axPos(2) + axPos(4)./6;
    if m==1
        insetX = axPos(1) + axPos(3)./15; 
        insetY = axPos(2) + axPos(3).*3./4;
    end

    axInset = axes('Position', [insetX insetY insetWidth insetHeight]);
    set(axInset,'XMinorTick','on','YMinorTick','on','box','on','FontSize',round(axes_fontsize/2)+1,'FontName','Times');
    hold on
    plot(axInset,t,abs(heq),'LineWidth',my_linewidth,'DisplayName','$\iota = 0^\circ$','LineStyle','-','Color','r')
    plot(axInset,t,abs(hr),'LineWidth',my_linewidth,'DisplayName','$\iota = 30^\circ$','LineStyle','--','Color','b')
    xlim([x1_rect,x2_rect])
    ylim auto

    y_rect = ylim;
    rectangle('Position', [x1_rect,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent', ax);

    nexttile(tl,5)
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,abs(abs(heq)-abs(hr))./abs(heq),'LineWidth',my_linewidth,'Color','k')
    xlim([t_start,t_end])
    ylabel('$| \Delta A |/ A_{eq}$', 'Interpreter', 'latex')
    xlabel('$u$','Interpreter','latex')

    ymax = ylim;
    ylim([0,1.1*ymax(2)])

    grid on


    axr = nexttile(tl,2,[2,1]);  
    set(axr,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(axr,t,freq(heq,t),'LineWidth',my_linewidth,'DisplayName','$\iota = 0^\circ$','LineStyle','-','Color','r')
    plot(axr,t,freq(hr,t),'LineWidth',my_linewidth,'DisplayName','$\iota = 30^\circ$ (rotated)','LineStyle','--','Color','b')
    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    ylabel(sprintf('$\\omega_{%d%d}$',l,m), 'Interpreter', 'latex')

    axrPos = get(axr, 'Position');  % [x y w h]
    insetWidth = 0.4 * axrPos(3);  
    insetHeight = 0.4 * axrPos(4); 
    insetX = axrPos(1) + axrPos(3)./10; 
    insetY = axrPos(2) + 3*axrPos(4)./10;

    axrInset = axes('Position', [insetX insetY insetWidth insetHeight]);
    set(axrInset,'XMinorTick','on','YMinorTick','on','box','on','FontSize',round(axes_fontsize/2)+1,'FontName','Times');
    hold on
    plot(axrInset,t,freq(heq,t),'LineWidth',my_linewidth,'DisplayName','$\iota = 0^\circ$','LineStyle','-','Color','r')
    plot(axrInset,t,freq(hr,t),'LineWidth',my_linewidth,'DisplayName','$\iota = 30^\circ$ (rotated)','LineStyle','--','Color','b')
    xlim([x1_rect,x2_rect])
    ylim auto

    y_rect = ylim;
    rectangle('Position', [x1_rect,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent', axr);

    nexttile(tl,6)
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,abs(freq(heq,t)-freq(hr,t))./freq(heq,t),'LineWidth',my_linewidth,'Color','k')
    xlim([t_start,t_end])
    ylabel('$| \Delta \omega |/  \omega_{eq} $', 'Interpreter', 'latex')
    xlabel('$u$','Interpreter','latex')

    ymax = ylim;
    ylim([0,1.1*ymax(2)])

    grid on

    A4Width(f);
return
