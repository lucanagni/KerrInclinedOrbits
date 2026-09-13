function S_rtp(s)
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    dyn = s.dyn;
    t = dyn.t;
    r = dyn.r;
    theta = dyn.th;
    phi = dyn.phi;

    tLR = tLR_splined(dyn);
    t_peak = tLR;
    t_in = t_peak-230;
    t_end = t_peak + 30;

    f1 = figure;
    tl = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');

%{
    ax1 = nexttile(tl,1);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(t,r,'LineWidth',my_linewidth,'Color',[0 0 0 0.9])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    ylabel('$r$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([t_in,t_end])
    %}

    ax2 = nexttile(tl,1);
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(t,phi,'LineWidth',my_linewidth,'Color',[0 0 0 0.9])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    ylabel('$\varphi$','FontSize',labels_fontsize,'Interpreter','latex')
    yticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
    yticklabels({'0','2\pi','4\pi','6\pi','8\pi','10\pi'})
    %ylim([pi,10*pi])
    xlim([t_in,t_end])

    ax3 = nexttile(tl,2);
    set(ax3,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(t,theta,'LineWidth',my_linewidth,'Color',[0 0 0 0.9])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    yline(pi/2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',1,'HandleVisibility','off')
    ylim([pi/3-0.1,2*pi/3+0.1])
    yticks([pi/3 pi/2 2*pi/3])
    yticklabels({'\pi/3','\pi/2','2\pi/3'})
    ylabel('$\theta$','FontSize',labels_fontsize,'Interpreter','latex')
    xlabel('$t$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([t_in,t_end])

    
    %load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat');
    %dyn_eq = s.dyn;

    %f2 = figure;
    %plot3(dyn.x,dyn.y,dyn.z,'LineWidth',my_linewidth,)
    %set(ax3,'XMinorTick','on','box','off','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');


return 