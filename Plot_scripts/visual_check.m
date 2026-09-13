function visual_check(spin)
    angles = [90,60,45,30];
    for i=1:length(angles)
        dir = sprintf('~/waveforms/K/plunge/th0_%d/a%s/wf.mat',angles(i),spin);
        if isfile(dir)
            load(dir,'s');
            wfs(i) = s;
        end
    end

    s90 = wfs(1);
    s60 = wfs(2);
    s45 = wfs(3);
    s30 = wfs(4);
    a = s90.dyn.chi1(3);

    f1 = figure;
    ax1 = axes('Parent', f1);

    f2 = figure;
    ax2 = axes('Parent', f2);

    f3 = figure;
    ax3 = axes('Parent', f3);

    f1.Position = [100 500 560 420];
    f2.Position = [700 500 560 420];
    f3.Position = [1300 500 560 420];

    [h2_90,t_90] = DB_hpc(s90,0,0,'ellmax',2);
    [h2_60,t_60] = DB_hpc(s60,0,0,'ellmax',2);
    [h2_45,t_45] = DB_hpc(s45,0,0,'ellmax',2);
    [h2_30,t_30] = DB_hpc(s30,0,0,'ellmax',2);

    h3_90 = DB_hpc(s90,0,0,'ellmax',3);
    h3_60 = DB_hpc(s60,0,0,'ellmax',3);
    h3_45 = DB_hpc(s45,0,0,'ellmax',3);
    h3_30 = DB_hpc(s30,0,0,'ellmax',3);

    h4_90 = DB_hpc(s90,0,0,'ellmax',4);
    h4_60 = DB_hpc(s60,0,0,'ellmax',4);
    h4_45 = DB_hpc(s45,0,0,'ellmax',4);
    h4_30 = DB_hpc(s30,0,0,'ellmax',4);

    t_60 = t_60 - find_deltat(h2_90,t_90,h2_60,t_60);
    t_45 = t_45 - find_deltat(h2_90,t_90,h2_45,t_45);
    t_30 = t_30 - find_deltat(h2_90,t_90,h2_30,t_30);
    [~,t_peak] = Apeak(h2_90,t_90);
    t_in = t_peak - 400;
    t_end = t_peak + 200;

    figure(f1)
    plot(ax1,t_90,freq(h2_90,t_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 2$',a))
    hold on
    plot(ax1,t_60,freq(h2_60,t_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 2$',a))
    plot(ax1,t_45,freq(h2_45,t_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 2$',a))
    plot(ax1,t_30,freq(h2_30,t_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 2$',a))

    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end    ])
    ylim([0,0.7])
    hold off

    figure(f2)
    plot(ax2,t_90,freq(h3_90,t_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 3$',a))
    hold on
    plot(ax2,t_60,freq(h3_60,t_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 3$',a))
    plot(ax2,t_45,freq(h3_45,t_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 3$',a))
    plot(ax2,t_30,freq(h3_30,t_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 3$',a))

    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end    ])
    ylim([0,0.7])


    figure(f3)
    plot(ax3,t_90,freq(h4_90,t_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 4$',a))
    hold on
    plot(ax3,t_60,freq(h4_60,t_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 4$',a))
    plot(ax3,t_45,freq(h4_45,t_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 4$',a))
    plot(ax3,t_30,freq(h4_30,t_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 4$',a))

    set(ax3,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end])
    ylim([0,0.7])   


    f4 = figure;
    ax4 = axes('Parent', f4);

    f5 = figure;
    ax5 = axes('Parent', f5);

    f6 = figure;
    ax6 = axes('Parent', f6);

    f4.Position = [100 50 560 420];
    f5.Position = [700 50 560 420];
    f6.Position = [1300 50 560 420];

    figure(f4)
    plot(ax4,t_90,abs(h2_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 2$',a))
    hold on
    plot(ax4,t_60,abs(h2_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 2$',a))
    plot(ax4,t_45,abs(h2_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 2$',a))
    plot(ax4,t_30,abs(h2_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 2$',a))

    set(ax4,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end    ])
    hold off

    figure(f5)
    plot(ax5,t_90,abs(h3_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 3$',a))
    hold on
    plot(ax5,t_60,abs(h3_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 3$',a))
    plot(ax5,t_45,abs(h3_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 3$',a))
    plot(ax5,t_30,abs(h3_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 3$',a))

    set(ax5,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end    ])

    figure(f6)
    plot(ax6,t_90,abs(h4_90),'LineWidth',1,'DisplayName',sprintf('$\\iota = 0, \\ a = %.1f, \\ell = 4$',a))
    hold on
    plot(ax6,t_60,abs(h4_60),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/6, \\ a = %.1f, \\ell = 4$',a))
    plot(ax6,t_45,abs(h4_45),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/4, \\ a = %.1f, \\ell = 4$',a))
    plot(ax6,t_30,abs(h4_30),'LineWidth',1,'DisplayName',sprintf('$\\iota = \\pi/3, \\ a = %.1f, \\ell = 4$',a))

    set(ax6,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u/M$','FontSize',14,'Interpreter','Latex')
    ylabel('$\omega$','Interpreter','Latex','Fontsize',14)
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([t_in,t_end])







return