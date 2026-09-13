function TEMP_PAplot
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;


    a1 = 0.5;
    a2 = 0.99;

    inputDB.verbose = 0;

    inputDB.chi1(3) = a1;
    iota1 = 75;
    inputDB.r0 = DB_LSSO(a1,iota1) + 1;
    inputDB.th0 = deg2rad(90-iota1);
    inputDB.PA = 0;
    inputDB.geodesics = 0;
    a05PA0 = DB_class(inputDB);

    inputDB.PA = 1;
    a05PA1 = DB_class(inputDB);

    inputDB.chi1(3) = a2;
    iota2 = 75;
    inputDB.r0 = DB_LSSO(a2,iota2) + 1;
    inputDB.th0 = deg2rad(90-iota2);
    inputDB.PA = 0;
    a09PA0 = DB_class(inputDB);

    inputDB.PA = 1;
    a09PA1 = DB_class(inputDB);

    f1 = figure;
    ax0 = gca;
    set(ax0,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(a05PA0.t,a05PA0.r,'LineWidth',my_linewidth,'DisplayName','Circular')
    plot(a05PA1.t,a05PA1.r,'LineWidth',my_linewidth,'DisplayName','Quasi-Circular')
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.6)
    xlabel('$t$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$r$','FontSize',labels_fontsize,'Interpreter','Latex');
    hold off

    ax1 = axes('Position',[.25,.22,.3,.3]);
    set(ax1,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
    hold on
    plot(a05PA0.t,a05PA0.r,'LineWidth',my_linewidth,'DisplayName','Circular')
    plot(a05PA1.t,a05PA1.r,'LineWidth',my_linewidth,'DisplayName','Quasi-Circular')    
    xlim([300,1000])
    y_rect = ylim;
    x_rect = xlim;
    
    %rectangle('Position', [x_rect(1),y_rect(1),abs(x_rect(2) - x_rect(1)),abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax0);
    text(ax0,x_rect(1),y_rect(1) - 0.5,sprintf('$a = %.2f,\\ \\iota = %d$',a1,iota1),'Interpreter','latex','FontSize',14)


    f2 = figure;
    ax2 = gca;
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(a09PA0.t,a09PA0.r,'LineWidth',my_linewidth,'DisplayName','Circular')
    plot(a09PA1.t,a09PA1.r,'LineWidth',my_linewidth,'DisplayName','Quasi-Circular')
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.6)
    xlabel('$t$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$r$','FontSize',labels_fontsize,'Interpreter','Latex');
    hold off

    ax3 = axes('Position',[.25,.22,.47,.47]);
    set(ax3,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
    hold on
    plot(a09PA0.t,a09PA0.r,'LineWidth',my_linewidth,'DisplayName','Circular')
    plot(a09PA1.t,a09PA1.r,'LineWidth',my_linewidth,'DisplayName','Quasi-Circular')    
    xlim([300,1000])
    y_rect = ylim;
    x_rect = xlim;
    
    %rectangle('Position', [x_rect(1),y_rect(1),abs(x_rect(2) - x_rect(1)),abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax2);
    %text(ax2,x_rect(1),y_rect(1) - 0.5,sprintf('$a = %.2f,\\ \\iota = %s$',a2,'5\\frac{5\\pi}{12}'),'Interpreter','latex','FontSize',14)
    text(ax2,157,4.7,sprintf('$a = %.2f,\\ \\iota = %s$',a2,'\frac{5\pi}{12}'),'Interpreter','latex','FontSize',14)



    A4Width(f2);