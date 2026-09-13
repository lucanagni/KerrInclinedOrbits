function TwoTraj
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    inputDB.verbose = 0;
    inputDB.r0 = 7;
    inputDB.chi1(3) = 0;
    inputDB.geodesics = 0;
    inputDB.th0 = pi/2;
    DBeq = DB_class(inputDB);

    inputDB.th0 = pi/3;
    DBp = DB_class(inputDB);

    RGB = orderedcolors("gem");
    C1 = RGB(1,:);
    C2 = RGB(2,:);

    lims = [-7 7];
    zlims = [-4,4];
    figure
    plot3(DBeq.x,DBeq.y,DBeq.z,'LineWidth',my_linewidth,'Color',C1)
    set(gca,'LineWidth',1,'FontName','Times','FontSize',axes_fontsize)
    x1 = xlabel('$x$','FontSize',labels_fontsize,'Interpreter','latex');
    y1 = ylabel('$y$','FontSize',labels_fontsize,'Interpreter','latex');
    z1 = zlabel('$z$','FontSize',labels_fontsize,'Interpreter','latex');
    xlim(lims)
    ylim(lims)
    zlim(zlims)
    Lx = abs(lims(2)+lims(1))/2;
    Ly = abs(lims(2)+lims(1))/2;
    Lz = abs(zlims(2)+zlims(1))/2;
    x1.Position(1:3) = [Lx,lims(1)-2,zlims(1)];
    y1.Position(1:3) = [lims(1)-2,Ly,zlims(1)];
    z1.Position(1:3) = [lims(1),lims(2)+2,Lz];
    grid on

    figure
    plot3(DBp.x,DBp.y,DBp.z,'LineWidth',my_linewidth,'Color',C2)
    set(gca,'LineWidth',1,'FontName','Times','FontSize',axes_fontsize)
    x2 = xlabel('$x$','FontSize',labels_fontsize,'Interpreter','latex');
    y2 = ylabel('$y$','FontSize',labels_fontsize,'Interpreter','latex');
    z2 = zlabel('$z$','FontSize',labels_fontsize,'Interpreter','latex');
    xlim(lims)
    ylim(lims)
    zlim(zlims)
    Lx = abs(lims(2)+lims(1))/2;
    Ly = abs(lims(2)+lims(1))/2;
    Lz = abs(zlims(2)+zlims(1))/2;
    x1.Position(1:3) = [Lx,lims(1)-2,zlims(1)];
    y1.Position(1:3) = [lims(1)-2,Ly,zlims(1)];
    z1.Position(1:3) = [lims(1),lims(2)+2,Lz];
    grid on
return
