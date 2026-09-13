function [DB,KOS] = DBvKOS(varargin)

    if isempty(varargin)
        inputDB.r0          = 5;
        inputDB.chi1        = [0 0 0.5000];
        inputDB.th0         = pi/2;
        inputDB.Tmax        = 1e5;
        inputDB.geodesics   = 1;

        inputKOS.Tmax = 1e5;
        inputKOS.a = 0.5;
        inputKOS.r0 = 5;
        inputKOS.abstol = 1.0000e-13;
        inputKOS.reltol = 1.0000e-13;

        DB = DB_class(inputDB);
        [~,KOS] = KerrOrbitSolver(inputKOS);
    else
        DB = varargin{1};
        KOS = varargin{2};
    end

    figure
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact')
    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    ylabel('$\Delta r_0$','Interpreter','Latex','Fontsize',18)
    hold on
    plot(DB.t,(DB.r-5)./5,'DisplayName','{\tt DBSolver}','LineWidth',1,'Color','k')
    plot(KOS.Dyn.time,(KOS.Dyn.r-5)./5,'DisplayName','{\tt BLSolver}','LineWidth',1,'Color','r')
    legend('Location','southwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    ylim([-5e-10,2e-10])
    grid on

    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    ylabel('$\Delta E_0$','Interpreter','Latex','Fontsize',18)
    xlabel('$t$','Interpreter','Latex','Fontsize',18)
    hold on
    plot(DB.t,(DB.Heff-DB.Heff(1))./DB.Heff(1),'DisplayName','{\tt DBSolver}','LineWidth',1,'Color','k')
    plot(KOS.Dyn.time,(KOS.Dyn.H-KOS.Dyn.H(1))./KOS.Dyn.H(1),'DisplayName','{\tt BLSolver}','LineWidth',1,'Color','r')
    legend('Location','southwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    ylim([-14e-12,  3e-12])
    grid on

    inputKOS.a = 0.5;
    inputKOS.abstol = 1.0000e-13;
    inputKOS.reltol = 1.0000e-13;
    inputKOS.Tmax = 1e4;
    inputKOS.sr0 = 8;
    inputKOS.e0 = 0.5;
    [~,KOSecc] = KerrOrbitSolver(inputKOS);

    pphi = KOSecc.Dyn.pph(1);
    r0 = KOSecc.Dyn.r(1);

    inputDB.Tmax        = 1e4;
    inputDB.DBvKOS      = pphi;
    inputDB.r0          = r0;
    inputDB.chi1        = [0 0 0.5000];
    inputDB.th0         = pi/2;
    inputDB.geodesics   = 1;
    DBecc = DB_class(inputDB);

    figure
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact')
    nexttile(1,[2,1])
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    ylabel('$r$','Interpreter','Latex','Fontsize',18)
    hold on
    plot(DBecc.t,DBecc.r,'DisplayName','{\tt DBSolver}','LineWidth',1,'Color','k')
    plot(KOSecc.Dyn.time,KOSecc.Dyn.r,'DisplayName','{\tt BLSolver}','LineWidth',1,'LineStyle','--','Color','r')
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.9)
    ylim([4,17])

    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    ylabel('$|\Delta r/r|$','Interpreter','Latex','Fontsize',18)
    xlabel('$t$','Interpreter','Latex','Fontsize',18)
    hold on
    plot(DBecc.t,abs(DBecc.r-spline(KOSecc.Dyn.time,KOSecc.Dyn.r,DBecc.t))./DBecc.r,'LineWidth',1,'Color','k')
    grid on
return