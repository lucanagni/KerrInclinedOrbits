function OmegaChange
    % ==========================================================================================================
    % Produce Fig. 6 of Paper I
    % ==========================================================================================================

    my_linewidth = 1.5;
    axes_fontsize = 10;
    legend_fontsize = 9;
    labels_fontsize = 18;

    a = 0.7;
    iota = 75;
    N = 20;
    span = linspace(-0.007,-0.002 ,N);

    inputDB.verbose = 0;
    inputDB.th0 = abs(deg2rad(90-iota));
    inputDB.chi1(3) = a;

    f=figure;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    %ylabel('$\Omega$','FontSize',labels_fontsize,'Interpreter','latex')
    xlabel('$t$','FontSize',labels_fontsize,'Interpreter','latex')
    legend('Interpreter','latex','FontSize',legend_fontsize,'Location','northwest')
    hold on

    C = gray(N+3);
    D = summer(N);
    for i=1:N
        inputDB.r0 = DB_LSSO(a,iota) + 0.5 + span(i);
        DB = DB_class(inputDB);
        if i==1
            tLR = tLR_splined(DB);
        end
        dt = tLR-tLR_splined(DB);

        plot(DB.t+dt,DB.Omg,'LineStyle','-','LineWidth',1,'Color',C(i,:),'DisplayName',sprintf('$\\theta_{\\rm end} = %2.f^\\circ$',rad2deg(DB.th(end))))
        plot(DB.t+dt,DB.Omg_orb,'LineStyle','-','LineWidth',1,'Color',D(N-i+1,:),'DisplayName',sprintf('$\\theta_{\\rm end} = %2.f^\\circ$',rad2deg(DB.th(end))),'HandleVisibility','off')
    end
    xlim([tLR-70,tLR+30])
    text(951,0.13,'$\Omega$','FontSize',14,'Interpreter','latex','Color','k')
    text(960,0.11,'$\Omega_{\rm orb}$','FontSize',14,'Interpreter','latex','Color','k')
    text(940,0.22,'${\tt a07i75}$','FontSize',14,'Interpreter','latex','Color','k')

    f.Position(3:4) = 1.05.*f.Position(3:4);

return


    