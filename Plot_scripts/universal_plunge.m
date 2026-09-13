function universal_plunge

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;
    C = [MyColors('a1');MyColors('a2');MyColors('a3');MyColors('a4');MyColors('a5')];
    %C = [MyColors('a5');MyColors('a4');MyColors('a3');MyColors('a2');MyColors('a1')];

    inputDB.th0 = pi/4;
    inputDB.verbose = 0;
    inputDB.chi1(3) = 0.5;
    inputDB.r0 = 5.5;
    inputDB.geodesics = 0;
    inputDB.PA = 1;

    q = 1e3;
    nu = DB_nuX1X2(q);

    fprintf('q = %d\n',q)
    inputDB.q = q;
    DB = DB_class(inputDB);
    tLR_ref = tLR_splined(DB)*0;
    %tLSSO = tLSSO_splined(DB);
    r_ref = DB.r;
    t_ref = DB.Omg;

    f = figure;
    ax = gca;
    set(ax,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$t-t_{\rm LR}$','FontSize',labels_fontsize,'Interpreter','latex')
    ylabel('$r$','FontSize',labels_fontsize,'Interpreter','latex')
    hold on
    plot(t_ref-tLR_ref,r_ref,'LineWidth',my_linewidth,'Color',C(1,:),'DisplayName',sprintf('$\\nu = 10^{-%d}$',log10(q)))
    legend('Location','southwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    %xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    %xline(tLSSO-tLR_ref,'Color',[.7 .7 .7],'LineStyle',':','LineWidth',my_linewidth+0.5,'HandleVisibility','off')
    yline(DB_LSSO(0.5,45),'HandleVisibility','off','Color',[.7 .7 .7],'LineStyle','-.')
    %xlim([tLSSO-tLR_ref-100,35])
    ylim([1.75,5])

    f.Units = 'centimeters';
    f.Position(3) = 21;

    %{
    ax1 = axes('Position',[0.35,0.22,0.35,0.35]);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','LineWidth',0.5,'box','on')
    hold on
    plot(t_ref-tLR_ref,r_ref,'LineWidth',my_linewidth,'Color',C(1,:),'DisplayName',sprintf('$\\nu = 10^{-%d}$',log10(q)))
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xlim([-75,25])
    %}

    for i=2:4
        clear DB
        q = q*10;
        nu = DB_nuX1X2(q);
        inputDB.q = q;
        inputDB.r0 = inputDB.r0-0.25;
        fprintf('q = %d\n',q)
        DB = DB_class(inputDB);
        tLR = tLR_splined(DB)*0;
        r = DB.r;
        t = DB.Omg;
        %dt = tLR_ref - tLR;
        %r2 = spline(t+dt,r,t_ref);
        p = plot(ax,t-tLR,r,'LineWidth',my_linewidth,'Color',C(i,:),'DisplayName',sprintf('$\\nu = 10^{-%d}$',log10(q)));
        uistack(p,'bottom',i-1)
        %r = plot(ax1,t-tLR,r,'LineWidth',my_linewidth/2,'Color',C(i,:),'DisplayName',sprintf('$\\nu = 10^{-%d}$',log10(q)));
        %uistack(r,'bottom',i-1)
        
    end
    uistack(ax1,'top')



    





return