function [rel_diff,t1] = dJdt_reldiff
    % ========================================================================================
    % for PDF fontsize 14 should be fine
    % input are [h,t] for equatorial and precessing dynamics (in theory at least)
    % ========================================================================================
    my_linewidth = 2;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    load('~/waveforms/S/fluxes/th0_90/wf.mat');
    s90 = s;
    load('~/waveforms/S/fluxes/th0_60/wf.mat');
    s60 = s;

    tLR1 = tLR_splined(s90.dyn);
    tLR2 = tLR_splined(s60.dyn);
    dt = tLR1 - tLR2;

    [J90v,t] = DB_dJdt(s90,'ellmax',3);
    [J60v,t60] = DB_dJdt(s60,'ellmax',3);

    J90 = sqrt(J90v(:,1).^2 + J90v(:,2).^2 + J90v(:,3).^2);
    J60 = sqrt(J60v(:,1).^2 + J60v(:,2).^2 + J60v(:,3).^2);

    J60 = spline(t60+dt,J60,t);

    t_start = tLR1-200;
    t_end = tLR1+50;

    f1 = figure;
    f1.Position(3:4) = 1.1.*f1.Position(3:4);

    ax1 = subplot(3,1,[1,2]);
    hold on
    plot(t,J90,'LineWidth',my_linewidth,'DisplayName','$\iota = 0$','LineStyle','-','Color',MyColors('b1'))
    plot(t,J60,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/6$','LineStyle','--','Color',MyColors('r1'))
    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    ylabel('$|\dot{\mathbf{J}}|$', 'Interpreter', 'latex','FontSize',labels_fontsize)
    ylim([0,0.075])
    
    subplot(3,1,3)
    plot(t,abs(J90-J60)./J90,'LineWidth',my_linewidth-1,'Color','k')
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    xlim([t_start,t_end])
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel('$ | \Delta \dot{\mathbf{J}} |/ | \dot{\mathbf{J}_{eq}} |$', 'Interpreter', 'latex','FontSize',labels_fontsize)

    ymax = ylim;
    ylim([0,1.1*ymax(2)])

    grid on

    f2 = figure;
    f2.Position(1:4) = [f1.Position(1)+f1.Position(3)+50,f1.Position(2),f1.Position(3),f1.Position(4)];
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    ylabel('$|\dot{{J_i}}|$', 'Interpreter', 'latex','FontSize',labels_fontsize)
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    hold on
    plot(t60+dt,J60v(:,3),'LineWidth',my_linewidth,'Color',MyColors('r2'),'DisplayName','$\dot{J_z}$')
    plot(t60+dt,J60v(:,2),'LineWidth',my_linewidth,'Color',MyColors('r3'),'DisplayName','$\dot{J_y}$')
    plot(t60+dt,J60v(:,1),'LineWidth',my_linewidth,'Color',MyColors('r1'),'DisplayName','$\dot{J_x}$')
    xlim([t_start,t_end])


    ymax = ylim;
    ylim([ymax(1),1.1*ymax(2)])

return