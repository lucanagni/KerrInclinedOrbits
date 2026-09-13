function delta = DB_compare_freq(dyn,varargin)
    % ==========================================================================================================
    % Compares different definitions of frequency in the precessing case. Not sure this is useful 
    % ==========================================================================================================
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    a = dyn.chi1(3);
    Omg_H = a./(2.*(1+sqrt(1-a.^2)));
    Omg = dyn.Omg;
    Omg_orb = dyn.Omg_orb;
    Omg_so = dyn.Omg_so;
    Omg_ph = dyn.Omg_phi;
    Omg_th = dyn.Omg_th;
    sth = sin(dyn.th);
    t = dyn.t;

    [~,max_orb] = OmgPeak(dyn,'orb');
    [~,max_tot] = OmgPeak(dyn,'tot');
    
    tLR = tLR_splined(dyn);

    %delta = tLR-max_orb;
    %delta = max_tot - max_orb;
    delta = tLR - max_orb;
    %fprintf('tLR - tOrb = %.5f\n',delta)

    if isempty(varargin)
        label = KerrLabel(dyn);
        f = figure;
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,Omg_th,'LineWidth',my_linewidth/2,'LineStyle','-','DisplayName','$\Omega_\theta$','Color',MyColors('b1'))
        %plot(t,sth.*Omg_th,'LineWidth',my_linewidth,'DisplayName','$\Omega_\theta\sin\theta$','Color',MyColors('b2'))
        plot(t,Omg_ph,'LineWidth',my_linewidth/2,'LineStyle','-','DisplayName','$\Omega_\varphi$','Color',MyColors('r1'))
        plot(t,Omg,'LineWidth',my_linewidth,'DisplayName','$\Omega$','Color','k')
        plot(t,Omg_orb,'LineWidth',my_linewidth,'DisplayName','$\Omega^{\rm orb}$','Color',MyColors('g1'))
        plot(t,Omg_so,'LineWidth',my_linewidth,'DisplayName','$\Omega^{\rm so}$','Color',MyColors('a1'))
        xlabel('$t$','FontSize',labels_fontsize,'Interpreter','latex')
        legend('Location','northwest','FontSize',legend_fontsize,'Interpreter','latex')
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(max_orb,'Color',[0 0 0],'LineStyle','-.','LineWidth',.5,'HandleVisibility','off')
        yline(Omg_H,'Color',[0 0 0 .9],'LineStyle','--','HandleVisibility','off')
        yline(0,'Color',[.9 .9 .9],'LineWidth',.25,'HandleVisibility','off')
        xlim([t(end)-200,t(end)+70])

        xlims = xlim;
        ylims = ylim;
        xlength = xlims(2) - xlims(1);
        ylength = ylims(2) - ylims(1);
        text(xlims(1) + xlength/3,ylims(2) - ylength/10,label,'Interpreter','latex','FontSize',14)
        text(t(end)+10,Omg_H+ylength/40,sprintf('$\\Omega_H = %.4f$',Omg_H),'FontSize',10,'Interpreter','latex')

        f2 = figure;
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,Omg_th.^2,'LineWidth',my_linewidth/2,'LineStyle','-','DisplayName','$\Omega_\theta^2$','Color',MyColors('b1'))
        %plot(t,sth.*Omg_th,'LineWidth',my_linewidth,'DisplayName','$\Omega_\theta\sin\theta$','Color',MyColors('b2'))
        plot(t,(sth.*Omg_ph).^2,'LineWidth',my_linewidth/2,'LineStyle','-','DisplayName','$\sin^2\theta\Omega_\varphi^2$','Color',MyColors('r1'))
        plot(t,Omg.^2,'LineWidth',my_linewidth,'DisplayName','$\Omega^2$','Color','k')
        plot(t,Omg_orb.^2,'LineWidth',my_linewidth,'DisplayName','$(\Omega^{\rm orb})^2$','Color',MyColors('g1'))
        xlabel('$t$','FontSize',labels_fontsize,'Interpreter','latex')
        legend('Location','northwest','FontSize',legend_fontsize,'Interpreter','latex')
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(max_orb,'Color',[0 0 0],'LineStyle','-.','LineWidth',.5,'HandleVisibility','off')
        %yline(Omg_H,'Color',[0 0 0 .9],'LineStyle','--','HandleVisibility','off')
        yline(0,'Color',[.9 .9 .9],'LineWidth',.25,'HandleVisibility','off')
        xlim([t(end)-150,t(end)+30])

        xlims = xlim;
        ylims = ylim;
        xlength = xlims(2) - xlims(1);
        ylength = ylims(2) - ylims(1);
        text(xlims(1) + xlength/3,ylims(2) - ylength/10,label,'Interpreter','latex','FontSize',14)
        text(t(end)+10,Omg_H+ylength/40,sprintf('$\\Omega_H = %.4f$',Omg_H),'FontSize',10,'Interpreter','latex')

        f3 = figure;
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,Omg,'LineWidth',my_linewidth/2,'LineStyle','-','DisplayName','$\Omega$','Color','k')
        plot(t,Omg_orb,'LineWidth',my_linewidth,'DisplayName','$\Omega^{\rm orb}$','Color',MyColors('g1'))
        plot(t,Omg_so,'LineWidth',my_linewidth,'DisplayName','$\Omega^{\rm so}$','Color',MyColors('a1'))

        xlabel('$t$','FontSize',labels_fontsize,'Interpreter','latex')
        legend('Location','northwest','FontSize',legend_fontsize,'Interpreter','latex')
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(max_orb,'Color',[0 0 0],'LineStyle','-.','LineWidth',.5,'HandleVisibility','off')
        %yline(0,'Color',[.9 .9 .9],'LineWidth',.25,'HandleVisibility','off')
        xlim([t(end)-150,t(end)+20])

        xlims = xlim;
        ylims = ylim;
        xlength = xlims(2) - xlims(1);
        ylength = ylims(2) - ylims(1);
        text(xlims(1) + xlength/3,ylims(2) - ylength/10,label,'Interpreter','latex','FontSize',14)
    end

return