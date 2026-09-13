function MultipolarSchw
    l = 2;
    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2));

    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+4;
    label_fontsize = 18;
    load('~/waveforms/S/plunge/th0_60/wf.mat');
    h22 = s.ell(2).emm(3).hlm./norm;
    h21 = s.ell(2).emm(2).hlm./norm;
    h20 = s.ell(2).emm(1).hlm./norm;
    h2m1 = s.ell(2).emminus(2).hlm./norm;
    h2m2 = s.ell(2).emminus(3).hlm./norm;
    t = s.ell(2).emm(3).t;

    [~,t_peak] = Apeak(h22,t);
    t_in = t_peak-250;
    t_end = t_peak + 150;

    tLR = tLR_splined(s.dyn);

    dyn = s.dyn;

    f1 = figure;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u$','Interpreter','latex','FontSize',label_fontsize)
    ylabel('$|\Psi_{\ell m}|/\nu$','Interpreter','latex','FontSize',label_fontsize)
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    xlim([t_in,t_end])
    hold on
    plot(t,abs(h22),'DisplayName','(2,2)','LineWidth',1,'Color',MyColors('r1'))
    plot(t,abs(h21),'DisplayName','(2,1)','LineWidth',1,'Color',MyColors('b1'))
    plot(t,abs(h20),'DisplayName','(2,0)','LineWidth',1,'Color',[.5 .5 .5])
    plot(t,abs(h2m1),'DisplayName','(2,-1)','LineWidth',1.5,'LineStyle',':','Color',MyColors('b1'))
    plot(t,abs(h2m2),'DisplayName','(2,-2)','LineWidth',1.5,'LineStyle',':','Color',MyColors('r1'))
    %plot(dyn.t,dyn.th/10,'Color',[.7 .7 .7],'DisplayName','$\theta (u)$')
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    %plot(dyn.t,mod(dyn.phi,2*pi)/20,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')

    
    f2 = figure;
    f2.Position(1) = f1.Position(1) + f1.Position(3) + 50;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u-u_{LR}$','Interpreter','latex','FontSize',label_fontsize)
    ylabel('$\omega_{\ell m}$','Interpreter','latex','FontSize',label_fontsize)
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','southwest')
    xlim([t_in-tLR,t_end-tLR])
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    hold on
    plot(t-tLR,freq(h22,t),'DisplayName','$(2,2)$','LineWidth',1,'Color',MyColors('r1'))
    plot(t-tLR,freq(h21,t),'DisplayName','$(2,1)$','LineWidth',1,'Color',MyColors('b1'))
    plot(t-tLR,freq(h20,t),'DisplayName','$(2,0)$','LineWidth',1.5,'Color',[.5 .5 .5])
    plot(t-tLR,freq(h2m1,t),'DisplayName','$(2,-1)$','LineWidth',1.5,'LineStyle',':','Color',MyColors('b1'))
    plot(t-tLR,freq(h2m2,t),'DisplayName','$(2,-2)$','LineWidth',1.5,'LineStyle',':','Color',MyColors('r1'))
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    text(t_end-40,2.5,'$\iota=\pi/3$','Interpreter','latex','FontSize',16)

    return


    