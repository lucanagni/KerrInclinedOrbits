function FitHughes
    load('~/waveforms/K/plunge/th0_30/a05/wf.mat');

    l = 2;
    m = 2;
    j = 0.5;
    A220p = 1.9629e+00;
    A220m = 2.2288e+00;
    A320p = 7.6702e-01;
    phi220p = 2.0825e-01;
    phi220m = 1.2353e+00;
    phi320p = 8.7733e-01;
    
    C220p = A220p.*exp(1i.*phi220p);
    C220m = A220m.*exp(1i.*phi220m);
    C320p = A320p.*exp(1i.*phi320p);

    sigma220p = kerr_FitKerrQNMs(l,abs(m),j,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma220m = kerr_FitKerrQNMs(l,-abs(m),j,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma320p = kerr_FitKerrQNMs(l+1,abs(m),j,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');

    h = s.ell(l).emm(m).hlm;
    t = s.ell(l).emm(m).t;

    tLR = tLR_splined(s.dyn);

    tin = tLR+10;
    tend = tLR+120;

    idx0 = find(t>tin,1);
    idx1 = find(t>tend,1);

    h_ringdown = h(idx0:idx1);
    t_ringdown = t(idx0:idx1);
    %t = t_ringdown;

    h_fit = C220p.*exp(-sigma220p.*t_ringdown) + C220m.*exp(-sigma220m.*t_ringdown) + C320p.*exp(-sigma320p.*t_ringdown);

    phi1 = -unwrap(angle(h_ringdown));
    phi2 = -unwrap(angle(h_fit));
    Dphi = (phi2(end)-phi1(end));
    ratio = h_ringdown./h_fit;

    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(t_ringdown,abs(exp(-1i.*Dphi).*h_ringdown),'LineWidth',1,'DisplayName','Numerical')
    plot(t_ringdown,abs(h_fit.*ratio(1)),'LineWidth',1,'DisplayName','Hughes Data')
    ylabel('$\Re[h_{22}]$','FontSize',14,'Interpreter','latex')
    xlabel('$u$','FontSize',14,'Interpreter','latex')
    legend('Interpreter','latex','FontSize',14,'BackgroundAlpha',0.7,'Location','northeast')   
    xlim([tin,tend])

    figure
    hold on
    plot(t_ringdown,-unwrap(angle(h_ringdown)))
    plot(t_ringdown,-unwrap(angle(h_fit)))
    xlim([tin,tend])
