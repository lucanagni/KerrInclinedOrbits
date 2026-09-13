function SchwarzchildOrgy
    l = 2;
    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat')
    wf_eq = s;
    teq = wf_eq.ell(2).emm(3).t;
    heq_22 = wf_eq.ell(2).emm(3).hlm./norm;
    heq_21 = wf_eq.ell(2).emm(2).hlm./norm;

    load('~/waveforms/S/plunge/th0_60/wf.mat')
    wf_30 = DB_mode_rotate_timedep(s,s.dyn,2,2);
    t30 = wf_30.ell(2).emm(3).t;
    h30_22 = wf_30.ell(2).emm(3).hlm./norm;
    wf_30 = DB_mode_rotate_timedep(s,s.dyn,2,1);
    h30_21 = wf_30.ell(2).emm(2).hlm./norm;
    load('~/waveforms/S/plunge/th0_45/wf.mat')
    wf_45 = DB_mode_rotate_timedep(s,s.dyn,2,2);
    t45 = wf_45.ell(2).emm(3).t;
    h45_22 = wf_45.ell(2).emm(3).hlm./norm;
    wf_45 = DB_mode_rotate_timedep(s,s.dyn,2,1);
    h45_21 = wf_45.ell(2).emm(2).hlm./norm;

    DA_30 = abs(abs(h30_22) - abs(spline(teq,heq_22,t30)))./abs(spline(teq,heq_22,t30));
    DA_45 = abs(abs(h45_22) - abs(spline(teq,heq_22,t45)))./abs(spline(teq,heq_22,t45));

    % compute phase difference phi_rot - phi_eq
    phieq = -unwrap(angle(heq_22));
    phi30 = -unwrap(angle(h30_22));
    phi45 = -unwrap(angle(h45_22));

    sphi30 = spline(teq,phieq,t30);
    sphi45 = spline(teq,phieq,t45);

    DeltaPhi30 = sphi30 - phi30;
    DeltaPhi45 = sphi45 - phi45;

    DBeq = wf_eq.dyn;
    %DB30 = wf_30.dyn;
    %B45 = wf_45.dyn;

    my_fontsize = 11;
    my_linewidth = 1;

    tLR = tLR_splined(s.dyn);
    t_start = tLR-150;
    t_end = tLR+100;

    f1 = figure;
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex',FontSize=my_fontsize+2)
    xlabel('$u$','Interpreter','latex','FontSize',my_fontsize+2)
    %ylabel('${\Psi}_{\ell m}$','Interpreter','latex')
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',my_fontsize,'FontName','Times');

    f2 = figure;
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex',FontSize=my_fontsize+2)
    xlabel('$u$','Interpreter','latex','FontSize',my_fontsize+2)
    %ylabel('$M\omega_{\ell m}$','Interpreter','latex')
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',my_fontsize,'FontName','Times');

    f3 = figure;

    tl = tiledlayout(2, 1,'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl,1);
    set(ax1,'XMinorTick','on','YMinorTick','on','box','on','FontSize',my_fontsize,'FontName','Times');
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex',FontSize=my_fontsize+2)
    xlabel('$u$','Interpreter','latex','FontSize',my_fontsize+2)
    xlim([t_start,t_end])
    ax2 = nexttile(tl,2);
    set(ax2,'XMinorTick','on','YMinorTick','on','box','on','FontSize',my_fontsize,'FontName','Times');
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex',FontSize=my_fontsize+2)
    xlabel('$u$','Interpreter','latex','FontSize',my_fontsize+2)
    %ylabel('$M\omega_{\ell m}$','Interpreter','latex')
    xlim([t_start,t_end])


    %teq = teq + Shift(0,7,0);
    %t30 = t30 + Shift(0,7,0);
    %t45 = t45 + Shift(0,6.5,0);

    dt22 = find_deltat(heq_22,teq,h45_22,t45);      
    dt21 = find_deltat(heq_21,teq,h45_21,t45);      
    dt21 = 0;
    dt22 = 0;




    r1 = [255,247,188]/255;
    r2 = [254,196,79]/255;
    r3 = [217,95,14]/255;

    b1 = [236,231,242]/255;
    b2 = [166,189,219]/255;
    b3 = [43,140,190]/255;

    figure(f1)
    hold on
    plot(teq,abs(heq_22),'DisplayName','$|{\Psi}_{22}|/\nu, \ \iota = 0$','LineWidth',my_linewidth,'Color',r1);
    plot(t30,abs(h30_22),'DisplayName','$|\tilde{\Psi}_{22}|/\nu, \ \iota = \pi/5$','LineWidth',my_linewidth,'LineStyle','--','Color',r2);
    plot(t45 - dt22,abs(h45_22),'DisplayName','$|\tilde{\Psi}_{22}|/\nu, \ \iota = \pi/4$','LineWidth',my_linewidth,'LineStyle','-.','Color',r3);

    plot(teq,abs(heq_21),'DisplayName','$|{\Psi}_{21}|/\nu, \ \iota = 0$','LineWidth',my_linewidth,'Color',b1);
    plot(t30,abs(h30_21),'DisplayName','$|\tilde{\Psi}_{21}|/\nu, \ \iota = \pi/6$','LineWidth',my_linewidth,'LineStyle','--','Color',b2);
    plot(t45 - dt21,abs(h45_21),'DisplayName','$|\tilde{\Psi}_{21}|/\nu, \ \iota = \pi/4$','LineWidth',my_linewidth,'LineStyle','-.','Color',b3);


    plot(DBeq.t,2*DBeq.Omg,'LineWidth',1,'Color','k','DisplayName','$2 \Omega$','Color',[.8 .8 .8])
    xline(tLR,'LineStyle','--','Color',[.5 .5 .5], 'HandleVisibility', 'off')
    xlim([t_start,t_end])
    ylim auto
    hold off

    figure(f2)
    hold on
    plot(teq,freq(heq_22,teq),'DisplayName','$\omega_{22}, \ \iota = 0$','LineWidth',my_linewidth,'Color',r1);
    plot(t30,freq(h30_22,t30),'DisplayName','$\tilde\omega_{22}, \ \iota = \pi/6$','LineWidth',my_linewidth,'LineStyle','--','Color',r3);
    plot(t45 - dt22,freq(h45_22,t45),'DisplayName','$\tilde\omega_{22}, \ \iota = \pi/4$','LineWidth',my_linewidth,'LineStyle','-.','Color',r3);

    plot(teq,freq(heq_21,teq),'DisplayName','$\omega_{21}, \ \iota = 0$','LineWidth',my_linewidth,'Color',b1);
    plot(t30,freq(h30_21,t30),'DisplayName','$\tilde\omega_{21}, \ \iota = \pi/6$','LineWidth',my_linewidth,'LineStyle','--','Color',b2);
    plot(t45 - dt21,freq(h45_21,t45),'DisplayName','$\tilde\omega_{21}, \ \iota = \pi/4$','LineWidth',my_linewidth,'LineStyle','-.','Color',b3);

    plot(DBeq.t,2*DBeq.Omg,'LineWidth',1,'Color','k','DisplayName','$2 \Omega$','Color',[.8 .8 .8])
    xline(tLR,'LineStyle','--','Color',[.5 .5 .5], 'HandleVisibility', 'off')
    xlim([t_start,t_end])
    ylim auto
    hold off

    figure(f3)
    hold(ax1,"on")
    plot(ax1,t30,DA_30,'DisplayName','$\iota = \pi/6$','LineWidth',my_linewidth,'Color',[0 0 0 0.9]);
    plot(ax1,t45,DA_45,'DisplayName','$\iota = \pi/4$','LineWidth',my_linewidth,'Color',[0 0 0 0.6]);

    return