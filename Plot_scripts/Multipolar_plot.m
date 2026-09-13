function Multipolar_plot(s,p,l,m)

    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    h1 = s.ell(l).emm(m+1).hlm./norm;
    t1 = s.ell(l).emm(m+1).t;
    h2 = p.ell(l).emm(m+1).hlm./norm;
    t2 = p.ell(l).emm(m+1).t;

    dt = 0*find_deltat(h1,t1,h2,t2); %compute dt if needed
    
    %amplitudes and frequencies of equatorial and tilted wf
    A1 = abs(h1); 
    A2 = abs(h2);
    Omg1 = freq(h1,t1);
    Omg2 = freq(h2,t2);

    tLR = tLR_splined(s.dyn);

    %compute rotated waveform with ampl and freq
    r = DB_mode_rotate_timedep(p,p.dyn,l,m);
    h2_r = r.ell(l).emm(m+1).hlm./norm;
    A2_r = abs(h2_r);
    Omg2_r = freq(h2_r,t2);
    
    % compute phase difference phi_rot - phi_eq
    phi1 = -unwrap(angle(h1));
    phi2 = -unwrap(angle(h2_r));

    sphi2 = spline(t2-dt,phi2,t1);
    DeltaPhi = sphi2 - phi1;

    %idx0 = find(t1>50,1);
    %idx1 = find(t1>100,1);
    %[~,idx_align] = max(real(h1(idx0:idx1)));
    %idx_align = idx0 + idx_align;

    %delta_phi = sphi2(idx_align) - phi1(idx_align);
    %h1 = h1.*exp(-1i*delta_phi);
    %DeltaPhi = sphi2 - phi1 - 0*delta_phi;

    % compute relative amplitude difference between rotated and eq wf
    rel_A_diff = abs(A1-spline(t2,A2_r,t1))./A1;

    %Make plots look beautiful
    axes_fs = 10;
    labels_fs = 14;
    legend_fs = 14;
    my_linewidth = 1;

    t_start = 80;
    t_end = tLR+110;

    f1 = figure; %First plot contains real part of the two waveforms
    f1.Position = [100 558 560 420];
    T = tiledlayout(2,2,'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(T,1,[1,2])
    ax1 = gca;
    set(ax1,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    ylabel('$\Re[\Psi_{22}]/\nu$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,real(h1),'Color','r','LineWidth',my_linewidth)
    text(350,0.35,'$\iota = 0$','Interpreter','latex','FontSize',16)
    xlim([0,t_end])
    lims = ylim;
    ylim([lims(1),-lims(1)])

    nexttile(T,3,[1,2])
    ax12 = gca;
    set(ax12,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    text(350,0.3,'$\iota = \pi/6$','Interpreter','latex','FontSize',16)
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)
    ylabel('$\Re[\Psi_{22}]/\nu$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t2,real(h2),'Color','b','LineStyle','--','LineWidth',my_linewidth)
    xlim([0,t_end])
    lims = ylim;
    ylim([lims(1),-lims(1)])

    %nexttile(T,3,[2,2])
    f2 = figure; %Second plot contains amplitudes 
    f2.Position = [100+560+10 558 560 420];
    ax2 = gca;
    set(ax2,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    legend('BackgroundAlpha',.7,'Location','northeast','Interpreter','latex','FontSize',legend_fs)
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,A1,'DisplayName','$|\Psi_{22}|/\nu \ (\iota = 0)$','LineWidth',my_linewidth,'Color','r')
    plot(t2 - dt,A2,'DisplayName','$|\Psi_{22}|/\nu \ (\iota = \pi/6)$','LineWidth',my_linewidth,'Color','b','LineStyle','--')
    plot(s.dyn.t,2*s.dyn.Omg,'DisplayName','$2\Omega$','LineStyle','-','LineWidth',my_linewidth,'Color',[.8 .8 .8])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle',':','HandleVisibility','off')
    xlim([tLR-100,t_end])

    f3 = figure; %Third plot contains frequencies
    f3.Position = [100+2*(560+10) 558 560 420];
    ax3 = gca;
    set(ax3,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    legend('BackgroundAlpha',.7,'Location','southeast','Interpreter','latex','FontSize',legend_fs)
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,Omg1,'DisplayName','$\omega_{22} \ (\iota = 0)$','Color','r','LineWidth',my_linewidth)
    plot(t2 - dt,Omg2,'DisplayName','$\omega_{22} \ (\iota = \pi/6)$','Color','b','LineWidth',my_linewidth,'LineStyle','--')
    plot(s.dyn.t,2*s.dyn.Omg,'DisplayName','$2\Omega$','LineStyle','-','LineWidth',my_linewidth,'Color',[.8 .8 .8])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle',':','HandleVisibility','off')
    xlim([tLR-100,t_end])

    f4 = figure; %Fourth plot contains rotated waveforms
    f4.Position = [100 558-420-20 560 420];
    ax4 = gca;
    set(ax4,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    legend('BackgroundAlpha',.7,'Location','southwest','Interpreter','latex','FontSize',legend_fs-1)
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,A1,'DisplayName','$|\Psi_{22}|/\nu \ (\iota = 0)$','LineWidth',my_linewidth,'Color','r')
    plot(t2 - dt,A2_r,'DisplayName','$|\tilde{\Psi}_{22}|/\nu \ (\iota = \pi/6)$','LineWidth',my_linewidth,'Color','b','LineStyle','--')
    plot(t1,Omg1,'DisplayName','$\omega_{22} \ (\iota = 0)$','Color','r','LineWidth',my_linewidth)
    plot(t2 - dt,Omg2_r,'DisplayName','$\tilde{\omega}_{22} \ (\iota = \pi/6)$','Color','b','LineWidth',my_linewidth,'LineStyle','--')
    plot(s.dyn.t,2*s.dyn.Omg,'DisplayName','$2\Omega$','LineStyle','-','LineWidth',my_linewidth,'Color',[.8 .8 .8])
    xline(tLR,'Color',[.7 .7 .7],'LineStyle',':','HandleVisibility','off')
    xlim([tLR-100,t_end])

    f5 = figure;
    f5.Position = [100+560+10 558-420-20 560 420];
    ax5 = gca;
    set(ax5,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    legend('BackgroundAlpha',.7,'Location','southeast','Interpreter','latex','FontSize',legend_fs)
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)
    ylabel('$\Re[\Psi_{22}]/\nu$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,real(h1),'Color','r','LineWidth',my_linewidth,'DisplayName','$\Psi_{22} \ (\iota = 0)$')
    plot(t2,real(h2_r),'Color','b','LineStyle','--','LineWidth',my_linewidth,'DisplayName','$\tilde{\Psi}_{22} \ (\iota = \pi/6)$')
    xlim([0,t_end])

    f6 = figure; %sixth plot shows phase diff and ampl rel diff
    f6.Position = [100+2*(560+10) 558-420-20 560 420];
    T2 = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(T2,1)
    ax6 = gca;
    set(ax6,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');

    hold on
    plot(t1,rel_A_diff,'LineWidth',my_linewidth,'Color','k')
    ylabel('$\Delta\mathcal{A}/A^{EQ}$','Interpreter','latex','FontSize',labels_fs)
    xlim([0,t_end])

    nexttile(T2,2)
    ax62 = gca;
    set(ax62,'XMinorTick','on','YMinorTick','on','FontSize',axes_fs,'FontName','Times','box','on');
    xlabel('$u/M$','Interpreter','latex','FontSize',labels_fs)

    hold on
    plot(t1,DeltaPhi,'LineWidth',my_linewidth,'Color','k')
    ylabel('$\Delta\phi_{22}$','Interpreter','latex','FontSize',labels_fs)
    xlim([0,t_end])

return