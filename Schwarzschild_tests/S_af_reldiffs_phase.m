function S_af_reldiffs_phase(l,m,incl)
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat');
    wf_eq = s;
    load(sprintf('~/waveforms/S/plunge/th0_%s/wf.mat',num2str(incl)));
    wf_prec = s;
    wf_rot = DB_mode_rotate_timedep(wf_prec,wf_prec.dyn,l,m);

    heq = wf_eq.ell(l).emm(m+1).hlm./norm;
    teq = wf_eq.ell(l).emm(m+1).t;
    
    hr = wf_rot.ell(l).emm(m+1).hlm./norm;
    tr = wf_rot.ell(l).emm(m+1).t;

    %dt = find_deltat(heq,teq,hr,tr);
    dt = 0;
    tr = tr - dt;

    heq = spline(teq,heq,tr);
    t = tr;

    % compute phase difference phi_rot - phi_eq
    phi1 = -unwrap(angle(heq));
    phi2 = -unwrap(angle(hr));

    DeltaPhi = phi2 - phi1;

    %[peak,t_peak] = Apeak(hr,t);
    t_peak = tLR_splined(s.dyn);

    t_start = t_peak-250;
    t_end = t_peak+100;
    %t_start = t_peak-130;
    %t_end = t_peak+70;

    f = figure;
    f.Position(4) = 1.3*f.Position(4);
    tl = tiledlayout(5,2,'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl,1,[3,2]);  
    set(ax,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(ax,t,abs(heq),'LineWidth',my_linewidth,'DisplayName',sprintf('$|\\Psi_{%d%d}|/\\mu$',l,m),'LineStyle','-','Color',[254,196,79]./255)
    plot(ax,t,abs(hr),'LineWidth',my_linewidth,'DisplayName',sprintf('$|\\tilde\\Psi_{%d%d}|/\\mu$',l,m),'LineStyle','--','Color',[217,95,14]./255)
    plot(ax,t,freq(heq,t),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\omega_{%d%d}$',l,m),'LineStyle','-','Color',[166,189,219]./255)
    plot(ax,t,freq(hr,t),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\tilde\\omega_{%d%d}$',l,m),'LineStyle','--','Color',[43,140,190]./255)
    plot(ax,t,real(heq),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[\\Psi_{%d%d}]/\\mu$',l,m),'Color',[0 0 0 0.1])
    plot(ax,t,real(hr),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[\\tilde\\Psi_{%d%d}]/\\mu$',l,m),'Color',[0 0 0 0.2],'LineStyle','--')
    %xline(tLSSO_splined(s.dyn),'Color',[.5 .5 .5],'LineStyle','--','HandleVisibility','off')
    %xline(tLR_splined(s.dyn),'Color',[.5 .5 .5],'LineStyle','--','HandleVisibility','off')

    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')   
    %ylabel(sprintf('$|\\Psi_{%d%d}|/\\mu$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)

    nexttile(tl,7,[1,2])
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,abs(abs(heq)-abs(hr))./abs(heq),'LineWidth',my_linewidth,'Color','k')
    xlim([t_start,t_end])
    ylabel(sprintf('$| \\Delta A_{%d%d} |$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)
    %xlabel('$u$','Interpreter','latex')

    ymax = ylim;
    ylim([0,1.2*ymax(2)])

    grid on

    nexttile(tl,9,[1,2])
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,DeltaPhi,'LineWidth',my_linewidth,'Color','k')
    xlim([t_start,t_end])
    ylabel(sprintf('$ \\Delta \\phi_{%d%d}$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)

    ymax = ylim;
    ylim([1.2*ymax(1),1.2*ymax(2)])

    ay = gca;
    ticks = ay.YTick;
    ay.YTick = ticks(1:2:end); 
    grid on

return
