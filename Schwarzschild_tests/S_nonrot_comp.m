function S_nonrot_comp(l,m)
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize;
    labels_fontsize = 14;

    % Normalization factor for RWZ function
    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat');
    wf_eq = s;
    load('~/waveforms/S/plunge/th0_60/wf.mat');
    wf_prec = s;

    addpath '/home/luca/repos/kerrorbitsolver/KerrDynamics' 
    omega0 = sign(0.5*(2*m-1))*imag(kerr_FitKerrQNMs(l,2,0,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));

    tLR = tLR_splined(wf_prec.dyn);

    heq = wf_eq.ell(l).emm(m+1).hlm./norm;
    teq = wf_eq.ell(l).emm(m+1).t;
    
    hp = wf_prec.ell(l).emm(m+1).hlm./norm;
    tp = wf_prec.ell(l).emm(m+1).t;

    heq = spline(teq,heq,tp);
    t = tp;

    [peak,t_peak] = Apeak(wf_prec.ell(2).emm(3).hlm,t);

    t_start = t_peak-200;
    t_end = t_peak+140;

    r1 = MyColors('r1');
    r2 = MyColors('r2');
    b1 = MyColors('b1');
    b2 = MyColors('b2');

    f = figure;
    f.Position(3:4) = [ 560   420];
    tl = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');
    axl = nexttile(tl,1);  
    set(axl,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(axl,t,abs(heq),'LineWidth',my_linewidth+0.3,'DisplayName','$\iota = 0$','LineStyle','-','Color',r2)
    plot(axl,t,abs(hp),'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/6$','LineStyle','--','Color',r1)
    xline(tLR,'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')
    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    %if m==1
    %    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','southeast')
    %end
    ylabel(sprintf('$|\\Psi_{%d%d}|/\\nu$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)

    axr = nexttile(tl,2);  
    set(axr,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(axr,t,freq(heq,t),'LineWidth',my_linewidth+0.3,'DisplayName','$\iota = 0$','LineStyle','-','Color',r2)
    plot(axr,t,freq(hp,t),'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/6$','LineStyle','--','Color',r1)
    xline(tLR,'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')
    xlim([t_start,t_end])
    yline(omega0,'LineStyle','--','HandleVisibility','off')
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')
    ylabel(sprintf('$\\omega_{%d%d}$',l,m), 'Interpreter', 'latex','FontSize',labels_fontsize)
    xlabel('$u$','Interpreter', 'latex','FontSize',labels_fontsize)
    xlength = t_end-t_start;
    lims = ylim;
    ylength = abs(lims(2)-lims(1));

    text(xlength/3+t_start,omega0 + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omega0),'Interpreter','latex','FontSize',14)

    %A4Width(f);
return