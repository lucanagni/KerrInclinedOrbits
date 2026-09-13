function S_hplus_plot
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    % Normalization factor for RWZ function
    %norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    load('~/waveforms/S/plunge/th0_30/hires/wf.mat')
    wf_red = s;
    
    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat');
    wf_eq = s;
    load('~/waveforms/S/plunge/th0_60/wf.mat');
    wf_60 = s;
    load('~/waveforms/S/plunge/th0_45/wf.mat');
    wf_45 = s;
    load('~/waveforms/S/plunge/th0_30/wf.mat');
    wf_30 = s;

    
    %if (s.dyn.th0-pi/6)<1e-5
    %    theta_legend = '\pi/12';
    %    phi_legend = '\pi';
    %    theta = pi/12;
    %    phi = pi;
    %    iota = '\pi/3';
    %elseif (s.dyn.th0-pi/4)<1e-5
    %    theta_legend = '0';
    %    phi_legend = '0';
    %    theta = 0;
    %    phi = 0;
    %    iota = '\pi/4';
    %elseif (s.dyn.th0-pi/3)<1e-5
        theta_legend = '\pi/12';
        phi_legend = '0';
        theta = pi/12;
        phi = 0;
        iota = '\pi/6';
    %end

    t_peak = tLR_splined(wf_60.dyn);
    t_start = t_peak-250;
    t_end = t_peak+100;

    [heq,teq] = DB_hpc(wf_eq,pi/4,0,'norm','none','ellmax',4);
    [h60,t] = DB_hpc(wf_60,pi/12,0,'norm','none','ellmax',4);
    [h45,t45] = DB_hpc(wf_45,0,0,'norm','none','ellmax',4);
    [h30,t30] = DB_hpc(wf_30,pi/12,pi,'norm','none','ellmax',4);
    [hred,tred] = DB_hpc(wf_red,pi/12,pi,'norm','none','ellmax',4);
    heq = spline(teq,heq,t);
    h45 = spline(t45,h45,t);
    h30 = spline(t30,h30,t);
    hred = spline(tred,hred,t);

    DeltaPhi60 = -unwrap(angle(h60)) + unwrap(angle(heq));
    DeltaPhi45 = -unwrap(angle(h45)) + unwrap(angle(heq));
    DeltaPhi30 = -unwrap(angle(h30)) + unwrap(angle(heq));
    DeltaPhired = -unwrap(angle(hred)) + unwrap(angle(heq)) + 2*pi;


    DeltaA60 = abs(abs(heq)-abs(h60))./abs(heq);
    DeltaA45 = abs(abs(heq)-abs(h45))./abs(heq);
    DeltaA30 = abs(abs(heq)-abs(h30))./abs(heq);
    DeltaAred = abs(abs(heq)-abs(hred))./abs(heq);
    

    f1 = figure;
    %f1.Position(4) = .75*f1.Position(4);
    %tl = tiledlayout(2,2,'TileSpacing', 'compact', 'Padding', 'compact');

    %ax = nexttile(tl,1,[2,1]);  
    ax = gca;
    set(ax,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on

    plot(ax,t,real(heq),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\iota = 0,\\quad \\ \\,(\\Theta = \\pi/4,\\Phi = 0)$'),'Color',[0 0 0 0.5])
    plot(ax,t,real(h60),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\iota = %s,\\ (\\Theta = %s,\\Phi=%s)$',iota,theta_legend,phi_legend),'LineStyle','--','Color',[0 0 0])
    %plot(ax,t,real(hdebug),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\iota = 0,\\quad \\ \\,(\\Theta = \\pi/4,\\Phi = 0)$ - EOB'),'LineStyle','-.','Color','r')
    xline(t_peak,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel('$D_L h_+/\mu$','FontSize',labels_fontsize,'Interpreter','latex')
    %ylim([-1,1])

    xlim([t_start,t_end])
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')   

    f2 = figure;
    tl = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl,1)
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,DeltaA60,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/6$','Color',[0 0 0 0.9],'LineStyle','-')
    plot(t,DeltaA45,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/4$','Color',[0 0 0 0.6],'LineStyle','--')
    plot(t,DeltaA30,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/3$','Color',[0 0 0 0.3],'LineStyle','-.')
    plot(t,DeltaAred,'LineWidth',my_linewidth/2,'DisplayName','$\iota = \pi/3^*$','Color','r','LineStyle','-');
    legend('Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'Location','northwest')   
    xlim([t_start,t_end])
    ylabel(sprintf('$| \\Delta A_h |$'), 'Interpreter', 'latex','FontSize',labels_fontsize)
    %xlabel('$u$','Interpreter','latex')

    %ymax = ylim;
    %ylim([0,1.2*ymax(2)])

    grid on

    nexttile(tl,2)
    set(gca,'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
    hold on
    plot(t,DeltaPhi60,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/6$','Color',[0 0 0 0.9],'LineStyle','-')
    plot(t,DeltaPhi45,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/4$','Color',[0 0 0 0.6],'LineStyle','--')
    plot(t,DeltaPhi30,'LineWidth',my_linewidth,'DisplayName','$\iota = \pi/3$','Color',[0 0 0 0.3],'LineStyle','-.')
    plot(t,DeltaPhired,'LineWidth',my_linewidth/2,'DisplayName','$\iota = \pi/3^*$','Color','r','LineStyle','-');
    xlim([t_start,t_end])
    ylabel(sprintf('$ \\Delta \\phi_h$'), 'Interpreter', 'latex','FontSize',labels_fontsize)
    xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize)
    %legend('Interpreter','latex','FontSize',legend_fontsize+4,'BackgroundAlpha',0.7,'Location','northwest')   
    %ymax = ylim;
    %ylim([1.2*ymax(1),1.2*ymax(2)])

    grid on



    %A4Width(f1)

return