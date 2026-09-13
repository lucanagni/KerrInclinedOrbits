function S_errors
    label_fontsize = 18;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize;
    my_linewidth = 1;

    % Normalization factor for RWZ function
    l = 2;
    norm = sqrt(factorial(l+2)./factorial(l-2));

    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat')
    wf_eq = s;
    teq = wf_eq.ell(2).emm(3).t;

    load('~/waveforms/S/plunge/th0_60/wf.mat')
    wf_30 = s;
    t30 = wf_30.ell(2).emm(3).t;

    load('~/waveforms/S/plunge/th0_45/wf.mat')
    wf_45 = s;
    t45 = wf_45.ell(2).emm(3).t;

    load('~/waveforms/S/plunge/th0_30/wf.mat')
    wf_60 = s;
    t60 = wf_60.ell(2).emm(3).t;

    load('~/waveforms/S/plunge/th0_30/hires/wf.mat')
    wf_red = s;
    t_red = wf_red.ell(2).emm(3).t;

    t = t30;

    tLR = tLR_splined(wf_eq.dyn);
    t_start = tLR-250;
    t_end = tLR+100;

    wf_30r = struct;
    wf_45r = struct;
    wf_60r = struct;
    wf_redr = struct;

    for m=1:2
        heq = wf_eq.ell(2).emm(m+1).hlm./norm;
        heq = spline(teq,heq,t);

        wf_30r = DB_mode_rotate_timedep(wf_30,wf_30.dyn,2,m,wf_30r);
        h30 = wf_30r.ell(2).emm(m+1).hlm./norm;
        h30 = spline(t30,h30,t);

        wf_45r = DB_mode_rotate_timedep(wf_45,wf_45.dyn,2,m,wf_45r);
        h45 = wf_45r.ell(2).emm(m+1).hlm./norm;
        h45 = spline(t45,h45,t);

        wf_60r = DB_mode_rotate_timedep(wf_60,wf_60.dyn,2,m,wf_60r);
        h60 = wf_60r.ell(2).emm(m+1).hlm./norm;
        h60 = spline(t60,h60,t);

        wf_redr = DB_mode_rotate_timedep(wf_red,wf_red.dyn,2,m,wf_redr);
        h60r = wf_redr.ell(2).emm(m+1).hlm./norm;
        h60r = spline(t_red,h60r,t);

        DA_30 = abs(abs(h30) - abs(heq))./abs(heq);
        DA_45 = abs(abs(h45) - abs(heq))./abs(heq);
        DA_60 = abs(abs(h60) - abs(heq))./abs(heq);
        DA_60r = abs(abs(h60r) - abs(heq))./abs(heq);


        % compute phase difference phi_rot - phi_eq
        phieq = -unwrap(angle(heq));
        phi30 = -unwrap(angle(h30));
        phi45 = -unwrap(angle(h45));
        phi60 = -unwrap(angle(h60));
        phi60r = -unwrap(angle(h60r));

        DeltaPhi30 = phi30 - phieq;
        DeltaPhi45 = phi45 - phieq;
        DeltaPhi60 = phi60 - phieq;
        DeltaPhi60r = phi60r - phieq + 2*pi;
        figure

        tl(m) = tiledlayout(2, 1,'TileSpacing', 'compact', 'Padding', 'compact');
        ax(1,m) = nexttile(tl(m),1);
        set(ax(1,m),'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
        legend(BackgroundAlpha=.7,Location='northeast',Interpreter='latex',FontSize=legend_fontsize)
        ylabel(sprintf('$|\\Delta A_{%d%d}|/A_{%d%d}$',l,m,l,m),'Interpreter','latex','FontSize',label_fontsize)
        xlim([t_start,t_end])
        ax(2,m) = nexttile(tl(m),2);
        set(ax(2,m),'XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
        %legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex',FontSize=legend_fontsize)
        xlabel('$u$','Interpreter','latex','FontSize',label_fontsize)
        ylabel(sprintf('$\\Delta \\phi_{%d%d}$',l,m),'Interpreter','latex','FontSize',label_fontsize)
        xlim([t_start,t_end])    
        
        ax1 = ax(1,m);
        hold(ax1,'on');
        plot(ax1,t,DA_30,'DisplayName','$\iota = \pi/6$','LineWidth',my_linewidth,'Color',[0 0 0 0.9],'LineStyle','-');
        plot(ax1,t,DA_45,'DisplayName','$\iota = \pi/4$','LineWidth',my_linewidth,'Color',[0 0 0 0.6],'LineStyle','--');
        plot(ax1,t,DA_60,'DisplayName','$\iota = \pi/3$','LineWidth',my_linewidth,'Color',[0 0 0 0.3],'LineStyle','-.');
        plot(ax1,t,DA_60r,'DisplayName','$\iota = \pi/3^*$','LineWidth',my_linewidth/2,'Color','r','LineStyle','-');
        


        ax2 = ax(2,m);
        hold(ax2,'on');
        plot(ax2,t,DeltaPhi30,'DisplayName','$\iota = \pi/6$','LineWidth',my_linewidth,'Color',[0 0 0],'LineStyle','-');
        plot(ax2,t,DeltaPhi45,'DisplayName','$\iota = \pi/4$','LineWidth',my_linewidth,'Color',[0 0 0 0.7],'LineStyle','--');
        plot(ax2,t,DeltaPhi60,'DisplayName','$\iota = \pi/3$','LineWidth',my_linewidth,'Color',[0 0 0 0.3],'LineStyle','-.');
        plot(ax2,t,DeltaPhi60r,'DisplayName','$\iota = \pi/3^*$','LineWidth',my_linewidth/2,'Color','r','LineStyle','-');
        grid on

        if m==1
            % Force scientific notation for y-axis
            ax1.YAxis.Exponent = -2;  % Set exponent explicitly
            ax1.YAxis.TickLabelFormat = '%.0f';  % Show only one significant digit
            ax2.YAxis.Exponent = -2;  % Set exponent explicitly
            ax2.YAxis.TickLabelFormat = '%.0f';  % Show only one significant digit
        end 

        %ax1.YAxis.Exponent = 0; 
        %ax2.YAxis.Exponent = 0;
        %tickformat('%.1e');  
        %ax1.TickLabelInterpreter = 'none'; 
        %ax2.TickLabelInterpreter = 'none'; 

    end



    %figure
    %abs(h60)
    %plot(t,abs(h60))
    %hold on
    %plot(t,abs(heq))

return