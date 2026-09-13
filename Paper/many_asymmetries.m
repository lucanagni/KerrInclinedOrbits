function many_asymmetries
    % ==========================================================================================================
    % Produce Fig. 16 of Paper I
    % Assumes Schwarzschild waveforms are located in <waveformdir>/S/plunge/th0_<value>/wf.mat
    % Assumes Kerr waveforms are located in <waveformdir>/K/plunge/th0_<value>/a0<value>/wf.mat
    % wf.mat is a struct containing the dynamics (element of DB_class) as wf.dyn and modes as wf.ell(l).emm(m+1).hlm (emminus for negative m)
    % ==========================================================================================================

    waveformdir = '/home/luca/waveforms';

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = .8;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;


    a00i0 = load(sprintf('%s/S/plunge/th0_90/more_timesteps/wf.mat',waveformdir));
    a00i30 = load(sprintf('%s/S/plunge/th0_60/wf.mat',waveformdir));
    a00i45 = load(sprintf('%s/S/plunge/th0_45/wf.mat',waveformdir));
    a00i60 = load(sprintf('%s/S/plunge/th0_30/wf.mat',waveformdir));

    a02i0 = load(sprintf('%s/K/plunge/th0_90/a02/wf.mat',waveformdir));
    a02i30 = load(sprintf('%s/K/plunge/th0_60/a02/wf.mat',waveformdir));
    a02i45 = load(sprintf('%s/K/plunge/th0_45/a02/wf.mat',waveformdir));
    a02i60 = load(sprintf('%s/K/plunge/th0_30/a02/wf.mat',waveformdir));
    a02i120 = load(sprintf('%s/K/plunge/th0_30/a-02/wf.mat',waveformdir));
    a02i135 = load(sprintf('%s/K/plunge/th0_45/a-02/wf.mat',waveformdir));
    a02i150 = load(sprintf('%s/K/plunge/th0_60/a-02/wf.mat',waveformdir));
    a02i180 = load(sprintf('%s/K/plunge/th0_90/a-02/wf.mat',waveformdir));

    a05i0 = load(sprintf('%s/K/plunge/th0_90/a05/wf.mat',waveformdir));
    a05i30 = load(sprintf('%s/K/plunge/th0_60/a05/wf.mat',waveformdir));
    a05i45 = load(sprintf('%s/K/plunge/th0_45/a05/wf.mat',waveformdir));
    a05i60 = load(sprintf('%s/K/plunge/th0_30/a05/wf.mat',waveformdir));
    a05i120 = load(sprintf('%s/K/plunge/th0_30/a-05/wf.mat',waveformdir));
    a05i135 = load(sprintf('%s/K/plunge/th0_45/a-05/wf.mat',waveformdir));
    a05i150 = load(sprintf('%s/K/plunge/th0_60/a-05/wf.mat',waveformdir));
    a05i180 = load(sprintf('%s/K/plunge/th0_90/a-05/wf.mat',waveformdir));

    a09i0 = load(sprintf('%s/K/plunge/th0_90/a09/wf.mat',waveformdir));
    a09i30 = load(sprintf('%s/K/plunge/th0_60/a09/wf.mat',waveformdir));
    a09i45 = load(sprintf('%s/K/plunge/th0_45/a09/wf.mat',waveformdir));
    a09i60 = load(sprintf('%s/K/plunge/th0_30/a09/wf.mat',waveformdir));
    a09i120 = load(sprintf('%s/K/plunge/th0_30/a-09/wf.mat',waveformdir));
    a09i135 = load(sprintf('%s/K/plunge/th0_45/a-09/wf.mat',waveformdir));
    a09i150 = load(sprintf('%s/K/plunge/th0_60/a-09/wf.mat',waveformdir));
    a09i180 = load(sprintf('%s/K/plunge/th0_90/a-09/wf.mat',waveformdir));

    configs = [a00i0 a00i30 a00i45 a00i60];%...
                %a02i0 a02i30 a02i45 a02i60 a02i120 a02i135 a02i150 a02i180 ...
                %a05i0 a05i30 a05i45 a05i60 a05i120 a05i135 a05i150 a05i180 ...
                %a09i0 a09i30 a09i45 a09i60 a09i120 a09i135 a09i150 a09i180];
    f1 = figure;
    Pos1 = f1.Position;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlabel('$t-t_{\rm LR}$','FontSize',labels_fontsize,'Interpreter','latex')
    ylabel('$\mathcal{A}$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([-230,100])
    ylim([0,1])
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')

    C = gray(4);
    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(5-i,:));
    end

    configs = [a02i0 a02i30 a02i45 a02i60 a02i120 a02i135 a02i150 a02i180];
    f2 = figure;
    f2.Position(1) = Pos1(1) - Pos1(3) - 50;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlabel('$t-t_{\rm LR}$','FontSize',labels_fontsize,'Interpreter','latex')
    ylabel('$\mathcal{A}$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([-350,100])
    ylim([.3,.7])
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')

    C1 = autumn(4);
    C2 = winter(4);
    C = [C2;C1];
    %C = spring(8);
    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');
        %if i==length(configs)
        %    plot(t-tLR,a,'LineWidth',my_linewidth,'LineStyle','--','DisplayName',label,'Color',C(length(configs)+1-i,:));
        %    break
        %end  
        
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end
    ax1 = axes('Position',[0.39,0.48,0.32,0.38]);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','LineWidth',0.5,'box','on');
    xlim([0,100])
    hold on
    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');
        %if i==length(configs)
        %    plot(t-tLR,a,'LineWidth',my_linewidth,'LineStyle','--','DisplayName',label,'Color',C(length(configs)+1-i,:));
        %    break
        %end  
        
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end

    configs = [a05i0 a05i30 a05i45 a05i60 a05i120 a05i135 a05i150 a05i180];
    f3 = figure;
    f3.Position(1) = Pos1(1) + Pos1(3) + 50;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlabel('$t-t_{\rm LR}$','FontSize',labels_fontsize,'Interpreter','latex')
    ylabel('$\mathcal{A}$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([-350,100])
    %ylim([0,1])
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')

    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');  
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end
    
    ax2 = axes('Position',[0.4,0.56,0.31,0.31]);
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','LineWidth',0.5,'box','on');
    xlim([-250,-200])
    hold on
    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');  
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end

    configs = [a09i0 a09i30 a09i45 a09i60 a09i120 a09i135 a09i150 a09i180];
    f4 = figure;
    f4.Position(2) = Pos1(2) - Pos1(4) - 50;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlabel('$t-t_{\rm LR}$','FontSize',labels_fontsize,'Interpreter','latex')
    ylabel('$\mathcal{A}$','FontSize',labels_fontsize,'Interpreter','latex')
    xlim([-350,100])
    %ylim([0,1])
    xline(0,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')

    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end

    ax3 = axes('Position',[0.39,0.59,0.31,0.29]);
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize/2,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    xlim([-250,-200])
    for i=1:length(configs)
        s = configs(i).s;
        label = KerrLabel(s.dyn);
        tLR = tLR_splined(s.dyn);

        [a,t] = DB_asymmetry(s,'nodyn');
        plot(t-tLR,a,'LineWidth',my_linewidth,'DisplayName',label,'Color',C(length(configs)+1-i,:));
    end

return



