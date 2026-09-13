function DB_plot_amplitude_freq(wf1,wf2,varargin)
    % ===================================================================================================
    % varargin: ellmode, mmode, alignement criterion (to compute arbitrary timeshift)
    % ===================================================================================================

    usage_flag = 'pdf';
    insets=0;
    l = 2;
    m = 2;
    align = 'none';

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'l'
                    i      = i + 1;
                    l = varargin{i};
                case 'm'
                    i      = i + 1;
                    m   = varargin{i};
                case 'align'
                    i      = i + 1;
                    align   = varargin{i};
            end
        end
        i = i + 1;
    end

    h1 = wf1.ell(l).emm(m+1).hlm;
    t1 = wf1.ell(l).emm(m+1).t;
    h2 = wf2.ell(l).emm(m+1).hlm;
    t2 = wf2.ell(l).emm(m+1).t;

    iota1 =     90 - rad2deg(wf1.dyn.th0);
    iota2 = 90 - rad2deg(wf2.dyn.th0);

    if strcmp(usage_flag,'pdf')
        my_fontsize = 14;
        my_linewidth = 1;
    elseif strcmp(usage_flag,'slides')
        my_fontsize = 18;
        my_linewidth = 2;
    end

    f1 = figure;
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex')
    fontsize(legend,my_fontsize,'points')
    xlabel('$u/M$','Interpreter','latex')
    ylabel('$|h_{22}|/\nu$','Interpreter','latex')
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times','box','on');

    f2 = figure;
    f2.Position(3:4) = f1.Position(3:4);
    legend(BackgroundAlpha=.7,Location='northwest',Interpreter='latex')
    fontsize(legend,my_fontsize,'points')
    xlabel('$u/M$','Interpreter','latex')
    ylabel('$\omega_{22}$','Interpreter','latex')
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times','box','on');

    [~,t_peak] = Apeak(h1,t1);

    if strcmp(align,'none')
        disp('Alignment criterion not specified: setting dt = 0.');
        dt = 0;
        t_align = 0;
    elseif strcmp(align,'peak')
        dt = find_deltat(h1,t1,h2,t2);
        t_align = t_peak;
    elseif strcmp(align,'LR')
        [t1LR,t2LR,dt] = tLR_splined(wf1.dyn,wf2.dyn);
        t_align = t1LR;
    end    

    %t_start = max(t1(1),t2(1)-dt)+20;
    t_start = t_peak-200;
    t_end = t_peak+100;

    label1 = sprintf('$\\iota = %d^\\circ$, \\ a = %.2f ',iota1,wf1.dyn.chi1(3));
    label2 = sprintf('$\\iota = %d^\\circ$, \\ a = %.2f ',iota2,wf2.dyn.chi1(3));
    %label1 = 'KOS';
    %label2 = 'DB';

    figure(f1);
    ax1 = gca;
    hold on
    plot(t1,abs(h1),'DisplayName',label1,'LineWidth',my_linewidth);
    plot(t2 - dt,abs(h2),'DisplayName',label2,'LineWidth',my_linewidth,'LineStyle','-');
    xlim([t_start,t_end])
    ylim auto
    %ylim([0,0.55])

    if insets==1
        %x1_inset = t_start + 100;
        %x2_inset = t_start + 350;
        x1_inset = t_peak - 30;
        x2_inset = t_peak + 30;
        
        x_length = x2_inset - x1_inset;
        ax2 = axes('Position',[.25 .22 .3 .3]);
        box on
        hold on
        xlim([x1_inset,x2_inset])
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(my_fontsize/2)+1,'FontName','Times');
        plot(t1,abs(h1),'DisplayName',label1,'LineWidth',my_linewidth);
        plot(t2 - dt,abs(h2),'DisplayName',label2,'LineWidth',my_linewidth);
        ylim auto
        %ylim([0.6,0.79])

        y_rect = ylim;
        rectangle('Position', [x1_inset,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax1);
    end

    figure(f2)
    ax3 = gca;
    hold on
    plot(t1,freq(h1,t1),'DisplayName',label1,'LineWidth',my_linewidth)
    plot(t2 - dt,freq(h2,t2),'DisplayName',label2,'LineWidth',my_linewidth,'LineStyle','-');
    hold off
    xlim([t_start,t_end])
    ylim auto
    %ylim([0.15,0.5])
    
    if insets==1
        ax4 = axes('Position',[.25 .42 .3 .3]);
        box on
        hold on

        x1_inset = t_end - 105;
        x2_inset = t_end - 10;
        x_length = x2_inset - x1_inset; 

        xlim([x1_inset,x2_inset])
        %ylim([y1_rect,y2_rect])
        plot(t1,freq(h1,t1),'DisplayName',label1,'LineWidth',my_linewidth)
        plot(t2 - dt,freq(h2,t2),'DisplayName',label2,'LineWidth',my_linewidth);
        %xlabel('$t/M$','FontSize',10,'Interpreter','Latex');
        %ylabel('$\Re [h_{22}]/\mu$','FontSize',10,'Interpreter','Latex');
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(my_fontsize/2)+1,'FontName','Times');
        ylim auto

        y_rect = ylim;
        rectangle('Position', [x1_inset,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent',ax3);
    end

    if ~isempty(varargin)
        figure(f1)
        xline(t_align,'LineStyle','--','Color',[.7 .7 .7], 'HandleVisibility', 'off')
        xline(ax1,t_align,'LineStyle','--','Color',[.7 .7 .7], 'HandleVisibility', 'off')


        figure(f2)
        xline(t_align,'LineStyle','--','Color',[.7 .7 .7], 'HandleVisibility', 'off')
        xline(ax3,t_align,'LineStyle','--','Color',[.7 .7 .7], 'HandleVisibility', 'off')
    end
return