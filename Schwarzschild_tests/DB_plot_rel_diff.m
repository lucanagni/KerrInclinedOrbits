function [rel_diff,t1] = DB_plot_rel_diff(heq,teq,hp,tp,mycase)
    % ========================================================================================
    % for PDF fontsize 14 should be fine
    % input are [h,t] for equatorial and precessing dynamics (in theory at least)
    % ========================================================================================


    usage_flag = 'pdf';
    inset = 1;

    if strcmp(usage_flag,'pdf')
        my_fontsize = 12;
        my_linewidth = 2;
    elseif strcmp(usage_flag,'slides')
        my_fontsize = 18;
        my_linewidth = 2;
    end

    % if what I'm plotting has a peak, align at peak
    %[~,peak1] = Apeak(heq,teq);
    %[~,peak2] = Apeak(hp,tp);
    %dt = peak2 - peak1;
    dt = find_deltat(heq,teq,hp,tp);
    %dt = 0;

    tp = tp - dt;

    %spline to the smaller t(end) so that we can work with only one time
    if teq(end)>(tp(end)+dt)
        heq = spline(teq,heq,tp);
        teq = tp;
    else 
        hp = spline(tp,hp,teq);
        tp = teq;
    end
    
        
    %{
    absh = abs(heq);
    peak = max(findpeaks(absh(t_junkend:end)));
    peak_pos = find(absh==peak,1);
    t_start = teq(peak_pos)-1000;
    t_end = teq(peak_pos)+100;
    t_start_pos = find(teq>t_start,1);
    %}

    %t_start = 100;
    %t_end = 400;

    [peak,t_peak] = Apeak(heq,teq);

    t_start = t_peak-200;
    t_end = t_peak+100;

     
    %UNCOMMENT TO GET PHASE DIFFERENCE
        Dphi = -unwrap(angle(hp))+unwrap(angle(heq));
        figure
        plot(teq,Dphi+2*pi,'LineWidth',1)
        xlim([t_start,t_end])
        xlabel('$u/M$','Interpreter','latex')
        %legend('Interpreter','latex','FontSize',my_fontsize,'BackgroundAlpha',0.7,'Location','southeast')
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times');
        ylabel('$|\phi_{eq}-\phi_{tilt}|$','Interpreter','latex')
        %hp = hp.*exp(+1i.*Dphi);
    
    


    if strcmp(mycase,'ampl')
        array1 = abs(heq);
        t1 = teq;
        array2 = abs(hp);
        t2 = tp;
        rel_diff = abs(array2-array1)./array1;
    elseif strcmp(mycase,'freq')
        array1 = freq(heq,teq);
        t1 = teq;
        array2 = freq(hp,tp);
        t2 = tp;
        rel_diff = abs(array2-array1)./array1;
    elseif strcmp(mycase,'real')
        array1 = real(heq);
        t1 = teq;
        array2 = real(hp);
        t2 = tp;
        rel_diff = abs(array2-array1);
    elseif strcmp(mycase,'imag')
        array1 = imag(heq);
        t1 = teq;
        array2 = imag(hp);
        t2 = tp;
        rel_diff = abs(array2-array1);
    else 
        disp('Unrecognised case')
        array1 = heq;
        array2 = spline(tp,hp,teq);
        t1 = teq;
        t2 = t1;
        rel_diff = abs(array2-array1); % ACHTUNG: This is absolute difference!
    end


    figure

    ax1 = subplot(3,1,[1,2]);
    hold on
    plot(t1,array1,'LineWidth',my_linewidth,'DisplayName','$\iota = 0^\circ$','LineStyle','-')
    plot(t2,array2,'LineWidth',my_linewidth,'DisplayName','$\iota = 30^\circ$   ','LineStyle','--')
    xlim([t_start,t_end])
    %xlabel('$u/M$','Interpreter','latex')
    legend('Interpreter','latex','FontSize',my_fontsize+2,'BackgroundAlpha',0.7,'Location','northwest')
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times');
    %ylabel('$|h_{21}|/\nu$','Interpreter','latex')
    ylabel('$\left| d \mathbf{J}/dt \right|$', 'Interpreter', 'latex')
    
    if inset==1
        x1_rect = t_start + 10;
        x2_rect = t_start + 150;
        x_length = x2_rect - x1_rect;

        ax2 = axes('Position',[.28 .6 .3 .3]);
        box on
        hold on
        xlim([x1_rect,x2_rect])
        %ylim([y1_rect,y2_rect])
        ylim auto
        plot(t1,array1,'LineWidth',my_linewidth,'LineStyle','-')
        plot(t2,array2,'LineWidth',my_linewidth,'LineStyle','-')
        %xlabel('$t/M$','FontSize',10,'Interpreter','Latex');
        %ylabel('$\Re [h_{22}]/\mu$','FontSize',10,'Interpreter','Latex');
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(my_fontsize/2)+1,'FontName','Times');


        y_rect = ylim;
        rectangle('Position', [x1_rect,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent', ax1);
    end

    subplot(3,1,3)
    plot(t1,rel_diff,'LineWidth',my_linewidth)
    xlim([t_start,t_end])
    xlabel('$u/M$','Interpreter','latex')
    %ylabel('$\Delta\mathcal{A}/\mathcal{A}^{eq}$','Interpreter','Latex');
    ylabel('$ | \Delta \dot{\mathbf{J}} |/ | \dot{\mathbf{J}_{eq}} |$', 'Interpreter', 'latex')
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times');

    ymax = ylim;
    ylim([0,1.1*ymax(2)])

    grid on

return