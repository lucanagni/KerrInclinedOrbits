function Twist_plot(spin)
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    l = 2;
    norm = sqrt((l+2)*(l+1)*l*(l-1));

    if spin==0.2
        load('~/waveforms/K/plunge/th0_60/a023/wf.mat');
        sp = s;
        load('~/waveforms/K/plunge/th0_90/a02/wf.mat');
        se = s;

        eq_label = '$\Psi_{22}$ - {\tt a02i0}';
        rot_label = '$\tilde\Psi_{22}$ - {\tt a023i30}';
        prec_label = '$\Psi_{22}$ - {\tt a023i30}';
    elseif spin==0.5
        load('~/waveforms/K/plunge/th0_45/a05/wf.mat');
        sp = s;
        load('~/waveforms/K/twist/i0/wf.mat');
        se = s;

        load('~/waveforms/K/plunge/th0_90/a05/wf.mat');
        hextra = s.ell(2).emm(3).hlm./norm;
        textra = s.ell(2).emm(3).t;
        
        eq_label = '$\Psi_{22}$ - {\tt a035i0}';
        rot_label = '$\tilde\Psi_{22}$ - {\tt a05i45}';
    else
        error('spin can be either 0.2 or 0.5')
    end
    

    if spin==0.2
        
    elseif spin==0.5
        
    else
        error('spin can be either 0.2 or 0.5')
    end

    

    sr = struct;
    sr = DB_mode_rotate_timedep(sp,sp.dyn,2,2,sr);

    dt = tLR_splined(sp.dyn) - tLR_splined(se.dyn);
    if spin==0.5
        dt2 = tLR_splined(s.dyn) - tLR_splined(se.dyn);
    end
    tLR = tLR_splined(se.dyn);

    dyn_end = se.dyn.t(end);

    heq = se.ell(2).emm(3).hlm./norm;
    teq = se.ell(2).emm(3).t;
    hr = sr.ell(2).emm(3).hlm./norm;
    tr = sp.ell(2).emm(3).t;

    if spin==0.2
        t_start = tLR - 800;
        t_end = tLR + 130;
    else
        t_start = tLR - 800;
        t_end = tLR + 80;
    end

    f1 = figure;
    f1.Position = [635 295 550 426];
    ax1 = gca;
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(teq,abs(heq),'LineWidth',my_linewidth,'Color',MyColors('b1'),'DisplayName',eq_label)
    plot(tr-dt,abs(hr),'LineWidth',my_linewidth,'Color',MyColors('r1'),'DisplayName',rot_label)
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xline(dyn_end,'Color',[.7 .7 .7],'LineStyle',':','HandleVisibility','off')
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    ylabel('$|\Psi_{22}|/\nu$','FontSize',labels_fontsize,'Interpreter','Latex');
    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
    xlim([t_start,t_end])
    axis manual
    hold on
    if spin==0.5
        plot(textra-dt2,abs(hextra),'LineWidth',my_linewidth,'Color',MyColors('g1'),'DisplayName','$\Psi_{22}$ - {\tt a05i0}')
    else
        plot(tr-dt,abs(sp.ell(2).emm(3).hlm./norm),"Color",[0 0 0 .1], 'LineWidth',my_linewidth,'DisplayName',prec_label)
    end

    x1_rect = t_start + 100;
    x2_rect = t_start + 450;
    x_length = x2_rect - x1_rect;

    ax_inset1 = axes('Position',[.3 .2 .3 .3],'FontName','Times');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(axes_fontsize/2)+1,'FontName','Times','box','on');    
    hold on
    xlim([x1_rect,x2_rect])
    %ylim([y1_rect,y2_rect])
    ylim auto
    plot(teq,abs(heq),'LineWidth',my_linewidth,'Color',MyColors('b1'),'DisplayName',eq_label)
    plot(tr-dt,abs(hr),'LineWidth',my_linewidth,'Color',MyColors('r1'),'DisplayName',rot_label)
    y_rect = ylim;
    rectangle('Position', [x1_rect,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent', ax1);
    axis manual 
    hold on
    if spin==0.5
        plot(textra-dt2,abs(hextra),'LineWidth',my_linewidth,'Color',MyColors('g1'),'DisplayName','$\Psi_{22}$ - {\tt a05i0}')
    else
        plot(tr-dt,abs(sp.ell(2).emm(3).hlm./norm),"Color",[0 0 0 .1], 'LineWidth',my_linewidth,'DisplayName',prec_label)
    end


    f2 = figure;
    f2.Position = [635 295 550 426];
    ax2 = gca;
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(teq,freq(heq,teq),'LineWidth',my_linewidth,'Color',MyColors('b1'),'DisplayName',eq_label)
    plot(tr-dt,freq(hr,tr),'LineWidth',my_linewidth,'Color',MyColors('r1'),'DisplayName',rot_label)

    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xline(dyn_end,'Color',[.7 .7 .7],'LineStyle',':','HandleVisibility','off')
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$\omega_{22}$','FontSize',labels_fontsize,'Interpreter','Latex');
    xlim([t_start,t_end])
    axis manual 
    hold on
    if spin==0.5
        plot(textra-dt2,freq(hextra,textra),'LineWidth',my_linewidth,'Color',MyColors('g1'),'DisplayName','$\Psi_{22}$ - {\tt a05i0}')
        ylim([0.15,0.5])
    else
        plot(tr-dt,freq(sp.ell(2).emm(3).hlm./norm,tr),"Color",[0 0 0 .1], 'LineWidth',my_linewidth,'DisplayName',prec_label)
    end

    ax_inset2 = axes('Position',[.3 .4 .3 .3],'FontName','Times');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(axes_fontsize/2)+1,'FontName','Times','box','on');    
    hold on
    xlim([x1_rect,x2_rect])
    %ylim([y1_rect,y2_rect])
    ylim auto
    plot(teq,freq(heq,teq),'LineWidth',my_linewidth,'Color',MyColors('b1'),'DisplayName',eq_label)
    plot(tr-dt,freq(hr,tr),'LineWidth',my_linewidth,'Color',MyColors('r1'),'DisplayName',rot_label)
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',round(axes_fontsize/2)+1,'FontName','Times');

    y_rect = ylim;
    rectangle('Position', [x1_rect,y_rect(1),x_length,abs(y_rect(2) - y_rect(1))],'EdgeColor', 'k', 'LineWidth', .5, 'Parent', ax2);

    axis manual
    hold on
    if spin==0.5
        plot(textra-dt2,freq(hextra,textra),'LineWidth',my_linewidth,'Color',MyColors('g1'),'DisplayName','$\omega_{22}$ - {\tt a05i0}')
    else
        plot(tr-dt,freq(sp.ell(2).emm(3).hlm./norm,tr),"Color",[0 0 0 .1], 'LineWidth',my_linewidth,'DisplayName',prec_label)
    end

    %{
    hdot = DB_D1(hr,tr,4);
    newfreq = -imag((hdot./hr));

    j = find(tr>sp.dyn.t(end),1);
    j = j-1;
    %tr(j-10:j+10)

    figure
    hold on
    plot(tr-dt,-unwrap(angle(hr)))
    plot(tr-dt,DB_D1(-unwrap(angle(hr)),tr,4))
    plot(tr-dt,newfreq)
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--')
    xline(tr(j)-dt,'Color','r','LineWidth',2,'LineStyle',':')
    xline(tr(j-1)-dt,'Color','b','LineWidth',2,'LineStyle',':')
    xline(tr(j-2)-dt,'Color','m','LineWidth',2,'LineStyle',':')

    legend
    ylim([0.12197      0.45558])
    xlim([779.92       986.18])

    figure
    hold on
    plot(tr-dt,DB_D1(real(hr),tr,2))
    plot(tr-dt,DB_D1(imag(hr),tr,2))
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xline(tr(j)-dt,'Color','r','LineWidth',2,'LineStyle',':')

    figure
    hold on
    plot(tr-dt,real(hr))
    plot(tr-dt,imag(hr))
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
    xline(tr(j)-dt,'Color','r','LineWidth',2,'LineStyle',':')
    %}
return