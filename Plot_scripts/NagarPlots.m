function NagarPlots(l,m)
    % ==========================================================================================================================
    % teq = timesteps for the equatorial waveform
    % heq = equatorail waveform
    % tprec = timesteps for the precssing waveform
    % hprec = precessing waveform (inertial frame)
    % hrot = rotated precessing waveform (coprecessing frame)
    % ==========================================================================================================================
    old = 0; %0 -> single plot containing real part + ampl/phase difference. 1 -> two separate plots
    emminus = 0; %if 1 plots the difference between +-2 modes

    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    load('/home/luca/waveforms/K/geod/eq/config/wf1.mat')
    seq = s;
    %load('/home/luca/waveforms/K/geod/precessing/config/wf1.mat')
    load('/home/luca/waveforms/K/twist/i30/wf1.mat')  %same spin-spin
    sp = s;
    sr = struct;
    sr = DB_mode_rotate_timedep(sp,sp.dyn,l,m,sr,'iter',1); %iter 1 uses L, iter 3 uses L_N
    %for m=-l:l
    %    sr = DB_mode_rotate_timedep(sp,sp.dyn,l,m,sr,'iter',1);
    %end
    %sr = DB_mode_rotate_timedep(sr,sr.dyn,l,m,sr,'iter',2);
    %sr = DB_mode_rotate_timedep(sp,sp.dyn,l,m,sr,'iter',1);

    norm = sqrt((l+2)*(l+1)*l*(l-1));
    
    teq = seq.ell(l).emm(m+1).t;
    heq = seq.ell(l).emm(m+1).hlm/norm;

    tprec = sp.ell(l).emm(m+1).t;
    hprec = sp.ell(l).emm(m+1).hlm/norm;

    hrot = sr.ell(l).emm(m+1).hlm/norm;

    [heq_clean,teq_clean] = CleanJunk(heq,teq);
    [hrot_clean,tprec_clean] = CleanJunk(hrot,tprec);

    DA = abs(abs(heq_clean) - abs(spline(tprec_clean,hrot_clean,teq_clean)))./abs(heq_clean);

    i = 0;
    DP_mean = mean(ph(heq_clean)-spline(tprec_clean,ph(hrot_clean),teq_clean));
    %if DP_mean < 0
    %    while abs(DP_mean) > 2*pi
    %        DP_mean = DP_mean + i*pi;
    %        i = i + 1;
    %    end
    %else
    %    while abs(DP_mean) > 2*pi
    %        DP_mean = DP_mean + i*pi;
    %        i = i - 1;
    %    end
    %end
    
    DP = ph(heq_clean)-spline(tprec_clean,ph(hrot_clean),teq_clean) + i*pi;

    % REAL PART OF EQUATORIAL WAVEFORM
    f1 = figure;
    %f1.Position = [100 560 600 600];
    %tiledlayout(2,1)
    plot(teq,real(heq),'LineWidth',my_linewidth+.5,'DisplayName',' $\iota = 0^\circ,\ a = 0.2$')
    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$\Re [\Psi_{22}]/\nu$','FontSize',labels_fontsize,'Interpreter','Latex');
    set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    %legend('BackgroundAlpha',.7,'FontSize',legend_fontsize,'Interpreter','latex')
    text(400,0.4,'{\tt a02i0}','FontSize',14,'Interpreter','latex')
    ylim([-.5,.5])
    xlim([200,2500])

    % REAL PART OF PRECESSING WAVEFORM
    f2 = figure;
    plot(tprec,real(hprec),'LineWidth',my_linewidth+.5,'DisplayName',' $\iota = 30^\circ,\ a = 0.2309$','Color','#D95319')
    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$\Re [\Psi_{22}]/\nu$','FontSize',labels_fontsize,'Interpreter','Latex');
    set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
    %legend('BackgroundAlpha',.7,'FontSize',legend_fontsize,'Interpreter','latex')
    text(400,0.4,'{\tt a023i30}','FontSize',14,'Interpreter','latex')
    ylim([-.5,.5])
    xlim([200,2500])

    if old
        % REAL PART OF EQUATORIAL AND ROTATED WAVEFORM (WITH INSETS)
        f3 = figure;
        %f2.Position = [100 100 600 600];
        hold on
        box on
        plot(teq,real(heq),'LineWidth',my_linewidth+.5,'DisplayName','$\Psi_{22}$ - {\tt a02i0}')
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'DisplayName','$\tilde\Psi_{22}$ - {\tt a023i30}','LineStyle','--')
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        ylabel('$\Re [\Psi_{22}]/\nu$','FontSize',labels_fontsize,'Interpreter','Latex');
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
        %title('Overlap', 'FontSize',16, 'FontName', 'Times');
        legend('BackgroundAlpha',.7,'FontSize',legend_fontsize,'Interpreter','latex','Location','southeast')
        ylim([-.5,.5])
        xlim([200,2500])

        [peaks,peaks_pos] = findpeaks(real(heq_clean));
        x_length = (teq_clean(peaks_pos(2))-teq_clean(peaks_pos(1)))./3;
        y_length = peaks(1)./13;

        x1_rect = teq_clean(peaks_pos(1)) - x_length./2;
        y1_rect = peaks(1) - .75*y_length;
        rectangle('Position', [x1_rect,y1_rect,x_length,y_length],'EdgeColor', 'k', 'LineWidth', .5);

        x2_rect = teq_clean(peaks_pos(end))-x_length./2;
        y2_rect = peaks(end) - .75*y_length;
        rectangle('Position', [x2_rect, y2_rect,x_length,y_length],'EdgeColor', 'k', 'LineWidth', .5);

        axes('Position',[.23 .72 .18 .18]);
        box on
        hold on
        xlim([x1_rect,x1_rect+x_length])
        ylim([y1_rect,y1_rect+y_length])
        plot(teq,real(heq),'LineWidth',my_linewidth+.5)
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'LineStyle','--')
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
        
        axes('Position',[.7 .72 .18 .18]);
        box on
        hold on
        xlim([x2_rect,x2_rect+x_length])
        ylim([y2_rect,y2_rect+y_length])
        plot(teq,real(heq),'LineWidth',my_linewidth+.5)
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'LineStyle','--')
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
        

        % RELATIVE AMPLITUDE DIFFERENCE AND DEPHASING
        f3 = figure;
        %f3.Position = [100 100 600 600];
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact')

        nexttile
        plot(teq_clean,DA,'DisplayName','Amplitude Relative Difference','LineWidth',my_linewidth,'Color','k')
        %xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        ylabel('$\Delta A/{A}^{eq}$','FontSize',labels_fontsize,'Interpreter','Latex');
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
        xlim([teq_clean(1),teq_clean(end)])
        ylim(ylim + [0, 1] * 0.1 * range(ylim));
        grid on

        nexttile
        plot(teq_clean,DP,'DisplayName','Phase Difference','LineWidth',my_linewidth,'Color','k')
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        ylabel('$\Delta\phi$','FontSize',labels_fontsize,'Interpreter','Latex');
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
        xlim([teq_clean(1),teq_clean(end)])
        grid on
    else
        % REAL PART OF EQUATORIAL AND ROTATED WAVEFORM (WITH INSETS) + ERRORS
        f3 = figure;
        LW = f3.Position(3:4);
        f3.Position(3:4) = LW*1.2;
        tl = tiledlayout(5,1,'TileSpacing','tight','Padding','tight');

        ax = nexttile(tl,[3,1]);
        hold on
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','box','on','FontSize',axes_fontsize,'FontName','Times');
        ylabel('$\Re [\Psi_{22}]/\nu$','FontSize',labels_fontsize,'Interpreter','Latex');
        legend('BackgroundAlpha',.7,'FontSize',legend_fontsize-1,'Interpreter','latex','Location','southeast')
        plot(teq,real(heq),'LineWidth',my_linewidth+.5,'DisplayName','$\Psi_{22}$ - {\tt a02i0}')
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'DisplayName','$\tilde\Psi_{22}$ - {\tt a023i30}','LineStyle','--')
        ylim([-.4,.6])
        xlim([200,2500])

        [peaks,peaks_pos] = findpeaks(real(heq_clean));
        x_length = (teq_clean(peaks_pos(2))-teq_clean(peaks_pos(1)))./3;
        y_length = peaks(1)./13;

        x1_rect = teq_clean(peaks_pos(1)) - x_length./2;
        y1_rect = peaks(1) - .75*y_length;
        rectangle('Position', [x1_rect,y1_rect,x_length,y_length],'EdgeColor', 'k', 'LineWidth', .5);

        x2_rect = teq_clean(peaks_pos(end))-x_length./2;
        y2_rect = peaks(end) - .75*y_length;
        rectangle('Position', [x2_rect, y2_rect,x_length,y_length],'EdgeColor', 'k', 'LineWidth', .5);

        axes('Position',[.21 .8 .18 .16]); %[.76 .78 .18 .18] if tiledlayout(3,1)
        box on
        hold on
        xlim([x1_rect,x1_rect+x_length])
        ylim([y1_rect,y1_rect+y_length])
        plot(teq,real(heq),'LineWidth',my_linewidth+.5)
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'LineStyle','--')
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
        
        axes('Position',[.76 .8 .18 .16]); %[.76 .78 .18 .18] if tiledlayout(3,1)
        box on
        hold on
        xlim([x2_rect,x2_rect+x_length])
        ylim([y2_rect,y2_rect+y_length])
        plot(teq,real(heq),'LineWidth',my_linewidth+.5)
        plot(tprec,real(hrot),'LineWidth',my_linewidth+.5,'LineStyle','--')
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
        
        nexttile(tl)
        hold on
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','box','on');
        %xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        plot(teq_clean,DA,'LineStyle','-','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
        %plot(teq_clean,DP,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)

        xlim([teq_clean(1),teq_clean(end)])
        ylim(ylim + [0, 1] * 0.1 * range(ylim));
        grid on

        nexttile(tl)
        hold on
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','box','on');
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        %plot(teq_clean,DA,'LineStyle','-','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
        plot(teq_clean,DP,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)

        xlim([teq_clean(1),teq_clean(end)])
        ylim(ylim + [0, 1] * 0.1 * range(ylim));
        grid on
    end

    if emminus==1
        sr = DB_mode_rotate_timedep(sp,sp.dyn,l,-m,sr);
        hrot2 = sr.ell(l).emminus(m+1).hlm/norm;
        %[hrot_clean,tprec_clean] = CleanJunk(hrot,tprec);

        % REAL PART OF +2 AND -2 MODES WITH AMPLITUDES (ROTATED WAVEFORM)
        ff = figure;
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times');
        xlim([1000,2000])
        ylim([-.22,.22])
        hold on
        plot(tprec,real(hrot),'LineWidth',my_linewidth,'LineStyle','-','Color',MyColors('r1'),'DisplayName','$\Re[\tilde \Psi_{22}]/\nu$')
        plot(tprec,real(hrot2),'LineWidth',my_linewidth,'LineStyle','--','Color',MyColors('b1'),'DisplayName','$\Re[\tilde \Psi_{2-2}]/\nu$')
        plot(tprec,abs(hrot),'LineWidth',my_linewidth,'LineStyle','-','Color',MyColors('r2'),'DisplayName','$|\tilde \Psi_{22}|/\nu$')
        plot(tprec,abs(hrot2),'LineWidth',my_linewidth,'LineStyle','--','Color',MyColors('b2'),'DisplayName','$|\tilde \Psi_{2-2}|/\nu$')
        legend('BackgroundAlpha',.85,'FontSize',legend_fontsize,'Interpreter','latex','Location','northwest')
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');

        x_rect = 1700;
        x_length = 130;
        y_rect = 0.166;
        y_length = 0.018;
        rectangle('Position', [x_rect, y_rect,x_length,y_length],'EdgeColor', 'k', 'LineWidth', .5);

        axes('Position',[0.25,0.19,0.4,0.4]);
        set(gca,'box','on','XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize/2+1,'FontName','Times');
        xlim([x_rect,x_rect+x_length])
        ylim([y_rect,y_rect+y_length])
        yticks([0.17,0.18])
        xticks([1700,1760,1820])
        hold on
        plot(tprec,real(hrot),'LineWidth',my_linewidth,'LineStyle','-','Color',MyColors('r1'),'DisplayName','$\Re[\tilde \Psi_{22}]/\nu$')
        plot(tprec,real(hrot2),'LineWidth',my_linewidth,'LineStyle','--','Color',MyColors('b1'),'DisplayName','$\Re[\tilde \Psi_{2-2}]/\nu$')
        plot(tprec,abs(hrot),'LineWidth',my_linewidth,'LineStyle','-','Color',MyColors('r2'),'DisplayName','$|\tilde \Psi_{22}|/\nu$')
        plot(tprec,abs(hrot2),'LineWidth',my_linewidth,'LineStyle','--','Color',MyColors('b2'),'DisplayName','$|\tilde \Psi_{2-2}|/\nu$')
        

        A4Width(ff)
    end
return


function x = ph(x)
    x = -unwrap(angle(x));
return

function [h_c,t_c] = CleanJunk(h,t)
    tin_pos = find(t>200,1);
    
    t_c = t(tin_pos:end);
    h_c = h(tin_pos:end);
return