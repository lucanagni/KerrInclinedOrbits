function p = DB_TEMP_build_analytical_structure(a,iota,theta,phi,varargin)

    % ============================================================================================================================================
    % Loads configuration structures according to input and creates a struct p containing the modes computed analytically with the multipoles. 
    % Useful to compute hplus and hcross analytically
    % TEMP because ideally we would like to compute the modes directly with the dynamics (as in KOS) but this is not the case yet
    % ============================================================================================================================================

    if isempty(varargin)
        p = struct;

        for l=2:4
            for m=-l:l
                p = DB_waveforms_test(a,iota,l,m,'output',p);
                close
            end
        end
    else
        p = varargin{1};
    end

    mydir = finddir(a,iota);          
    dir = sprintf('~/waveforms/K/plunge/%s/wf.mat',mydir);
    load(dir)

    [h_num,t_num] = DB_hpc(s,theta,phi);
    [h_an,t_an] = DB_hpc(p,theta,phi);
    DB = s.dyn;

    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    T = t_an-tLR_splined(DB);
    sh_num= spline(t_num,h_num,t_an);
    A_num = abs(sh_num);
    A_an = abs(h_an);

    DeltaA = (A_an-A_num)./A_num;
    DeltaPhi = -unwrap(angle(h_an)) + unwrap(angle(sh_num));
    M = mean(DeltaPhi(200:end-200));
    while abs(M)>pi/2
        if M>0
            DeltaPhi = DeltaPhi - pi;
        else
            DeltaPhi = DeltaPhi + pi;
        end
        M = mean(DeltaPhi(200:end-200));
    end

    
    f = figure;
    f.Position(3:4) = f.Position(3:4)*1.2;
    tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
    t_in = -tLR_splined(DB)+200;
    t_end = 0;
    label = KerrLabel(DB);
    wf_label = 'boh';

    nexttile(tl,[2,1])
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
    xlim([t_in,t_end])
    hold on
    plot(T,real(sh_num),'LineWidth',my_linewidth,'DisplayName',sprintf('$h_+^{\\rm exact}$'))
    plot(T,real(h_an),'LineWidth',my_linewidth,'DisplayName',sprintf('$h_+^{\\rm EOB}$'),'LineStyle','--') 
    plot(T,A_num,'LineWidth',my_linewidth,'DisplayName',sprintf('$|h|^{\\rm exact}$'),'Color',[1 0 0 .2])
    plot(T,A_an,'LineWidth',my_linewidth,'DisplayName',sprintf('$|h|^{\\rm exact}$'),'LineStyle','--','Color',[1 0 0 .2])
    ys = ylim;
    xs = xlim;
    text(xs(1)+abs(xs(2)-xs(1))./3,ys(2) - abs(ys(2)-ys(1))./8,label,'Interpreter','latex','FontSize',16)

    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlim([t_in,t_end])
    hold on

    %if flags.logscale_switch==1
    %    plot(T,abs(DeltaPhi),'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
    %    plot(T,abs(DeltaA),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
%
    %    yscale log
    %    yticks([1e-4,1e-3,1e-2,1e-1,1])
    %    grid off
    %    grid on
    %else
        plot(T,DeltaPhi,'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
        plot(T,DeltaA,'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')

        xlim([t_in+50,t_end-100])
        ylim auto
        ys = ylim;
        xlim([t_in,t_end])
        ylim(ys*1.5)
        grid on
    %end
        %yticks(-0.5:0.25:0.5)

    xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')

return

function dir = finddir(a,iota)
    if iota>90
        sign = '-';
    else
        sign = '';
    end
    th = abs(90-iota);
    th0 = num2str(th);

    dir = sprintf('th0_%s/a%s0%d',th0,sign,a);
return