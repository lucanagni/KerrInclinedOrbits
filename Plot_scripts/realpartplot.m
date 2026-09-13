function realpartplot(wf1,wf2,varargin)
% ========================================================================================
% plots real part of waveform nicely
% varargin: ellmode, mmode, alignement criterion (to compute arbitrary timeshift)
% ========================================================================================
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
    h2 = wf2.ell(l).emm(m+1).hlm;
    t1 = wf1.ell(l).emm(m+1).t;
    t2 = wf2.ell(l).emm(m+1).t;
    phi1 = -unwrap(angle(h1));
    phi2 = -unwrap(angle(h2));
    dyn1 = wf1.dyn;
    dyn2 = wf2.dyn;

    name1 = sprintf('$\\iota = %d^\\circ, \\ a = %.2f$',90 - rad2deg(dyn1.th0),dyn1.chi1(3));
    name2 = sprintf('$\\iota = %d^\\circ, \\ a = %.2f$',90 - rad2deg(dyn2.th0),dyn2.chi1(3));
    my_linewidth = 1;
    my_fontsize = 14;

    if strcmp(align,'none')
        disp('Alignment criterion not specified: setting dt = 0.');
        dt  = 0;
        sphi2 = spline(t2-dt,phi2,t1);

        idx0 = find(t1>50,1);
        idx1 = find(t1>100,1);
        [~,idx_align] = max(real(h1(idx0:idx1)));
        idx_align = idx0 + idx_align;

        delta_phi = sphi2(idx_align) - phi1(idx_align);
        h1 = h1.*exp(-1i*delta_phi);
        DeltaPhi = sphi2 - phi1 - delta_phi;

        t_start = max(t2(1)-dt,t1(1));
        t_end = min(t1(end),t2(end)-dt)-100;
    elseif strcmp(align,'LR')
        [t1LR,~,dt] = tLR_splined(dyn1,dyn2);
        sphi2 = spline(t2-dt,phi2,t1);

        idx_align = find(t1>t1LR,1);

        delta_phi = sphi2(idx_align) - phi1(idx_align);
        h1 = h1.*exp(-1i*delta_phi);
        DeltaPhi = sphi2 - phi1 - delta_phi;
    elseif strcmp(align,'peak')
        dt = find_deltat(h1,t1,h2,t2);
        sphi2 = spline(t2-dt,phi2,t1);

        idx_align = find(t1>Apeak(h1,t1),1);

        delta_phi = sphi2(idx_align) - phi1(idx_align);
        h1 = h1.*exp(-1i*delta_phi);
        DeltaPhi = sphi2 - phi1 - delta_phi;
    else 
        error('Alignment criterion not recognized');
    end    
    t_align = t1(idx_align);
    if dt~=0
        t_start = t_align-100;
        t_end = t_align + 150;
    end

    f1 = figure;
    subplot(3,1,[1 2])
    ax1 = gca;
    legend(BackgroundAlpha=.7,Location='northeast',Interpreter='latex')
    fontsize(legend,my_fontsize,'points')
    %xlabel('$u/M$','Interpreter','latex')
    ylabel('$\Re [h_{22}]/\nu$','Interpreter','latex')
    set(ax1,'XMinorTick','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times','box','on');
    xlim([t_start,t_end]);

    subplot(3,1,3)
    ax2 = gca;
    xlabel('$u/M$','Interpreter','latex')
    ylabel('$\Delta \phi_{22}$','Interpreter','latex')
    set(ax2,'XMinorTick','on','YMinorTick','on','FontSize',my_fontsize,'FontName','Times','box','on');
    xlim([t_start,t_end]);
    ylim([-pi/2,pi/2])

    hold(ax1,'on')
    plot(ax1,t1,real(h1),'DisplayName',name1,'LineWidth',my_linewidth)
    plot(ax1,t2-dt,real(h2),'DisplayName',name2,'LineWidth',my_linewidth,'LineStyle','--');    
    xline(ax1,t_align,'Color',[0.7 0.7 0.7],'LineStyle','--','HandleVisibility','off')

    hold(ax2,'on')
    plot(ax2,t1,DeltaPhi,'LineWidth',my_linewidth,'HandleVisibility','off')
    xline(ax2,t_align,'Color',[0.7 0.7 0.7],'LineStyle','--','HandleVisibility','off')

return