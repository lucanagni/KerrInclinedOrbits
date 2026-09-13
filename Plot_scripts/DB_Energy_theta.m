function DB_Energy_theta(j,varargin)
    ellmax = 4;

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'ellmax'
                    i      = i + 1;
                    ellmax = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    load(sprintf('~/waveforms/K/plunge/th0_90/a%s/wf.mat',j));
    s0 = s;
    load(sprintf('~/waveforms/K/plunge/th0_60/a%s/wf.mat',j));
    s30 = s;
    load(sprintf('~/waveforms/K/plunge/th0_45/a%s/wf.mat',j));
    s45 = s;
    load(sprintf('~/waveforms/K/plunge/th0_30/a%s/wf.mat',j));
    s60 = s;
    load(sprintf('~/waveforms/K/plunge/th0_30/a-%s/wf.mat',j));
    s120 = s;
    load(sprintf('~/waveforms/K/plunge/th0_45/a-%s/wf.mat',j));
    s135 = s;
    load(sprintf('~/waveforms/K/plunge/th0_60/a-%s/wf.mat',j));
    s150 = s;
    load(sprintf('~/waveforms/K/plunge/th0_90/a-%s/wf.mat',j));
    s180 = s;

    spin = s0.dyn.chi1(3);

    [D_0,Dlm_0] = integrate(s0,ellmax);
    [D_30,Dlm_30] = integrate(s30,ellmax);
    [D_45,Dlm_45] = integrate(s45,ellmax);
    [D_60,Dlm_60] = integrate(s60,ellmax);
    [D_120,Dlm_120] = integrate(s120,ellmax);
    [D_135,Dlm_135] = integrate(s135,ellmax);
    [D_150,Dlm_150] = integrate(s150,ellmax);
    [D_180,Dlm_180] = integrate(s180,ellmax);

    l2(1,:) = [Dlm_0.ell(2).emm(3).DeltaE/D_0, Dlm_30.ell(2).emm(3).DeltaE/D_30, Dlm_45.ell(2).emm(3).DeltaE/D_45, Dlm_60.ell(2).emm(3).DeltaE/D_60, NaN, Dlm_120.ell(2).emm(3).DeltaE/D_120, Dlm_135.ell(2).emm(3).DeltaE/D_135, Dlm_150.ell(2).emm(3).DeltaE/D_150, Dlm_180.ell(2).emm(3).DeltaE/D_180];
    l2(2,:) = [Dlm_0.ell(2).emm(2).DeltaE/D_0, Dlm_30.ell(2).emm(2).DeltaE/D_30, Dlm_45.ell(2).emm(2).DeltaE/D_45, Dlm_60.ell(2).emm(2).DeltaE/D_60, NaN, Dlm_120.ell(2).emm(2).DeltaE/D_120, Dlm_135.ell(2).emm(2).DeltaE/D_135, Dlm_150.ell(2).emm(2).DeltaE/D_150, Dlm_180.ell(2).emm(2).DeltaE/D_180];
    l2(3,:) = [Dlm_0.ell(2).emm(1).DeltaE/D_0, Dlm_30.ell(2).emm(1).DeltaE/D_30, Dlm_45.ell(2).emm(1).DeltaE/D_45, Dlm_60.ell(2).emm(1).DeltaE/D_60, NaN, Dlm_120.ell(2).emm(1).DeltaE/D_120, Dlm_135.ell(2).emm(1).DeltaE/D_135, Dlm_150.ell(2).emm(1).DeltaE/D_150, Dlm_180.ell(2).emm(1).DeltaE/D_180];
    l2(4,:) = [Dlm_0.ell(2).emminus(2).DeltaE/D_0, Dlm_30.ell(2).emminus(2).DeltaE/D_30, Dlm_45.ell(2).emminus(2).DeltaE/D_45, Dlm_60.ell(2).emminus(2).DeltaE/D_60, NaN, Dlm_120.ell(2).emminus(2).DeltaE/D_120, Dlm_135.ell(2).emminus(2).DeltaE/D_135, Dlm_150.ell(2).emminus(2).DeltaE/D_150, Dlm_180.ell(2).emminus(2).DeltaE/D_180];
    l2(5,:) = [Dlm_0.ell(2).emminus(3).DeltaE/D_0, Dlm_30.ell(2).emminus(3).DeltaE/D_30, Dlm_45.ell(2).emminus(3).DeltaE/D_45, Dlm_60.ell(2).emminus(3).DeltaE/D_60, NaN, Dlm_120.ell(2).emminus(3).DeltaE/D_120, Dlm_135.ell(2).emminus(3).DeltaE/D_135, Dlm_150.ell(2).emminus(3).DeltaE/D_150, Dlm_180.ell(2).emminus(3).DeltaE/D_180];
    if ellmax>2
        l3(1,:) = [Dlm_0.ell(3).emm(4).DeltaE/D_0, Dlm_30.ell(3).emm(4).DeltaE/D_30, Dlm_45.ell(3).emm(4).DeltaE/D_45, Dlm_60.ell(3).emm(4).DeltaE/D_60, NaN, Dlm_120.ell(3).emm(4).DeltaE/D_120, Dlm_135.ell(3).emm(4).DeltaE/D_135, Dlm_150.ell(3).emm(4).DeltaE/D_150, Dlm_180.ell(3).emm(4).DeltaE/D_180];
        l3(2,:) = [Dlm_0.ell(3).emm(3).DeltaE/D_0, Dlm_30.ell(3).emm(3).DeltaE/D_30, Dlm_45.ell(3).emm(3).DeltaE/D_45, Dlm_60.ell(3).emm(3).DeltaE/D_60, NaN, Dlm_120.ell(3).emm(3).DeltaE/D_120, Dlm_135.ell(3).emm(3).DeltaE/D_135, Dlm_150.ell(3).emm(3).DeltaE/D_150, Dlm_180.ell(3).emm(3).DeltaE/D_180];
        l3(3,:) = [Dlm_0.ell(3).emm(2).DeltaE/D_0, Dlm_30.ell(3).emm(2).DeltaE/D_30, Dlm_45.ell(3).emm(2).DeltaE/D_45, Dlm_60.ell(3).emm(2).DeltaE/D_60, NaN, Dlm_120.ell(3).emm(2).DeltaE/D_120, Dlm_135.ell(3).emm(2).DeltaE/D_135, Dlm_150.ell(3).emm(2).DeltaE/D_150, Dlm_180.ell(3).emm(2).DeltaE/D_180];
        l3(4,:) = [Dlm_0.ell(3).emm(1).DeltaE/D_0, Dlm_30.ell(3).emm(1).DeltaE/D_30, Dlm_45.ell(3).emm(1).DeltaE/D_45, Dlm_60.ell(3).emm(1).DeltaE/D_60, NaN, Dlm_120.ell(3).emm(1).DeltaE/D_120, Dlm_135.ell(3).emm(1).DeltaE/D_135, Dlm_150.ell(3).emm(1).DeltaE/D_150, Dlm_180.ell(3).emm(1).DeltaE/D_180];
        l3(5,:) = [Dlm_0.ell(3).emminus(2).DeltaE/D_0, Dlm_30.ell(3).emminus(2).DeltaE/D_30, Dlm_45.ell(3).emminus(2).DeltaE/D_45, Dlm_60.ell(3).emminus(2).DeltaE/D_60, NaN, Dlm_120.ell(3).emminus(2).DeltaE/D_120, Dlm_135.ell(3).emminus(2).DeltaE/D_135, Dlm_150.ell(3).emminus(2).DeltaE/D_150, Dlm_180.ell(3).emminus(2).DeltaE/D_180];
        l3(6,:) = [Dlm_0.ell(3).emminus(3).DeltaE/D_0, Dlm_30.ell(3).emminus(3).DeltaE/D_30, Dlm_45.ell(3).emminus(3).DeltaE/D_45, Dlm_60.ell(3).emminus(3).DeltaE/D_60, NaN, Dlm_120.ell(3).emminus(3).DeltaE/D_120, Dlm_135.ell(3).emminus(3).DeltaE/D_135, Dlm_150.ell(3).emminus(3).DeltaE/D_150, Dlm_180.ell(3).emminus(3).DeltaE/D_180];
        l3(7,:) = [Dlm_0.ell(3).emminus(4).DeltaE/D_0, Dlm_30.ell(3).emminus(4).DeltaE/D_30, Dlm_45.ell(3).emminus(4).DeltaE/D_45, Dlm_60.ell(3).emminus(4).DeltaE/D_60, NaN, Dlm_120.ell(3).emminus(4).DeltaE/D_120, Dlm_135.ell(3).emminus(4).DeltaE/D_135, Dlm_150.ell(3).emminus(4).DeltaE/D_150, Dlm_180.ell(3).emminus(4).DeltaE/D_180];
    end

    if ellmax>3
        l4(1,:) = [Dlm_0.ell(4).emm(5).DeltaE/D_0, Dlm_30.ell(4).emm(5).DeltaE/D_30, Dlm_45.ell(4).emm(5).DeltaE/D_45, Dlm_60.ell(4).emm(5).DeltaE/D_60, NaN, Dlm_120.ell(4).emm(5).DeltaE/D_120, Dlm_135.ell(4).emm(5).DeltaE/D_135, Dlm_150.ell(4).emm(5).DeltaE/D_150, Dlm_180.ell(4).emm(5).DeltaE/D_180];
        l4(2,:) = [Dlm_0.ell(4).emm(4).DeltaE/D_0, Dlm_30.ell(4).emm(4).DeltaE/D_30, Dlm_45.ell(4).emm(4).DeltaE/D_45, Dlm_60.ell(4).emm(4).DeltaE/D_60, NaN, Dlm_120.ell(4).emm(4).DeltaE/D_120, Dlm_135.ell(4).emm(4).DeltaE/D_135, Dlm_150.ell(4).emm(4).DeltaE/D_150, Dlm_180.ell(4).emm(4).DeltaE/D_180];
        l4(3,:) = [Dlm_0.ell(4).emm(3).DeltaE/D_0, Dlm_30.ell(4).emm(3).DeltaE/D_30, Dlm_45.ell(4).emm(3).DeltaE/D_45, Dlm_60.ell(4).emm(3).DeltaE/D_60, NaN, Dlm_120.ell(4).emm(3).DeltaE/D_120, Dlm_135.ell(4).emm(3).DeltaE/D_135, Dlm_150.ell(4).emm(3).DeltaE/D_150, Dlm_180.ell(4).emm(3).DeltaE/D_180];
        l4(4,:) = [Dlm_0.ell(4).emm(2).DeltaE/D_0, Dlm_30.ell(4).emm(2).DeltaE/D_30, Dlm_45.ell(4).emm(2).DeltaE/D_45, Dlm_60.ell(4).emm(2).DeltaE/D_60, NaN, Dlm_120.ell(4).emm(2).DeltaE/D_120, Dlm_135.ell(4).emm(2).DeltaE/D_135, Dlm_150.ell(4).emm(2).DeltaE/D_150, Dlm_180.ell(4).emm(2).DeltaE/D_180];
        l4(5,:) = [Dlm_0.ell(4).emm(1).DeltaE/D_0, Dlm_30.ell(4).emm(1).DeltaE/D_30, Dlm_45.ell(4).emm(1).DeltaE/D_45, Dlm_60.ell(4).emm(1).DeltaE/D_60, NaN, Dlm_120.ell(4).emm(1).DeltaE/D_120, Dlm_135.ell(4).emm(1).DeltaE/D_135, Dlm_150.ell(4).emm(1).DeltaE/D_150, Dlm_180.ell(4).emm(1).DeltaE/D_180];
        l4(6,:) = [Dlm_0.ell(4).emminus(2).DeltaE/D_0, Dlm_30.ell(4).emminus(2).DeltaE/D_30, Dlm_45.ell(4).emminus(2).DeltaE/D_45, Dlm_60.ell(4).emminus(2).DeltaE/D_60, NaN, Dlm_120.ell(4).emminus(2).DeltaE/D_120, Dlm_135.ell(4).emminus(2).DeltaE/D_135, Dlm_150.ell(4).emminus(2).DeltaE/D_150, Dlm_180.ell(4).emminus(2).DeltaE/D_180];
        l4(7,:) = [Dlm_0.ell(4).emminus(3).DeltaE/D_0, Dlm_30.ell(4).emminus(3).DeltaE/D_30, Dlm_45.ell(4).emminus(3).DeltaE/D_45, Dlm_60.ell(4).emminus(3).DeltaE/D_60, NaN, Dlm_120.ell(4).emminus(3).DeltaE/D_120, Dlm_135.ell(4).emminus(3).DeltaE/D_135, Dlm_150.ell(4).emminus(3).DeltaE/D_150, Dlm_180.ell(4).emminus(3).DeltaE/D_180];
        l4(8,:) = [Dlm_0.ell(4).emminus(4).DeltaE/D_0, Dlm_30.ell(4).emminus(4).DeltaE/D_30, Dlm_45.ell(4).emminus(4).DeltaE/D_45, Dlm_60.ell(4).emminus(4).DeltaE/D_60, NaN, Dlm_120.ell(4).emminus(4).DeltaE/D_120, Dlm_135.ell(4).emminus(4).DeltaE/D_135, Dlm_150.ell(4).emminus(4).DeltaE/D_150, Dlm_180.ell(4).emminus(4).DeltaE/D_180];
        l4(9,:) = [Dlm_0.ell(4).emminus(5).DeltaE/D_0, Dlm_30.ell(4).emminus(5).DeltaE/D_30, Dlm_45.ell(4).emminus(5).DeltaE/D_45, Dlm_60.ell(4).emminus(5).DeltaE/D_60, NaN, Dlm_120.ell(4).emminus(5).DeltaE/D_120, Dlm_135.ell(4).emminus(5).DeltaE/D_135, Dlm_150.ell(4).emminus(5).DeltaE/D_150, Dlm_180.ell(4).emminus(5).DeltaE/D_180];
    end

    iotas = [0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi];

    C0 = 'k';
    C1 = 'r';
    C2 = 'b';
    C3 = [0.4940 0.1840 0.5560];
    C4 = 'g';

    C0 = 'k';
    C1 = MyColors('r1');
    C2 = MyColors('b1');
    C3 = MyColors('r2');
    C4 = MyColors('b2');

    lbl = KerrLabel(s.dyn);

    f1 = figure;
    f1.Position(3:4) = f1.Position(3:4)*1.1;
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',14,'FontName','Times','LineWidth',0.5,'box','on','TickLabelInterpreter','latex');
    hold on
    plot(iotas,l2(1,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C2,'DisplayName',sprintf('$(2,2)$'))
    plot(iotas,l2(2,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C1,'DisplayName',sprintf('$(2,1)$'))
    plot(iotas,l2(3,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C0,'DisplayName',sprintf('$(2,0)$'))
    plot(iotas,l2(4,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C1,'DisplayName',sprintf('$(2,-1)$'))
    plot(iotas,l2(5,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C2,'DisplayName',sprintf('$(2,-2)$'))
    xlabel('$\iota$','FontSize',18,'Interpreter','Latex');
    xticks(iotas);
    xticklabels({'0' '$\frac{\pi}{6}$' '$\frac{\pi}{4}$' '$\frac{\pi}{3}$' '$\frac{\pi}{2}$' '$\frac{2\pi}{3}$' '$\frac{3\pi}{4}$' '$\frac{5\pi}{6}$' '$\pi$'})
    %xtickangle(45)
    ylabel('$\left[\Delta E_{\ell m}/\Delta E\right]^\mathrm{LSSO}$','FontSize',18,'Interpreter','Latex');
    legend('Location','best','Interpreter','latex','FontSize',12,'BackgroundAlpha',0.6)
    lims = ylim;
    text(pi/12,0.92*lims(2),sprintf('$a=%.1f$',spin),'FontSize',14,'Interpreter','latex')

    if ellmax>2
        f2 = figure;
        f2.Position(3:4)  = f1.Position(3:4);
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',14,'FontName','Times','LineWidth',0.5,'box','on','TickLabelInterpreter','latex');
        hold on
        plot(iotas,l3(1,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C3,'DisplayName',sprintf('$(3,3)$'))
        plot(iotas,l3(2,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C2,'DisplayName',sprintf('$(3,2)$'))
        plot(iotas,l3(3,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C1,'DisplayName',sprintf('$(3,1)$'))
        plot(iotas,l3(4,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C0,'DisplayName',sprintf('$(3,0)$'))
        plot(iotas,l3(5,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C1,'DisplayName',sprintf('$(3,-1)$'))
        plot(iotas,l3(6,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C2,'DisplayName',sprintf('$(3,-2)$'))
        plot(iotas,l3(7,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C3,'DisplayName',sprintf('$(3,-3)$'))
        xlabel('$\iota$','FontSize',18,'Interpreter','Latex');
        xticks(iotas);
        xticklabels({'0' '$\frac{\pi}{6}$' '$\frac{\pi}{4}$' '$\frac{\pi}{3}$' '$\frac{\pi}{2}$' '$\frac{2\pi}{3}$' '$\frac{3\pi}{4}$' '$\frac{5\pi}{6}$' '$\pi$'})
        %xtickangle(45)
        ylabel('$\left[\Delta E_{\ell m}/\Delta E\right]^\mathrm{LSSO}$','FontSize',18,'Interpreter','Latex');
        legend('Location','best','Interpreter','latex','FontSize',12,'BackgroundAlpha',0.6)
        lims = ylim;
        text(pi/12,0.92*lims(2),sprintf('$a=%.1f$',spin),'FontSize',14,'Interpreter','latex')
    end

    if ellmax>3
        f3 = figure;
        f3.Position(3:4)  = f1.Position(3:4);
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',14,'FontName','Times','LineWidth',0.5,'box','on','TickLabelInterpreter','latex');
        hold on
        plot(iotas,l4(1,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C4,'DisplayName',sprintf('$(4,4)$'))
        plot(iotas,l4(2,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C3,'DisplayName',sprintf('$(4,3)$'))
        plot(iotas,l4(3,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C2,'DisplayName',sprintf('$(4,2)$'))
        plot(iotas,l4(4,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C1,'DisplayName',sprintf('$(4,1)$'))
        plot(iotas,l4(5,:),'Marker','square','LineWidth',1.5,'LineStyle','-','Color',C0,'DisplayName',sprintf('$(4,0)$'))
        plot(iotas,l4(6,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C1,'DisplayName',sprintf('$(4,-1)$'))
        plot(iotas,l4(7,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C2,'DisplayName',sprintf('$(4,-2)$'))
        plot(iotas,l4(8,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C3,'DisplayName',sprintf('$(4,-3)$'))
        plot(iotas,l4(9,:),'Marker','square','LineWidth',1.5,'LineStyle','--','Color',C4,'DisplayName',sprintf('$(4,-4)$'))
        xlabel('$\iota$','FontSize',18,'Interpreter','Latex');
        xticks(iotas);
        xticklabels({'0' '$\frac{\pi}{6}$' '$\frac{\pi}{4}$' '$\frac{\pi}{3}$' '$\frac{\pi}{2}$' '$\frac{2\pi}{3}$' '$\frac{3\pi}{4}$' '$\frac{5\pi}{6}$' '$\pi$'})
        %xtickangle(45)
        ylabel('$\left[\Delta E_{\ell m}/\Delta E\right]^\mathrm{LSSO}$','FontSize',18,'Interpreter','Latex');
        legend('Location','best','Interpreter','latex','FontSize',12,'BackgroundAlpha',0.6)
        lims = ylim;
        text(pi/12,0.92*lims(2),sprintf('$a=%.1f$',spin),'FontSize',14,'Interpreter','latex')
    end

return

function [DeltaE_LSSO,DeltaE_lm] = integrate(s,ellmax)
    dyn = s.dyn;
    [dEdt,t,dEdt_lm] = DB_dEdt(s,'ellmax',ellmax);

    tLSSO = tLSSO_splined(dyn);

    idx0 = find(t>=tLSSO-10, 1, 'first');
    t_dense = t(idx0):0.01:t(end);

    dEdt_splined = spline(t,dEdt,t_dense);
    idx1 = find(t_dense>=tLSSO,1,'first');

    DeltaE_vec = cumtrapz(t_dense(idx1:end),dEdt_splined(idx1:end));

    DeltaE_LSSO = DeltaE_vec(end);

    for l=2:ellmax
        for m=-l:l
            if m<0
                dEdt_lm_0 = dEdt_lm.ell(l).emminus(-m+1).flux;
            else
                dEdt_lm_0 = dEdt_lm.ell(l).emm(m+1).flux;
            end
            dEdt_lm_splined = spline(t,dEdt_lm_0,t_dense);
            DeltaE_lm_vec = cumtrapz(t_dense(idx1:end),dEdt_lm_splined(idx1:end));
            if m<0
                DeltaE_lm.ell(l).emminus(-m+1).DeltaE = DeltaE_lm_vec(end);
            else
                DeltaE_lm.ell(l).emm(m+1).DeltaE = DeltaE_lm_vec(end);
            end
        end
    end


return
    