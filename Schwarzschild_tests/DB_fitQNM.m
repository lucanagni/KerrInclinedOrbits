function a = DB_fitQNM(s,l,m)
    if m<0
        h = s.ell(l).emminus(-m+1).hlm;
        t = s.ell(l).emminus(-m+1).t;
    else
        h = s.ell(l).emm(m+1).hlm;
        t = s.ell(l).emm(m+1).t;
    end

    omg = freq(h,t);
    tLR = tLR_splined(s.dyn);

    tin = tLR+30;
    tend = tLR+80;

    idx0 = find(t>tin,1);
    idx1 = find(t>tend,1);

    omg_ringdown = omg(idx0:idx1);
    t_ringdown = t(idx0:idx1);

    omega0 = sign(0.5*(2*m-1))*imag(kerr_FitKerrQNMs(l,2,0,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));

    omg_QNM = @(C,t_ringdown)((1-exp(2*C(1))).*omega0)./(1+exp(2*C(1))+2*exp(C(1)).*cos(2*omega0.*t_ringdown+C(2)));

    mdl = fitnlm(t_ringdown,omg_ringdown,omg_QNM,[log(0.5),1]);
    disp(mdl)
    a = exp(mdl.Coefficients{1,1});
    theta = mdl.Coefficients{2,1};
    fprintf('a_lmn = %.5f\n',a)

    fit_omg = omg_QNM([log(a),theta],t_ringdown);

    figure
    tl = tiledlayout(3,1,'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl,1,[2,1]);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on 
    plot(t_ringdown-tLR,omg_ringdown,'DisplayName','Numerical','LineWidth',1.3,'LineStyle','-','Color',MyColors('b1'))
    plot(t_ringdown-tLR,fit_omg,'DisplayName','Fit','LineWidth',1,'LineStyle','--','Color',[.2 .2 .2]);
    ylabel(sprintf('$\\omega_{%d%d}$',l,m),'Interpreter','latex','FontSize',18);
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([tin-tLR,tend-tLR])
    %ylim([0.15,.85])
    lims = ylim;
    xlength = tend-tin;
    ylength = abs(lims(2)-lims(1));
    yline(omega0,'LineStyle','--','Color','k','HandleVisibility','off')
    text((xlength)/50+tin-tLR,omega0 + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omega0),'Interpreter','latex','FontSize',12)


    ax2 = nexttile(tl,3);
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    grid on
    plot(t_ringdown-tLR,abs(omg_ringdown-fit_omg)./omg_ringdown,'LineWidth',1,'Color','k');
    ylabel(sprintf('$\\Delta \\omega_{%d%d}/\\omega_{%d%d}$',l,m,l,m),'Interpreter','latex','FontSize',18);
    xlabel('$u-u_{LR}$','Interpreter','latex','FontSize',18);
    xlim([tin-tLR,tend-tLR])
    %yticks(ax2,[0,0.02,0.04,0.06])
    %ylim([0,0.07])


return