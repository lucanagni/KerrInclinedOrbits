function mechanical_exact_flux(s,which)
    %if abs(s.dyn.chi1(3)) < 0.9
    %    ellmax = 4;
    %else
    %    ellmax = 3;
    %end
    ellmax = 4;
    [dEdt,t] = DB_dEdt(s,'ellmax',ellmax);
    [dJdt,t1] = DB_dJdt(s,'ellmax',ellmax);

    dyn = s.dyn;
    T = dyn.t;
    n = length(T);

    dEdt_mech = zeros(n,1);
    F_mech = zeros(3,n);

    name = KerrLabel(s.dyn);

    for i=1:n
        X = [dyn.x(i); dyn.y(i); dyn.z(i)];
        P = [dyn.px(i); dyn.py(i); dyn.pz(i)];
        dHdp = dyn.dHdp(i,:)';
        if which==1
            [f,dEdt_mech(i)] = DB_flux(X,P,dHdp,dyn.q,dyn.chi1);
        elseif which==2
            [f,dEdt_mech(i)] = DB_flux2(X,P,dHdp,dyn.q,dyn.chi1);
        else
            error('choose RR either 1 or 2')
        end
        %[~,~,~,fsr,fsphi,fsth] = DB_coords_cart2spherical(dyn.x(i),dyn.y(i),dyn.z(i),f(1),f(2),f(3));
        F_mech(:,i) = cross(X,f);
        %F_mech_norm(i) = sqrt(fsth.^2+(1./sin(dyn.th(i)).^2).*fsphi.^2)./DB_nuX1X2(dyn.q);
    end
    dEdt_mech = -dEdt_mech./DB_nuX1X2(dyn.q);
    F_mech_norm = sqrt(F_mech(1,:).^2 + F_mech(2,:).^2 + F_mech(3,:).^2);
    dJdt_norm = sqrt(dJdt(:,1).^2 + dJdt(:,2).^2 + dJdt(:,3).^2).*DB_nuX1X2(dyn.q);

    %[~,t_peak] = Apeak(s.ell(2).emm(3).hlm,s.ell(2).emm(3).t);
    t_peak = tLR_splined(dyn);
    t_LSSO = tLSSO_splined(dyn);

    %figure
    %set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    %xlabel('$u/M$','FontSize',14,'Interpreter','Latex');
%
    %hold on
    %plot(t,dEdt,'DisplayName','$\dot{E}$','LineWidth',1)
    %plot(T,dEdt_mech,'DisplayName','$\hat{\mathcal{F}}^\mathrm{mech}/\nu$','LineWidth',1,'LineStyle','-')
    %legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    %xlim([t_LSSO-35,t_peak + 100])
    %xline(t_LSSO,'LineStyle','-.','Color',[.7 .7 .7],'HandleVisibility','off')
    %xline(tLR_splined(dyn),'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')

    if abs(F_mech_norm(end))>10e-3
        y2 = max(dJdt_norm(round(length(dJdt_norm)./2):end))*1.5;
    end

    f = figure;
    f.Position(3:4) = f.Position(3:4)*1.1;
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,[2,1]);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    %xlabel('$u$','FontSize',18,'Interpreter','Latex');

    hold on
    plot(t1,dJdt_norm,'DisplayName','$\dot{J}/\nu$','LineWidth',1)
    plot(T,F_mech_norm,'DisplayName','$\hat{\mathcal{F}}^\mathrm{mech}$','LineWidth',1,'LineStyle','-')
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)
    xlim([300,t_peak + 100])
    xline(t_LSSO,'LineStyle','-.','Color',[.7 .7 .7],'HandleVisibility','off')
    xline(tLR_splined(dyn),'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')

    if abs(F_mech_norm(end))<=10e-3
        mylims = ylim;
        y2 = mylims(2).*1.1;
    end

    ylim([0,y2])

    mylims = ylim;
    text(t_peak-180,mylims(2)-(mylims(2)-mylims(1))/10,name,'Interpreter','latex','FontSize',16)

    axes('Position',[0.36,0.57,0.34,0.29]);   
    set(gca,'LineWidth',0.5,'FontName','Times','FontSize',7,'XMinorTick','on','YMinorTick','on','ZMinorTick','on','box','on')

    hold on
    plot(t1,dJdt_norm,'DisplayName','$\dot{J}/\nu$','LineWidth',1)
    plot(T,F_mech_norm,'DisplayName','$\hat{\mathcal{F}}^\mathrm{mech}$','LineWidth',1,'LineStyle','-')
    xlim([t_peak - 50,t_peak + 50])
    %xline(t_LSSO,'LineStyle','-.','Color',[.7 .7 .7],'HandleVisibility','off')
    xline(tLR_splined(dyn),'LineStyle','--','Color',[.7 .7 .7],'HandleVisibility','off')

    ylim([0,y2])

    %Compute relative difference |F^mech - Jdot|/Jdot
    reldiff = abs(spline(T,F_mech_norm,t1) - dJdt_norm)./dJdt_norm;

    ax2 = nexttile(tl);
    hold on
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',12,'FontName','Times','LineWidth',0.5,'box','on');
    xlabel('$u$','FontSize',18,'Interpreter','Latex');
    ylabel('$\Delta$','FontSize',18,'Interpreter','Latex')
    grid on
    plot(t1,reldiff,'Color','k')
    xline(t_LSSO,'LineStyle','-.','Color','k','HandleVisibility','off')
    xline(tLR_splined(dyn),'LineStyle','--','Color','k','HandleVisibility','off')
    xlim([300,t_peak + 100])

    if abs(s.dyn.th0-pi/2)>1e-5
        ylim([-0.02,0.1])
    else
        ylim([-0.01,0.04])
    end
    ylim([-0.1,.25])


return