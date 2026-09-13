function a = DB_fitKerrQNM(s,l,m,varargin)
    % =============================================================================
    % Used for fits in PaperI and various tests
    % Poorly written: uses kerr_FitQNMs from KerrOrbitSolver
    % =============================================================================
    
    j = s.dyn.chi1(3);
    model = 1;
    lbl = KerrLabel(s.dyn);
    if contains(lbl,'45')
        extra = 0;
    else
        extra = 5;
    end

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'model1'
                    model = 1;
                case 'model2'
                    model = 2;
                case 'model3'
                    model = 3;
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    if m<0
        h = s.ell(l).emminus(-m+1).hlm;
        t = s.ell(l).emminus(-m+1).t;
    else
        h = s.ell(l).emm(m+1).hlm;
        t = s.ell(l).emm(m+1).t;
    end

    omg = freq(h,t);
    tLR = tLR_splined(s.dyn);

    tin = tLR+30+extra;
    tend = tLR+75+extra;

    idx0 = find(t>tin,1);
    idx1 = find(t>tend,1);

    omg_ringdown = omg(idx0:idx1);
    t_ringdown = t(idx0:idx1);
    t_old = t;
    t = t_ringdown;
   
    %if model==1
    % ================================================================================================================================================================
    % working: fits interfeterence between positive and negative freq at given ml
        omega1 = sign(m)*imag(kerr_FitKerrQNMs(l,sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
        omega2 = sign(m)*imag(kerr_FitKerrQNMs(l,-sign(j)*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
        tau1 = real(kerr_FitKerrQNMs(l,sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
        tau2 = real(kerr_FitKerrQNMs(l,-sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));

        omg_QNM = @(C,t_ringdown)(-((exp(C(1)).^2.*exp(1).^(2.*t.*tau1)+exp(1).^(2.*t.*tau2)+2.*exp(C(1)).*exp(1) ...
                                    .^(t.*(tau1+tau2)).*cos((omega1+omega2).*t+exp(C(2)))).^(-1).*((-1).* ...
                                    exp(1).^(2.*t.*tau2).*omega1+exp(C(1)).^2.*exp(1).^(2.*t.*tau1).*omega2+ ...
                                    (-1).*exp(C(1)).*exp(1).^(t.*(tau1+tau2)).*(omega1+(-1).*omega2).*cos(( ...
                                    omega1+omega2).*t+exp(C(2)))+exp(C(1)).*exp(1).^(t.*(tau1+tau2)).*(tau1+(-1) ...
                                    .*tau2).*sin((omega1+omega2).*t+exp(C(2))))));

        mdl = fitnlm(t_ringdown,omg_ringdown,omg_QNM,[log(2),pi]);
        disp(mdl)
        a = exp(mdl.Coefficients{1,1});
        theta = mdl.Coefficients{2,1};
        fprintf('a_lmn = %.5f\n',a)
        
        fit_omg0 = omg_QNM([log(a),theta],t_ringdown);
    %else
    % ================================================================================================================================================================
    % try to include (l+1,m)

    sigma220p = kerr_FitKerrQNMs(l,sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma220m = kerr_FitKerrQNMs(l,-sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma320p = kerr_FitKerrQNMs(l+1,sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma320m = kerr_FitKerrQNMs(l+1,-sign(j).*abs(m),abs(j),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');

    omega220p = sign(m)*imag(sigma220p);
    tau220p = real(sigma220p);
    omega220m = -sign(m)*imag(sigma220m);
    tau220m = real(sigma220m);
    omega320p = sign(m)*imag(sigma320p);
    tau320p = real(sigma320p);
    omega320m = -sign(m)*imag(sigma320m);
    tau320m = real(sigma320m);

    omg_QNM = @(C,t)(exp(1).^((-1).*t.*(2.*tau220m+(-4).*tau220p+tau320m+tau320p)).*( ...
                        1  +exp(1).^(2.*t.*((-1).*tau220m+tau220p)).*C(1).^2+exp(1).^(2.*t.* ...
                        (  tau220p+(-1).*tau320p)).*C(2).^2+exp(1).^(2.*t.*(tau220p+(-1).* ...
                        tau320m  )).*C(3).^2+2.*exp(1).^(t.*((-1).*tau220m+tau220p)).*C(1).* ...
                        cos  (((-1).*omega220m+omega220p).*t+C(4))+2.*exp(1).^(t.*(tau220p+( ...
                        -1).*tau320p)).*C(2).*cos(omega220p.*t+(-1).*omega320p.*t+C(5))+ ...
                        2.  *exp(1).^((-1).*t.*(tau220m+(-2).*tau220p+tau320p)).*C(1).*C(2) ...
                        .*cos(omega220m.*t+(-1).*omega320p.*t+(-1).*C(4)+C(5))+2.*exp(1) ...
                        .^(t.*(2.*tau220p+(-1).*tau320m+(-1).*tau320p)).*C(2).*C(3).*cos( ...
                        omega320m  .*t+(-1).*omega320p.*t+C(5)+(-1).*C(6))+2.*exp(1).^(t.*( ...
                        tau220p  +(-1).*tau320m)).*C(3).*cos(omega220p.*t+(-1).*omega320m.* ...
                        t  +C(6))+2.*exp(1).^((-1).*t.*(tau220m+(-2).*tau220p+tau320m)).*C( ...
                        1  ).*C(3).*cos(omega220m.*t+(-1).*omega320m.*t+(-1).*C(4)+C(6))).^( ...
                        -1).*(exp(1).^(t.*(2.*tau220m+(-4).*tau220p+tau320m+tau320p)).* ...
                        omega220p  +exp(1).^(t.*((-2).*tau220p+tau320m+tau320p)).* ...
                        omega220m  .*C(1).^2+exp(1).^(t.*(2.*tau220m+(-2).*tau220p+tau320m+( ...
                        -1).*tau320p)).*omega320p.*C(2).^2+exp(1).^(t.*(2.*tau220m+(-2).* ...
                        tau220p  +(-1).*tau320m+tau320p)).*omega320m.*C(3).^2+exp(1).^(t.*( ...
                        tau220m  +(-3).*tau220p+tau320m+tau320p)).*(omega220m+omega220p).*C( ...
                        1  ).*cos(((-1).*omega220m+omega220p).*t+C(4))+exp(1).^(t.*(2.* ...
                        tau220m  +(-3).*tau220p+tau320m)).*(omega220p+omega320p).*C(2).*cos( ...
                        omega220p  .*t+(-1).*omega320p.*t+C(5))+exp(1).^(t.*(tau220m+(-2).* ...
                        tau220p  +tau320m)).*omega220m.*C(1).*C(2).*cos(omega220m.*t+(-1).* ...
                        omega320p  .*t+(-1).*C(4)+C(5))+exp(1).^(t.*(tau220m+(-2).*tau220p+ ...
                        tau320m  )).*omega320p.*C(1).*C(2).*cos(omega220m.*t+(-1).* ...
                        omega320p  .*t+(-1).*C(4)+C(5))+exp(1).^(2.*t.*(tau220m+(-1).* ...
                        tau220p  )).*omega320m.*C(2).*C(3).*cos(omega320m.*t+(-1).* ...
                        omega320p  .*t+C(5)+(-1).*C(6))+exp(1).^(2.*t.*(tau220m+(-1).* ...
                        tau220p  )).*omega320p.*C(2).*C(3).*cos(omega320m.*t+(-1).* ...
                        omega320p  .*t+C(5)+(-1).*C(6))+exp(1).^(t.*(2.*tau220m+(-3).* ...
                        tau220p  +tau320p)).*omega220p.*C(3).*cos(omega220p.*t+(-1).* ...
                        omega320m  .*t+C(6))+exp(1).^(t.*(2.*tau220m+(-3).*tau220p+tau320p)) ...
                        .*omega320m.*C(3).*cos(omega220p.*t+(-1).*omega320m.*t+C(6))+exp( ...
                        1  ).^(t.*(tau220m+(-2).*tau220p+tau320p)).*omega220m.*C(1).*C(3).* ...
                        cos  (omega220m.*t+(-1).*omega320m.*t+(-1).*C(4)+C(6))+exp(1).^(t.*( ...
                        tau220m  +(-2).*tau220p+tau320p)).*omega320m.*C(1).*C(3).*cos( ...
                        omega220m  .*t+(-1).*omega320m.*t+(-1).*C(4)+C(6))+exp(1).^(t.*( ...
                        tau220m  +(-3).*tau220p+tau320m+tau320p)).*tau220m.*C(1).*sin(((-1) ...
                        .*omega220m+omega220p).*t+C(4))+(-1).*exp(1).^(t.*(tau220m+(-3).* ...
                        tau220p  +tau320m+tau320p)).*tau220p.*C(1).*sin(((-1).*omega220m+ ...
                        omega220p  ).*t+C(4))+(-1).*exp(1).^(t.*(2.*tau220m+(-3).*tau220p+ ...
                        tau320m  )).*tau220p.*C(2).*sin(omega220p.*t+(-1).*omega320p.*t+C(5) ...
                        )+exp(1).^(t.*(2.*tau220m+(-3).*tau220p+tau320m)).*tau320p.*C(2).* ...
                        sin  (omega220p.*t+(-1).*omega320p.*t+C(5))+(-1).*exp(1).^(t.*( ...
                        tau220m  +(-2).*tau220p+tau320m)).*tau220m.*C(1).*C(2).*sin( ...
                        omega220m  .*t+(-1).*omega320p.*t+(-1).*C(4)+C(5))+exp(1).^(t.*( ...
                        tau220m  +(-2).*tau220p+tau320m)).*tau320p.*C(1).*C(2).*sin( ...
                        omega220m  .*t+(-1).*omega320p.*t+(-1).*C(4)+C(5))+(-1).*exp(1).^( ...
                        2.  *t.*(tau220m+(-1).*tau220p)).*tau320m.*C(2).*C(3).*sin( ...
                        omega320m  .*t+(-1).*omega320p.*t+C(5)+(-1).*C(6))+exp(1).^(2.*t.*( ...
                        tau220m  +(-1).*tau220p)).*tau320p.*C(2).*C(3).*sin(omega320m.*t+( ...
                        -1).*omega320p.*t+C(5)+(-1).*C(6))+(-1).*exp(1).^(t.*(2.*tau220m+( ...
                        -3).*tau220p+tau320p)).*tau220p.*C(3).*sin(omega220p.*t+(-1).* ...
                        omega320m  .*t+C(6))+exp(1).^(t.*(2.*tau220m+(-3).*tau220p+tau320p)) ...
                        .*tau320m.*C(3).*sin(omega220p.*t+(-1).*omega320m.*t+C(6))+(-1).* ...
                        exp  (1).^(t.*(tau220m+(-2).*tau220p+tau320p)).*tau220m.*C(1).*C(3) ...
                        .*sin(omega220m.*t+(-1).*omega320m.*t+(-1).*C(4)+C(6))+exp(1).^( ...
                        t  .*(tau220m+(-2).*tau220p+tau320p)).*tau320m.*C(1).*C(3).*sin( ...
                        omega220m  .*t+(-1).*omega320m.*t+(-1).*C(4)+C(6))));

        omg_QNM2 = @(C2,t)(omg_QNM([C2(1),C2(2),0,C2(3),C2(4),0],t));

        if model==2
            mdl = fitnlm(t_ringdown,omg_ringdown,omg_QNM2,[log(.1),log(.5),pi,pi/2]);
            disp(mdl)
            A1 = exp(mdl.Coefficients{1,1});
            A2 = exp(mdl.Coefficients{2,1});
            theta1 = mdl.Coefficients{3,1};
            theta2 = mdl.Coefficients{4,1};
            fprintf('A1 = %.5f\n A2 = %.5f\n',A1,A2)
            fit_omg = omg_QNM2([log(A1),log(A2),theta1,theta2],t);
        elseif model==3
            mdl = fitnlm(t_ringdown,omg_ringdown,omg_QNM,[log(.7),log(.1),log(10),pi,2*pi,pi/3]);
            disp(mdl)
            A1 = exp(mdl.Coefficients{1,1});
            A2 = exp(mdl.Coefficients{2,1});
            A3 = exp(mdl.Coefficients{3,1});
            theta1 = mdl.Coefficients{4,1};
            theta2 = mdl.Coefficients{5,1};
            theta3 = mdl.Coefficients{6,1};
            fprintf('A1 = %.5f\n A2 = %.5f\n A3 = \n%.5f\n',A1,A2,A3)

            fit_omg = omg_QNM([log(A1),log(A2),log(A3),theta1,theta2,theta3],t);
        end
        omega1 = omega220p;
    %end
    % ================================================================================================================================================================
    figure
    tl = tiledlayout(3,1,'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl,1,[2,1]);
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on 
    plot(t_old-tLR,omg,'DisplayName',sprintf('Numerical'),'LineWidth',1.3,'LineStyle','-','Color',MyColors('b1'))
    plot(t_ringdown-tLR,fit_omg0,'DisplayName','Fit  $\sigma_{220} + \sigma_{2-20}^*$','LineWidth',1.5,'LineStyle','--','Color',[.5 .5 .5]);
    plot(t_ringdown-tLR,fit_omg,'DisplayName','Fit  $\sigma_{220} + \sigma_{2-20}^* + \sigma_{320}$','LineWidth',1.5,'LineStyle','--','Color','r');
    ylabel(sprintf('$\\omega_{%d%d}$',l,m),'Interpreter','latex','FontSize',18);
    legend('Location','southeast','Interpreter','latex','FontSize',12   ,'BackgroundAlpha',0.7 + extra*0.04)
    xlim([tin-tLR,tend-tLR])
    %ylim([0.15,.85])
    lims = ylim;
    xlength = tend-tin;
    ylength = abs(lims(2)-lims(1));
    yline(omega1,'LineStyle','-.','Color',[.7 .7 .7],'HandleVisibility','off')
    %text((xlength)/50+tin-tLR,omega1 + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omega1),'Interpreter','latex','FontSize',12)
    text(.43*xlength + tin - tLR, lims(2) - ylength/10,lbl,'Interpreter','latex','FontSize',12)


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