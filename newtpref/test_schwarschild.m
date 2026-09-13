function test_schwarzschild(m)
    % ==========================================================================================================================================================================================
    % Some checks on Schwarzschild: shows that rotating the equatorial Newtonian multipoles one retrieves the correct Newtonian tilted multipole. 
    % Also shows that by attaching to each Ulm/Vlm the correct corrections the tilted waveform is more accurate than the one computed using the tilted dynamucs (which in a sense is obvious)
    % ==========================================================================================================================================================================================
    
    l=2;
    load("~/waveforms/S/fluxes/th0_90/wf.mat");
    seq = s;
    load("~/waveforms/S/fluxes/th0_60/wf.mat");
    stilt = s;

    %srot = struct;
    %srot = DB_mode_rotate_timedep(seq,stilt.dyn,l,m,srot,'direction','backward');

    beta = pi/6;
    d22 = wigner_D_function(l,2,m,0,-beta,0);
    d12 = wigner_D_function(l,1,m,0,-beta,0);
    dm12 = wigner_D_function(l,-1,m,0,-beta,0);
    dm22 = wigner_D_function(l,-2,m,0,-beta,0);

    flags.Newt = 0;
    flags.mass = 2;
    flags.curr = 1;
    [U22,V21] = Multipoles(seq.dyn,flags); %Compute multipoles from equatorial dynamics (add flag 'Newt' for Newtonian approx)

    flags.Newt = 1;
    flags.mass = m;
    flags.curr = m;
    [U22t,V22t] = Multipoles(stilt.dyn,flags); %Compute multipoles from tilted dynamics (add flag 'Newt' for Newtonian approx)

    %h22_tilt = U22.*(d22+dm22) - 1i.*V21.*(d12 - dm12);
    h22_tilt = U22.*d22 + conj(U22).*dm22 - 1i.*(V21.*d12 - dm12.*conj(V21));
    h22_tilt_orig = U22t - 1i.*V22t;


    % ========================================================
    % FOR THE PLOTS
    % ========================================================
    tLR = tLR_splined(stilt.dyn);
    T = stilt.ell(l).emm(m+1).t;
    Hlm = stilt.ell(l).emm(m+1).hlm;

    T_an = seq.dyn.t;
    Hlm_an = h22_tilt;

    %T_rot = srot.ell(l).emm(m+1).t;
    %Hlm_rot = srot.ell(l).emm(m+1).hlm;

    Hlm_an_s = spline(T_an,Hlm_an,T);
    %Hlm_rot_s = spline(T_rot,Hlm_rot,T);

    DeltaA_an = (abs(Hlm)-abs(Hlm_an_s))./abs(Hlm);
    %DeltaA_rot = (abs(Hlm)-abs(Hlm_rot_s))./abs(Hlm);

    DeltaPhi_an = (-unwrap(angle(Hlm_an_s))+unwrap(angle(Hlm)));
    %DeltaPhi_rot = (-unwrap(angle(Hlm_rot_s))+unwrap(angle(Hlm)));

    Hlm_an_tilt_s = spline(stilt.dyn.t,h22_tilt_orig,T);

    M = mean(DeltaPhi_an(200:end-8000));
    while abs(M)>pi/2
        if M>0
            DeltaPhi_an = DeltaPhi_an - pi;
        else
            DeltaPhi_an = DeltaPhi_an + pi;
        end
        M = mean(DeltaPhi_an(200:end-8000));
    end
    %N = mean(DeltaPhi_rot(200:end-200));
    %while abs(N)>pi/2
    %    if N>0
    %        DeltaPhi_rot = DeltaPhi_rot - pi;
    %    else
    %        DeltaPhi_rot = DeltaPhi_rot + pi;
    %    end
    %    N = mean(DeltaPhi_rot(200:end-200));
    %end
    T = T-tLR;

    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    f = figure;
    f.Position(3:4) = f.Position(3:4)*1.2;
    tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
    t_in = -620;
    t_end = -0;
    label = KerrLabel(stilt.dyn);
    wf_label = sprintf('h_{%d%d}',l,m);

    nexttile(tl,[2,1])
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
    xlim([t_in,t_end])
    hold on
    plot(T,real(Hlm),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s]$',wf_label))
    %plot(T,real(Hlm_rot_s),'LineWidth',my_linewidth,'DisplayName',sprintf('Numerical Rotated'),'LineStyle',':','LineWidth',1.25)
    plot(T,real(Hlm_an_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s^{\\rm EOB}]$ (equatorial rotated)',wf_label),'LineStyle','--') 
    %plot(T,real(Hlm_an_tilt_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s^{\\rm EOB}]$ (equatorial rotated)',wf_label),'LineStyle','--') 

    plot(T,abs(Hlm),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s|$',wf_label),'Color',[1 0 0 .2],'HandleVisibility','off')
    plot(T,abs(Hlm_an_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s^{\\rm EOB}|$ - Equatorial Rotated',wf_label),'LineStyle','--','Color',[1 0 0 .2],'HandleVisibility','off')
    %plot(T,abs(Hlm_an_tilt_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s^{\\rm EOB}|$ - Equatorial Rotated',wf_label),'LineStyle','--','Color',[1 0 0 .2],'HandleVisibility','off')

    ys = ylim;
    xs = xlim;
    text(xs(1)+abs(xs(2)-xs(1))./6,ys(2) - abs(ys(2)-ys(1))./8,label,'Interpreter','latex','FontSize',16)
    
    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlim([t_in,t_end])
    hold on
    plot(T,abs(DeltaPhi_an),'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$');%_{\rm analytical}$')
    %plot(T,abs(DeltaPhi_rot),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi_{\rm numerical}$')
    plot(T,abs(DeltaA_an),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$[\Delta A/A]$');%_{\rm analytical}$')
    %plot(T,abs(DeltaA_rot),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$[\Delta A/A]_{\rm numerical}$')

    yscale log
    yticks([1e-4,1e-3,1e-2,1e-1,1])
    grid off
    grid on
    ylim([-.02,.02])

    %figure
    %hold on 
    %plot(stilt.ell(l).emm(m+1).t,real(stilt.ell(l).emm(m+1).hlm),'DisplayName','Numerical Tilted')
    %plot(seq.dyn.t,real(h22_tilt),'LineStyle','--','DisplayName','Analytical rotated')
    %plot(srot.ell(l).emm(m+1).t,real(srot.ell(l).emm(m+1).hlm),'DisplayName','Numerical Rotated')
    %xlim([200,1300])
    %legend

    f = figure;
    f.Position(3:4) = f.Position(3:4)*1.2;
    tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');

    nexttile(tl,[2,1])
    hold on
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
    xlim([t_in,t_end])

    plot(T,real(Hlm_an_tilt_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s^{\\rm Newt}]$ (tilted)',wf_label),'LineStyle','-') 
    plot(T,real(Hlm_an_s),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re[%s^{\\rm Newt}]$ (equatorial rotated)',wf_label),'LineStyle','--') 
    ys = ylim;
    xs = xlim;
    text(xs(1)+abs(xs(2)-xs(1))./6,ys(1) + abs(ys(2)-ys(1))./8,label,'Interpreter','latex','FontSize',16)

    nexttile
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
    xlim([t_in,t_end])
    hold on
    plot(T,-unwrap(angle(Hlm_an_tilt_s))+unwrap(angle(Hlm_an_s)),'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$');
    plot(T,(abs(Hlm_an_tilt_s)-abs(Hlm_an_s))./abs(Hlm_an_tilt_s),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$[\Delta A/A]$');%_{\rm analytical}$')
    grid off
    yscale log
    yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1])
    grid on

return

function out = wigner_D_function(l,m1,m,alpha,beta,gamma)  %as written, this computes D^l_{m1,m} (order of indices is important)
    cth  = cos(beta*0.5);
    sth  = sin(beta*0.5);
    norm = sqrt( (factorial(l+m1) * factorial(l-m1) * factorial(l+m) * factorial(l-m)) );
    ki   = max(0,m-m1);
    kf   = min(l+m,l-m1);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+m-k) * factorial(l-m1-k) * factorial(k-m+m1) );
        dWig = dWig + div*( (-1).^(k) * cth.^(2*l+m-m1-2*k) * sth.^(2*k+m1-m) );
    end
    dWig_normalized = dWig*norm;
    out = exp(1i.*alpha.*m1).*dWig_normalized.*exp(1i.*gamma.*m);
return

function [U,V] = Multipoles(DB,flags)
    l=2;

    a = DB.chi1(3);
    t = DB.t;
    r = DB.r;

    x = DB.x.*(1+a.*r.^(-3/2)).^(2/3);
    y = DB.y.*(1+a.*r.^(-3/2)).^(2/3);
    z = DB.z.*(1+a.*r.^(-3/2)).^(2/3);

    rOmg = r.*(1+DB.chi1(3).*DB.r.^(-3/2)).^(2/3);
    %[x,y,z] = DB_coords_spherical2cart(rOmg,DB.phi,DB.th,0,0,0); %Compute prefactor with rOmg instead of r. Is more accurate

    r0    = 1.213061319425267e+00;   % 2/sqrt(e);
    Omega = DB.Omg;
    Heff = DB.Heff;
    jhat = DB_jhat(DB);

    %jhat = DB.pphi./(rOmg.^2.*Omega);
    Tail22 = DB_Tail(l,2.*Omega,2.*Omega,r0);
    Tail21 = DB_Tail(l,Omega,Omega,r0);

    %a = DB.chi1(3)./sin(DB.th0);
    rho = DB_rholm(Omega,a,r);
    delta = DB_deltalm(Omega,a);

    xdot = DB_D1(x,t,4);
    x2dot = DB_D1(xdot,t,4);
    x3dot = DB_D1(x2dot,t,4);
    x4dot = DB_D1(x3dot,t,4);
    x5dot = DB_D1(x4dot,t,4);

    ydot = DB_D1(y,t,4);
    y2dot = DB_D1(ydot,t,4);
    y3dot = DB_D1(y2dot,t,4);
    y4dot = DB_D1(y3dot,t,4);
    y5dot = DB_D1(y4dot,t,4);

    zdot = DB_D1(z,t,4);
    z2dot = DB_D1(zdot,t,4);
    z3dot = DB_D1(z2dot,t,4);
    z4dot = DB_D1(z3dot,t,4);
    z5dot = DB_D1(z4dot,t,4);

    rho22 = rho.l2m2;
    delta22 = delta.l2m2;

    xPN = (rOmg.*Omega).^2;
    %rho21 = (1-0.72.*xPN);
    rho21 = rho.l2m1;
    delta21 = delta.l2m1;

    MassQuadrupole_22 = 2.*((2/15).*pi).^(1/2).*(xdot.^2+(sqrt(-1)*(-1)).*x2dot.*y+x.*( ...
            x2dot+(sqrt(-1)*(-1)).*y2dot)+(-1).*y.*y2dot+(sqrt(-1)*(-2)).* ...
            xdot.*ydot+(-1).*ydot.^2);

    MassQuadrupole_21 = (-2).*((2/15).*pi).^(1/2).*(x2dot.*z+(sqrt(-1)*(-1)).*(y2dot.*z+ ...
            sqrt(-1).*x.*z2dot+y.*z2dot+(sqrt(-1)*2).*xdot.*zdot+2.*ydot.* ...
            zdot));

    CurrQuadrupole_21 = ((2/15).*pi).^(1/2).*((sqrt(-1)*(-1)).*x3dot.*y.^2+(-1).*x.^2.* ...
            y3dot+sqrt(-1).*x3dot.*z.^2+y3dot.*z.^2+(-1).*y.*z.*z3dot+3.* ...
            y2dot.*z.*zdot+(-3).*y.*z2dot.*zdot+3.*x2dot.*(xdot.*y+(sqrt(-1)*( ...
            -1)).*y.*ydot+sqrt(-1).*z.*zdot)+x.*(x3dot.*y+(-3).*xdot.*y2dot+ ...
            sqrt(-1).*(y.*y3dot+3.*y2dot.*ydot+(-1).*z.*z3dot+(-3).*z2dot.* ...
            zdot)));

    CurrQuadrupole_22 = (-1).*((2/15).*pi).^(1/2).*(x3dot.*y.*z+3.*xdot.*y2dot.*z+(sqrt( ...
            -1)*(-1)).*y.*y3dot.*z+(sqrt(-1)*(-3)).*y2dot.*ydot.*z+3.*x2dot.*( ...
            sqrt(-1).*xdot+ydot).*z+(-3).*xdot.*y.*z2dot+(sqrt(-1)*3).*y.* ...
            ydot.*z2dot+(sqrt(-1)*(-1)).*x.^2.*z3dot+sqrt(-1).*y.^2.*z3dot+x.* ...
            (sqrt(-1).*x3dot.*z+y3dot.*z+(sqrt(-1)*(-3)).*xdot.*z2dot+(-3).* ...
            ydot.*z2dot+(-2).*y.*z3dot));

    if flags.Newt
        if flags.mass==2
            U = 2*sqrt(3).*MassQuadrupole_22.*(1./sqrt(2));
        elseif flags.mass==1
            U = 2*sqrt(3).*MassQuadrupole_21.*(1./sqrt(2));
        end
        if flags.curr==1
            V = -8./sqrt(3).*CurrQuadrupole_21.*(1./sqrt(2));
        elseif flags.curr==2
            fprintf('V22\n')
            V = -8./sqrt(3).*CurrQuadrupole_22.*(1./sqrt(2));
        else
            error('flag.m non valid')
        end
    else
        U = 2*sqrt(3).*MassQuadrupole_22.*(1./sqrt(2)).*Heff.*Tail22.*(rho22).^l.*exp(1i.*delta22);
        if flags.curr==1
            V = -8./sqrt(3).*CurrQuadrupole_21.*(1./sqrt(2)).*jhat.*Tail21.*(rho21).^l.*exp(1i.*delta21);
        elseif flags.curr==2
            V = -8./sqrt(3).*CurrQuadrupole_22.*(1./sqrt(2)).*jhat.*Tail22.*(rho22).^l.*exp(1i.*delta22);
        else
            error('flag.m non valid')
        end
    end

return


