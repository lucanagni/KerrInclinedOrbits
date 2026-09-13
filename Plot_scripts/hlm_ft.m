function hlm_ft(s,l,m)
    a = s.dyn.chi1(3);
    n = abs(m);
    t = s.ell(l).emm(n+1).t;
    sigm = 1./(1+exp((-t+150)./30));
    if m>0
        FT = fft(s.ell(l).emm(m+1).hlm.*sigm);
    else
        FT = fft(s.ell(l).emminus(-m+1).hlm.*sigm);
    end

    sigma = kerr_FitKerrQNMs(l,2,a,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    omg_QNM = imag(sigma);

    dt = t(2) - t(1);
    Fs = round(1/dt);

    L = length(t);

    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(Fs/L*(-L/2:L/2-1),abs(fftshift(FT)),'DisplayName','$\omega_{\ell m}>0$')
    plot(-Fs/L*(-L/2:L/2-1),abs(fftshift(FT)),'DisplayName','$\omega_{\ell m}<0$')
    xlabel(sprintf('$|\\omega_{%d%d}|/2\\pi$',l,m),'FontSize',18,'Interpreter','Latex');
    ylabel(sprintf('FT$[h_{%d%d}]$',l,m),'FontSize',18,'Interpreter','Latex');

    xline(omg_QNM/2/pi,'LineStyle','--','Color',[.7 .7 .7],'DisplayName',sprintf('$\\omega^\\mathrm{QNM}$'))

    legend('Location','northeast','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.6)

    xscale log
    yscale log

    % ======================================================================================================
    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat')
    [hlm,t] = DB_hpc(s,0,0,'ellmax',3);
    FT = fft(hlm);
    a = 0;
    sigma = kerr_FitKerrQNMs(3,3,a,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    sigma2 = kerr_FitKerrQNMs(l+1,2,a,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
    %f_QNM = omg_QNM/2/pi;
    omega = Fs/L*(-L/2:L/2-1)*2*pi;
    OMEGA = -1i*sigma(2);

    OMEGA2 = 1i.*sigma2;
    F = ((omega - OMEGA)./(omega - conj(OMEGA)))';
    F2 = ((omega - OMEGA2)./(omega - conj(OMEGA2)))';
    %F = F.*F2;
    inv = ifft(F.*fftshift(FT)).*exp(1i.*(omega(end)-omega(1))./2.*t);
    final_freq = DB_D1(-unwrap(angle(inv)),t,4);

    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',10,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    plot(t,final_freq,'LineWidth',1)




return  