function [a,t] = DB_asymmetry(s,varargin)
    % ========================================================================
    % Compute asymmetry of a given mode. Used in DB_many_asymmetries
    % ========================================================================
    
    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 18;

    label = KerrLabel(s.dyn);

    if abs(s.dyn.chi1(3))~=0.9
        ellmax = 4;
    else
        ellmax = 3;
    end
     
    num = 0;
    den = 0;

    L = 2*ellmax + 1;
    tend_vec = zeros(L,1);

    for m=-ellmax:ellmax
        if m<0
            tend = s.ell(ellmax).emminus(-m+1).t(end);
        else
            tend = s.ell(ellmax).emm(m+1).t(end);
        end

        tend_vec(m + ellmax + 1) = tend;
    end

    tstart = s.ell(2).emm(3).t(1);
    tend_ref = min(tend_vec);
    t = linspace(tstart,tend_ref,5e3)';


    s = spline_struct(s,t,ellmax);


    for l=2:ellmax
        for m = -l:l
            if m<0
                num = num + abs(s.ell(l).emminus(-m+1).hlm - (-1)^(l+m).*conj(s.ell(l).emm(-m+1).hlm)).^2;
                den = den + abs(s.ell(l).emminus(-m+1).hlm).^2;
            else
                num = num + abs(s.ell(l).emm(m+1).hlm - (-1)^(l+m).*conj(s.ell(l).emminus(m+1).hlm)).^2;
                den = den + abs(s.ell(l).emm(m+1).hlm).^2;
            end
        end
    end

    a = sqrt(num./(4.*den));
    %t = s.ell(2).emm(3).t;
    %r = a;

    if ~isempty(varargin)
        return
    end 
    
    %{
    if dyn.th0 ~= pi/2
        for l=2:ellmax
            for m=0:l
                sprintf('m = %d',m)
                p.ell(l).emm(m+1).hlm = DB_mode_rotate_timedep(s,dyn,l,m);
            end
            for n=1:l
                sprintf('m = %d',-n)
                p.ell(l).emminus(n+1).hlm = DB_mode_rotate_timedep(s,dyn,l,-n);
            end
            p.ell(l).emminus(1) = p.ell(l).emm(1);
        end

        num = 0;
        den = 0;

        for l=2:ellmax
            for m = 0:l
                num = num + abs(p.ell(l).emm(m+1).hlm - (-1)^(l+m).*conj(p.ell(l).emminus(m+1).hlm)).^2;
                den = den + abs(p.ell(l).emm(m+1).hlm).^2;
            end
            for n=1:l
                m=-n;
                num = num + abs(p.ell(l).emminus(n+1).hlm - (-1)^(l+m).*conj(p.ell(l).emm(n+1).hlm)).^2;
                den = den + abs(p.ell(l).emminus(n+1).hlm).^2;
            end
        end
    
        r = sqrt(num./(4.*den));
        %}

        figure
        %tiledlayout(2,1,'Padding','compact','TileSpacing','compact')
        
        %nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,a,'LineWidth',my_linewidth,'DisplayName',label)
        %plot(t,r,'DisplayName','rotated precessing')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
        xlabel('t','FontSize',labels_fontsize,'Interpreter','latex')
        ylabel('a','FontSize',labels_fontsize,'Interpreter','latex')
        xlim([200,t(end)])
        xline(tLR_splined(s.dyn),'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        %yscale log

        %{
        nexttile
        plot(t,abs(a-r)./a)
        grid on
        xlim([200,dyn.t(end)])
        xlabel('t')
        ylabel('\Deltaa')
        title('Relative Difference')

        nexttile(1)
        axes('Position',[.1 .7 .2 .2])
        box on
        plot(t,a)
        xlim([1500,1700])
        ylim([.135,.155])
        hold on
        plot(t,r)
    end
            %}

return

function [h,t] = rotation(w,t,dyn,l,M,varargin)
    mmin=0;
    dyn_time = dyn.t;

    L = cross([dyn.x,dyn.y,dyn.z],[dyn.px,dyn.py,dyn.pz]);
    [~,dyn_phi,dyn_theta] = DB_coords_cart2spherical(L(:,1),L(:,2),L(:,3),0,0,0);
    
    if dyn.geodesics == 0
        t_max = dyn.t(end);
        for i=1:length(t)
            if t(i)>t_max
                j=i-1;
                break
            end
        end
        t_inspiral = t(1:j);
    else
        t_inspiral = t;
    end

    interp_theta = spline(dyn_time,dyn_theta,t_inspiral);
    interp_phi = spline(dyn_time,dyn_phi,t_inspiral);

    beta  =  -interp_theta; %to match the two frames I have to rotate clockwise of theta
    gamma =  pi-interp_phi; %L starts at phi=pi, but the alpha angle of the Euler rotation is zero at t=0 (because I rotate of -beta)

    dBeta = DB_D1(beta,t_inspiral,4);
    integrand = gamma.*dBeta.*sin(beta);
    correction = cumtrapz(t_inspiral,integrand);

    %alpha = -gamma.*cos(beta);
    alpha = -gamma.*cos(beta) - correction;

    if dyn.geodesics == 0
        zeros_ringdown = t(j+1:end)*0;

        beta_ringdown = zeros_ringdown + beta(end);
        alpha_ringdown = zeros_ringdown + alpha(end);
        gamma_ringdown = zeros_ringdown + gamma(end);

        beta = [beta;beta_ringdown];
        alpha = [alpha;alpha_ringdown];
        gamma = [gamma;gamma_ringdown];
    end
    
    hre = 0;
    him = 0;
    dWigner = t*0;
    for m=mmin:l
        hlm    = w.ell(l).emm(m+1).hlm;
        philm  = unwrap(angle(hlm));
        Alm    = abs(hlm);
        
        cp       = cos(philm);
        sp       = sin(philm);

        for i=1:length(t)
            dWigner(i,:)  = wigner_d_function(l,m,M,-beta(i));
        end

        DWigner = exp(-1i*M*alpha).*dWigner.*exp(-1i*m*gamma);
        d_r = real(DWigner);
        d_i = imag(DWigner);
    
        hre    = hre + Alm.*(d_r.*cp - d_i.*sp);
        him    = him + Alm.*(d_r.*sp + d_i.*cp);
    end
    for n=1:l
        m=-n; 
        hlm    = w.ell(l).emminus(n+1).hlm;
        philm  = unwrap(angle(hlm));
        Alm    = abs(hlm);
        
        cp     = cos(philm);
        sp     = sin(philm);

        for i=1:length(t)
            dWigner(i,:)  = wigner_d_function(l,m,M,-beta(i));
        end

        DWigner = exp(-1i*M*alpha).*dWigner.*exp(-1i*m*gamma);
        d_r = real(DWigner);
        d_i = imag(DWigner);
        
        hre    = hre + Alm.*(d_r.*cp - d_i.*sp);
        him    = him + Alm.*(d_r.*sp + d_i.*cp);
    end


    h = hre + 1i.*him;
return

function out = wigner_d_function(l,m,s,beta)  %as written, this computes d^l_{m,s} (order of indices is important)
    cth  = cos(beta*0.5);
    sth  = sin(beta*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,s-m);
    kf   = min(l+s,l-m);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+s-k) * factorial(l-k-m) * factorial(k-s+m) );
        dWig = dWig + div*( (-1).^(k) * cth.^(2*l-m+s-2*k) * sth.^(2*k-s+m) );
    end
    out = dWig*norm;
return

function factor = fix_mode(type,l,m)
    if strcmp(type,'none')
        factor = 1;
    elseif strcmp(type,'eob')
        factor = (-1).^m;
    elseif strcmp(type, 'rwz')
        factor = (-1).^(l);
        if mod(l+m,2) % if odd
            factor = -factor*1i;
        end
    else
        error("type='%s' is not valid. Use 'rwz' or 'eob'", type)
    end
return 

function s = spline_struct(s,t,ellmax)
    for l=2:ellmax
        for m=-l:l
            if m<0
                s.ell(l).emminus(-m+1).hlm = spline(s.ell(l).emminus(-m+1).t,s.ell(l).emminus(-m+1).hlm,t);
                s.ell(l).emminus(-m+1).t = t;
            else
                s.ell(l).emm(m+1).hlm = spline(s.ell(l).emm(m+1).t,s.ell(l).emm(m+1).hlm,t);
                s.ell(l).emm(m+1).t = t;
            end
            s.ell(l).emminus(1).hlm = s.ell(l).emm(1).hlm;
            s.ell(l).emminus(1).t = s.ell(l).emm(1).t;
        end
    end

return