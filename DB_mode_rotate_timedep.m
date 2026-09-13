function p = DB_mode_rotate_timedep(s,dyn,l,M,p,varargin)
    % ==============================================================================================================================
    % NOTE
    % given an orientation of axes and some modes hlm computes in this reference frame, apply a rotation of euler angles alpha,beta
    % this script computes how the modes would look like in this new RF
    % https://arxiv.org/pdf/2005.05338 pag 8
    % ==============================================================================================================================

    iter = 0;
    direction = 'forward'; %can be either 'forward' (inertial to co-precessing) or 'backward' (co-precessing to inertial)
    
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'orbit'
                    i     = i + 1;
                    orbit = varargin{i};
                case 'r0'
                    i     = i + 1;
                    r0    = varargin{i};
                case 'obj'
                    i     = i + 1;
                    obj   = varargin{i};
                case 'iter'
                    i     = i + 1;
                    iter  = varargin{i};
                case 'direction'
                    i     = i + 1;
                    direction  = varargin{i}; 
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    dyn_time = dyn.t;
    %dt = 4565.2;
    %dyn_time = dyn.t - dt;
    %idx0 = find(dyn_time>0,1);
    %dyn_time = dyn_time(idx0:end);

    if iter==0
        L = dyn.L;
    elseif iter==1
        L = dyn.VOmg_schw;
    elseif iter==2
        L = dyn.VOmg_spin_only;
    elseif iter==3
        L = dyn.VOmg;
    end

    if dyn.geodesics==0
        tLR_dyn = tLR_splined(dyn);
        tLR_wf = tLR_splined(s.dyn);
        if tLR_dyn~=tLR_wf
            dt = tLR_dyn - tLR_wf;
            dyn_time = dyn.t - dt;
        end
    end
            
    %[enn,lamb,ell] = DB_coprecessing_basis(dyn); %attempt to use different co-precessing frame

    [~,dyn_phi,dyn_theta] = DB_coords_cart2spherical(L(:,1),L(:,2),L(:,3),0,0,0);
    
    %t  = s.ell(2).emm(3).t;
    t  = s.ell(l).emm(M+1).t; 
    if dyn.geodesics == 0
        j = find(t>dyn_time(end),1);
        %j = find(t>tLR_splined(dyn),1);

        j = j-1;
        t_inspiral = t(1:j);
    else
        t_inspiral = t;
    end

    interp_theta = spline(dyn_time,dyn_theta,t_inspiral);
    interp_phi = spline(dyn_time,dyn_phi,t_inspiral);

    beta  =  - interp_theta; %to match the two frames I have to rotate clockwise of theta
    gamma =  interp_phi - pi; %L starts at phi=pi, but the gamma angle of the Euler rotation is zero at t=0 (because I rotate of -beta)

    dBeta = DB_D1(beta,t_inspiral,4);
    integrand = gamma.*dBeta.*sin(beta);
    correction = cumtrapz(t_inspiral,integrand);

    alpha = - (gamma.*cos(beta) + correction);

    if strcmp(direction,'backward')
       
        temp = gamma;
        gamma = -alpha;
        beta = -beta;
        alpha = -temp;
    end


    % =========================================================================================================
    %DEBUG
    %{
    %figure
    %hold on
    %plot(dyn_time,dyn_theta,'DisplayName','\theta_L')
    %plot(dyn_time,dyn_phi,'DisplayName','\phi_L')
    %legend
    

    x = [1;0;0];
    y = [0;1;0];
    z = [0;0;1];

    %A = spline(t,alpha,dyn.t);
    %B = spline(t,beta,dyn.t);
    %C = spline(t,gamma,dyn.t);    
    A = alpha;
    B = beta;
    C = gamma;

    a = A(end);
    b = B(end);
    c = C(end);

    Ra = [[cos(a) -sin(a) 0]; [sin(a) cos(a) 0]; [0 0 1]];
    Rb = [[cos(b) 0 sin(b)]; [0 1 0]; [-sin(b) 0 cos(b)];];
    Rc = [[cos(c) -sin(c) 0]; [sin(c) cos(c) 0]; [0 0 1]];

    x1 = Rc*Rb*Ra*x;
    y1 = Rc*Rb*Ra*y;
    z1 = Rc*Rb*Ra*z;

    %disp(x1')
    %disp(enn(pos,:))
    %disp(y1')
    %disp(lamb(pos,:))
    disp(z1')

    Lx = spline(dyn_time,L(:,1),t_inspiral);
    Ly = spline(dyn_time,L(:,2),t_inspiral);
    Lz = spline(dyn_time,L(:,3),t_inspiral);
    disp([Lx(end) Ly(end) Lz(end)]./norm([Lx(end) Ly(end) Lz(end)]))
    
    %}
    % =========================================================================================================


    if dyn.geodesics == 0
        zeros_ringdown = t(j+1:end)*0;

        beta_ringdown = zeros_ringdown + beta(end);
        alpha_ringdown = zeros_ringdown + alpha(end);
        gamma_ringdown = zeros_ringdown + gamma(end);

        beta = [beta;beta_ringdown];
        alpha = [alpha;alpha_ringdown];
        gamma = [gamma;gamma_ringdown];
    end

    %plot_angles(t,alpha,beta,gamma);

    %idx0 = length(s.ell(4).emminus(5).t); %DEBUG
    %t = s.ell(4).emminus(5).t;
    %idx0 = length(t);
    h = 0;
    DWigner = t*0;
    for m=-l:l
        fprintf('m = %d ',m)
        if m<0
            hlm    = s.ell(l).emminus(-m+1).hlm;
        else
            hlm    = s.ell(l).emm(m+1).hlm;
        end

        for i=1:length(t)
            DWigner(i,:) = conj(wigner_D_function(l,m,M,-gamma(i),-beta(i),-alpha(i)));
            %DWigner(i,:) = conj(wigner_D_function(l,m,M,alpha(i),beta(i),gamma(i)));
            %DWigner(i,:) = conj(wigner_D_function(l,m,M,alpha(i),beta(i),gamma(i)));
            %DWigner(i,:) = (wigner_D_function(l,M,m,-alpha(i),-beta(i),-gamma(i)));
        end

        %h = h + hlm(1:idx0).*DWigner(1:idx0); %DEBUG
        h = h + hlm.*DWigner;
    end
    fprintf('\n')

    if M<0
        p.ell(l).emminus(-M+1).hlm = h;
        p.ell(l).emminus(-M+1).t = t;
    else
        p.ell(l).emm(M+1).hlm = h;
        p.ell(l).emm(M+1).t = t;
    end
        

    %wf.l = l;
    %wf.m = M;
    p.dyn = dyn;
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

function plot_angles(t,alpha,beta,gamma)
    figure
    hold on
    plot(t,alpha,'DisplayName','\alpha')
    plot(t,beta,'DisplayName','\beta')
    plot(t,gamma,'DisplayName','\gamma')
    legend
    xlabel('t')
    ylabel('rad')
    title('Euler Angles')
    hold off
return