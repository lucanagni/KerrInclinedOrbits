function out = DB_spin_correction(s,p,l,m)
    % s = precessing waveform and dunamics 
    % p = rotated equatorial waveform with equivalent spin projection (i think)
    % Implements PN precessing effects from Kidder Boyle paper 
    %NOT WORKING    

    h_num = s.ell(l).emm(m+1).hlm;
    t_num = s.ell(l).emm(m+1).t;
    dyn = p.dyn;

    max_iter = length(dyn.t);

    vd = abs(dyn.Omg).^(1/3);

    R = [dyn.x,dyn.y,dyn.z];

    n = zeros(max_iter,3);
    l = zeros(max_iter,3);
    lambda = zeros(max_iter,3);
    
    dndt = zeros(max_iter,3);
    v = zeros(max_iter,1);

    for i=1:max_iter
        n(i,:) = R(i,:)./norm(R(i,:));
    end
    for i=1:3
        dndt(:,i) = DB_D1(n(:,i),dyn.t,4);
    end 
    for i=1:max_iter
        lambda(i,:) = dndt(i,:)./norm(dndt(i,:));
        l(i,:) = cross(n(i,:),lambda(i,:));
        Omg_orb = cross(n(i,:),dndt(i,:));
        v(i,:) = norm(Omg_orb).^(1/3);
    end
    %   v = mean(v);
    v = vd;
    v2 = v.^2;
    v3 = v.^3;
    v4 = v.^4;

    % =========================================================================================
    % Attempt to use the angles of L instead of the ones in the paper 
%{
    [~,dyn_phi,dyn_theta] = DB_coords_cart2spherical(dyn.L(:,1),dyn.L(:,2),dyn.L(:,3),0,0,0);
    normL = norm([dyn.L(:,1),dyn.L(:,2),dyn.L(:,3)]);
    interp_theta = spline(dyn.t,dyn_theta,t_num);
    interp_phi = spline(dyn.t,dyn_phi,t_num);

    beta  =  -interp_theta; %to match the two frames I have to rotate clockwise of theta
    gamma =  pi-interp_phi; %L starts at phi=pi, but the alpha angle of the Euler rotation is zero at t=0 (because I rotate of -beta)

    dBeta = DB_D1(beta,t_num,4);
    integrand = gamma.*dBeta.*sin(beta);
    correction = cumtrapz(t_num,integrand);

    alpha = -gamma.*cos(beta) - correction;

    for i=1:max_iter
        normL = norm([dyn.L(i,1),dyn.L(i,2),dyn.L(i,3)]);

        a = alpha(i);
        b = beta(i);
        c = gamma(i);

        x = [1;0;0];
        y = [0;1;0];
        z = [0;0;1];

        Ra = [[cos(a) -sin(a) 0]; [sin(a) cos(a) 0]; [0 0 1]];
        Rb = [[cos(b) 0 sin(b)]; [0 1 0]; [-sin(b) 0 cos(b)];];
        Rc = [[1 0 0]; [0 cos(c) -sin(c)]; [0 sin(c) cos(c)];];

        n(i,:) = (Rc*Rb*Ra*x)';
        lambda(i,:) = (Rc*Rb*Ra*y)';
        %l(i,:) = (Rc*Rb*Ra*z)';
        l(i,:) = dyn.L(i,:)./normL;
    end
%}
    % =========================================================================================

    M = 1;
    [nu,M1,M2] = DB_nuX1X2(dyn.q);
    delta = (M2 - M1)./M;
    %S = dyn.chi1;    my_zeros = zeros(max_iter,1);

    Slambda = zeros(max_iter,1);
    Sn = zeros(max_iter,1);
    Sl = zeros(max_iter,1);

    for i=1:max_iter
        Sn(i,:) = dot(dyn.chi1,n(i,:));
        Slambda(i,:) = dot(dyn.chi1,lambda(i,:));
        Sl(i,:) = dot(dyn.chi1,l(i,:));
    end

    Sigman = Sn.*M./M1;
    Sigmalambda = Slambda.*M./M1;
    Sigmal = Sl.*M./M1;


    if m==2
        h_hat = ... %-(2.*v3.*(3.*Sl + delta.*Sigmal))./3
        - v2.*(Sigmalambda + 1i.*Sigman)./2 ...
        + v4.*(182.*1i.*delta.*Sn + 19.*delta.*Slambda + 14.*1i.*(7 - 20.*nu).*Sigman + (5-43.*nu).*Sigmalambda)./84;

        %{
        h_hat = - v2.*(Sigmalambda + 1i.*Sigman)./2 ...
        + v4.*(182i.*delta.*Sn + 19.*delta.*Slambda + 14i.*(7 - 20.*nu).*Sigman + ...
        (5-43.*nu).*Sigmalambda)./84;
        %}
    elseif m==1
        h_hat = ...%1i.*v2.*Sigmal./2 - 1i.*v4.*(86.*delta.*Sl + (79 - 139.*nu).*Sigmal)./42 + ...
        v3.*(25.*Sn + 4i.*Slambda + 13.*delta.*Sigman + 4i.*delta.*Sigmalambda)./6; 
    elseif m==0
        h_hat = 1i.*(v2.*Sigman./sqrt(6) + v4.*(255.*delta.*Sn + (45 - 506.*nu).*Sigman)./(21.*sqrt(6)));
    else 
        error('m=%d correction not implemented',m)
    end

    r_omg = dyn.r.*(1+dyn.chi1(3).*dyn.r.^(-3/2)).^(2/3);
    h_newt = 8*sqrt(pi/5).*(vd).^2.*exp(-2i.*dyn.phi);
    h_newt_asymm = h_newt + h_hat;
    %h_spin = DB_nuX1X2(1000).*(8.*v2).*sqrt(pi./5).*h_hat;
    h_spin = h_hat;

    %h_spin_interp = spline(dyn.t,h_spin,t_num);
    h_num_interp = spline(t_num,h_num,dyn.t);

    h_tot = abs(h_num_interp + h_spin);

    %figure
    %DB_plot_rel_diff(vd,dyn.t,v,dyn.t);
    out.htot = h_tot;
    out.hspin = h_spin;
    out.hnewt = h_newt;
    out.hasym = h_newt_asymm;
    out.t = dyn.t;
    
return