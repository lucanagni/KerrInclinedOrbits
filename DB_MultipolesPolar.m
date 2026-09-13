function [Ulm,Vlm] = DB_MultipolesPolar(DB,l,m,flags)
    % h21 only

    % For Amplitude
    Omega = DB.Omg;
    r = DB.r;
    iota = (pi/2-min(DB.th))
    phi = atan2(DB.y.*cos(iota),DB.x);
    v_phi = r.*Omega;
    x = v_phi.^2;

    AbsCurr = (4/3).*((1/5).*pi).^(1/2).*x.^(3/2).*(2+cos(2.*iota)+cos(4.*iota)+ ...
      (-2).*(1+2.*cos(2.*iota)).*cos(2.*phi).*sin(iota).^2).^(1/2);
    PhaseCurr_bound = atan(cos(2.*iota).*cot(phi).*sec(iota));
    %PhaseCurr_bound = atan2(cos(2.*iota),tan(phi).*cos(iota));
    %PhaseCurr = mod(pi/2-PhaseCurr_bound,2*pi);
    PhaseCurr =pi/2+unwrap(4.*(PhaseCurr_bound))./4;

      %{
    (sqrt(-1)*(-4/3)).*exp(1).^((sqrt(-1)*(-1)).*phi).*((1/5).*pi).^( ...
      1/2).*x.^(3/2).*((-1).*cos(iota)+exp(1).^((sqrt(-1)*2).*phi).*cos( ...
      iota)+(-1).*cos(2.*iota)+(-1).*exp(1).^((sqrt(-1)*2).*phi).*cos( ...
      2.*iota));
      ((1/5).*pi).^(1/2).*(x.^3.*(2+cos(2.*iota)+cos(4.*iota)+(-2).*(1+ ...
        2.*cos(2.*iota)).*cos(2.*phi).*sin(iota).^2)).^(1/2);
    %}

    CurrQuadrupole = AbsCurr.*exp(1i.*PhaseCurr);
    Vlm = sqrt(2).*CurrQuadrupole;
    Ulm = 0.*Vlm;

    figure
    plot(DB.t,PhaseCurr_bound)
    pause

    %figure
    %plot(DB.t,exp(1i.*ArgCurr))
return
