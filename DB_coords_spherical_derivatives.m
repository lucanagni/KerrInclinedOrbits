function [dX,dY,dZ,dPx,dPy,dPz] = DB_coords_spherical_derivatives(r,phi,th,pr,pphi,pth)

    N = length(r);
    my_zeros = zeros(N,1);

    dxdr = cos(phi).*sin(th);
    dxdth = r.*cos(th).*cos(phi);
    dxdphi = (-1).*r.*sin(phi).*sin(th);

    dydr = sin(phi).*sin(th);
    dydth = r.*cos(th).*sin(phi);
    dydphi = r.*cos(phi).*sin(th);

    dzdr = cos(th);
    dzdth = (-1).*r.*sin(th);
    dzdphi = my_zeros;

    dpxdr = (-1).*pth.*r.^(-2).*cos(phi).*cos(th)+pphi.*r.^(-2).*csc(th).*sin(phi);
    dpxdth = pr.*cos(phi).*cos(th)+pphi.*r.^(-1).*cot(th).*csc(th).*sin(phi)+( ...
                -1).*pth.*r.^(-1).*cos(phi).*sin(th);
    dpxdphi = (-1).*pphi.*r.^(-1).*cos(phi).*csc(th)+(-1).*pth.*r.^(-1).*cos(th) ...
                .*sin(phi)+(-1).*pr.*sin(phi).*sin(th);
    dpxdpr = cos(phi).*sin(th);
    dpxdpth = r.^(-1).*cos(phi).*cos(th);
    dpxdpphi = (-1).*r.^(-1).*csc(th).*sin(phi);

    dpydr = (-1).*pphi.*r.^(-2).*cos(phi).*csc(th)+(-1).*pth.*r.^(-2).*cos(th).*sin(phi);
    dpydth = (-1).*pphi.*r.^(-1).*cos(phi).*cot(th).*csc(th)+pr.*cos(th).*sin( ...
                phi)+(-1).*pth.*r.^(-1).*sin(phi).*sin(th);
    dpydphi = pth.*r.^(-1).*cos(phi).*cos(th)+(-1).*pphi.*r.^(-1).*csc(th).*sin( ...
                phi)+pr.*cos(phi).*sin(th);
    dpydpr = sin(phi).*sin(th);
    dpydpth = r.^(-1).*cos(th).*sin(phi);
    dpydpphi = r.^(-1).*cos(phi).*csc(th);

    dpzdr = pth.*r.^(-2).*sin(th);
    dpzdth = (-1).*pth.*r.^(-1).*cos(th)+(-1).*pr.*sin(th);
    dpzdphi = my_zeros;
    dpzdpr = cos(th);
    dpzdpth = (-1).*r.^(-1).*sin(th);
    dpzdpphi = my_zeros;

    dX.dr = dxdr;
    dX.dth = dxdth;
    dX.dphi = dxdphi;
    dX.dpr = my_zeros;
    dX.dpth = my_zeros;
    dX.dpphi = my_zeros;

    dY.dr = dydr;
    dY.dth = dydth;
    dY.dphi = dydphi;
    dY.dpr = my_zeros;
    dY.dpth = my_zeros;
    dY.dpphi = my_zeros;

    dZ.dr = dzdr;
    dZ.dth = dzdth;
    dZ.dphi = dzdphi;
    dZ.dpr = my_zeros;
    dZ.dpth = my_zeros;
    dZ.dpphi = my_zeros;

    dPx.dr = dpxdr;
    dPx.dth = dpxdth;
    dPx.dphi = dpxdphi;
    dPx.dpr = dpxdpr;
    dPx.dpth = dpxdpth;
    dPx.dpphi = dpxdpphi;

    dPy.dr = dpydr;
    dPy.dth = dpydth;
    dPy.dphi = dpydphi;
    dPy.dpr = dpydpr;
    dPy.dpth = dpydpth;
    dPy.dpphi = dpydpphi;

    dPz.dr = dpzdr;
    dPz.dth = dpzdth;
    dPz.dphi = dpzdphi;
    dPz.dpr = dpzdpr;
    dPz.dpth = dpzdpth;
    dPz.dpphi = dpzdpphi;
return