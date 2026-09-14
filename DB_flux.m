function [F,dEdt] = DB_flux(x,p,dHdp,q,chi1)

%==========================================================================
% Precessing flux
% NON-resummed
% Cartesian coordinates
% Omega?
% PN corrections do not work
% PR: 07/02/2019
%==========================================================================

% Useful variables
nu         = q/(1+q).^2;
eulergamma = 0.57721566490153286061;
X1         = q./(1+q);
X2         = 1./(1+q);
  
[r,phi,theta] = DB_coords_cart2spherical(x(1),x(2),x(3),p(1),p(2),p(3));

vr    = dHdp(1,:).*sin(theta).*cos(phi)+ dHdp(2,:).*sin(theta).*sin(phi) + dHdp(3,:).*cos(theta);
omega = 1./r.*sqrt(dot(dHdp,dHdp) - vr.^2);

%vph =-dHdp(1).*sin(phi).*csc(theta).*(1./r) + dHdp(2).*cos(phi).*csc(theta).*(1./r);
%vth = (1./r).*cos(theta).*(dHdp(1).*cos(phi)+dHdp(2).*sin(phi)) - dHdp(3).*(1./r).*sin(theta));

%omega = sqrt(vth.^2+sin(theta).^2.*vph.^2);
v_omg = omega.^(1/3);

% Angular momentum
L = cross(x,p);
l = L./sqrt(dot(L,L));

% Newtonian flux
dEdtN = -32./5.*nu.*v_omg.^10; %!!!!! NU SHOULD BE SQUARED !!!!!

f2    = -1247./336 -35./12.*nu;
f3    = 4.*pi;
f4    = -44711./9072 + 9271./504.*nu + 65./18.*nu.^2;
f5    = -(8191./672 + 583./24.*nu).*pi;
f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma...
    -(134543./7776 - 41./48.*pi.^2).*nu - 94403./3024.*nu.^2 - 775./324.*nu.^3;
fl6   = -1712./105;
f7    = -(16285./504 - 214745./1728.*nu - 193385./3024.*nu.^2).*pi;

f3so  = -0.25.*(11.*X1 + 5.*X2).*X1.*dot(l,chi1);

dEdt  = dEdtN.*(1 + f2.*v_omg.^2 + (f3 + f3so).*v_omg.^3 + f4.*v_omg.^4 ...
    + f5.*v_omg.^5 + (f6 + fl6.*log(4.*v_omg)).*v_omg.^6 + f7.*v_omg.^7);

%%%% FIXME
%dEdt = dEdtN;
dEdt = dEdtN.*(DB_pade(3,4,f2,f3,f4,f5,f6,fl6,f7,0,0,v_omg)+f3so.*v_omg.^3);

Fspin = (61.*X1 + 48.*X2).*X1.*dot(p,chi1);

% Total flux
F     = 1./(omega.*norm(L)).*dEdt.*p...
    + 8./15.*nu.*v_omg.^8./dot(L,L)./sqrt(dot(x,x)).*Fspin.*L;

end