function [F,dEdt] = DB_flux2(x,p,dHdp,q,chi1,chi2)

%==========================================================================
% Precessing flux
% NON-resummed
% Cartesian coordinates
% Omega?
% PN corrections do not work
% PR: 07/02/2019
%==========================================================================

% Useful variables
nu         = DB_nuX1X2(q);
eulergamma = 0.57721566490153286061;

[r,phi,theta] = DB_coords_cart2spherical(x(1),x(2),x(3),p(1),p(2),p(3));

%vr    = dHdp(1,:).*sin(theta).*cos(phi)+ dHdp(2,:).*sin(theta).*sin(phi) + dHdp(3,:).*cos(theta);
%omega = 1./r.*sqrt(dot(dHdp,dHdp) - vr.^2);

vph =-dHdp(1).*sin(phi).*csc(theta).*(1./r) + dHdp(2).*cos(phi).*csc(theta).*(1./r);
vth = (1./r).*cos(theta).*(dHdp(1).*cos(phi)+dHdp(2).*sin(phi)) - dHdp(3).*(1./r).*sin(theta);

Omg = sqrt(vth.^2+sin(theta).^2.*vph.^2);

%fprintf('Omg - omega = %f\n',abs(Omg-omega))
v_omg = Omg.^(1/3);

% Angular momentum
L = cross(x,p);
l = L./sqrt(dot(L,L)); %3x1 column vector

% Newtonian flux
dEdtN = -32./5.*nu.*Omg.^(10./3);

f2    = -1247./336;
f3    = 4.*pi;
f4    = -44711./9072;
f5    = -(8191./672).*pi;
f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma;
fl6   = -1712./105;
f7    = -(16285./504).*pi;

f3so  = -11./4.*dot(l,chi1);

%%%% FIXME
%dEdt = dEdtN;
dEdt = dEdtN.*(DB_pade(3,4,f2,f3,f4,f5,f6,fl6,f7,0,0,v_omg)+f3so.*v_omg.^3);
%a = sqrt(dot(chi2,chi2));

Fspin = 61.*dot(p,chi1);
% Total flux
F     = 1./(Omg.*norm(L)).*dEdt.*p + 8./15.*nu.*Omg.^(8./3)./dot(L,L)./norm(x).*Fspin.*L; %3x1 column vector

end
