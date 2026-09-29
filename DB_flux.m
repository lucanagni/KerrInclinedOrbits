function [F,dEdt] = DB_flux(x,p,dHdp,q,chi1)

%==========================================================================
% Precessing flux
% NON-resummed
% Cartesian coordinates
% Same as DB_flux2 but keeps the nu corrections
% PN corrections do not work
% PR: 07/02/2019
%==========================================================================

% Useful variables
[nu, X1, X2] = DB_nuX1X2(q);
eulergamma   = 0.57721566490153286061;

[r,phi,theta] = DB_coords_cart2spherical(x(1,:),x(2,:),x(3,:),p(1,:),p(2,:),p(3,:));

vph =-dHdp(1,:).*sin(phi).*csc(theta).*(1./r) + dHdp(2,:).*cos(phi).*csc(theta).*(1./r);
vth = (1./r).*cos(theta).*(dHdp(1,:).*cos(phi)+dHdp(2,:).*sin(phi)) - dHdp(3,:).*(1./r).*sin(theta);

Omg = sqrt(vth.^2+sin(theta).^2.*vph.^2);

v_omg = Omg.^(1/3);

% Angular momentum
L = cross(x,p,1);
l = L./sqrt(dot(L,L,1)); %3x1 column vector

% Newtonian flux
dEdtN = -32./5.*nu.*Omg.^(10./3);

f2    = -1247./336 -35./12.*nu;
f3    = 4.*pi;
f4    = -44711./9072 + 9271./504.*nu + 65./18.*nu.^2;
f5    = -(8191./672 + 583./24.*nu).*pi;
f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma...
    -(134543./7776 - 41./48.*pi.^2).*nu - 94403./3024.*nu.^2 - 775./324.*nu.^3;
fl6   = -1712./105;
f7    = -(16285./504 - 214745./1728.*nu - 193385./3024.*nu.^2).*pi;

% l and p are 3xN with chi1 a fixed 3x1 broadcast against them, so
% dot() (which needs matching sizes) is replaced by an explicit
% broadcast-multiply-then-sum.
f3so  = -0.25.*(11.*X1 + 5.*X2).*X1.*sum(l.*chi1,1);

%%%% FIXME
%dEdt = dEdtN;
dEdt = dEdtN.*(DB_pade(3,4,f2,f3,f4,f5,f6,fl6,f7,0,0,v_omg)+f3so.*v_omg.^3);

Fspin = (61.*X1 + 48.*X2).*X1.*sum(p.*chi1,1);

% Total flux
F     = 1./(Omg.*vecnorm(L,2,1)).*dEdt.*p + 8./15.*nu.*Omg.^(8./3)./dot(L,L,1)./vecnorm(x,2,1).*Fspin.*L; %3x1 column vector

end
