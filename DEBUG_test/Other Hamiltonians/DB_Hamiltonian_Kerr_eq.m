function [Heff,H,dHeff,dH] = DB_Hamiltonian_Kerr_eq(X,P,a)

%==========================================================================
% Kerr equatorial Hamiltonian
% Input in cartesian coordinates
%
% SA: 2/31/2024
% +++ To check/clean/improve +++
%==========================================================================

% Cartesian to polar 
[r,~,~,pr,pph,~] = DB_coords_cart2spherical(X(1),X(2),X(3),P(1),P(2),P(3));

% Kerr metric functions
[A, dA, B, dB] = DB_metric_Kerr_eq(X,a);

prstar = sqrt(A./B).*pr;

% shorthands
a2      = a.^2;
%a4      = a2.^2;
u       = 1./r;
u2      = u.^2;
u3      = u.*u2;
r2      = r.^2;
r3      = r.*r2;
pph2    = pph.^2;
prstar2 = prstar.^2;

%------------------------------------------
% Centrifugal radius. It is rc=r when r-> 0
%------------------------------------------
rc2 = r2 + a2.*(1+2*u);
uc2 = 1./rc2;
%uc3 = uc2.^(3/2);
rc  = sqrt(rc2);
drc = r./rc.*(1-a2.*u3);

% Gyro-gravitomagnetic ratio and its derivatives
Gs  =  2*u.*uc2;
dGs = -2*(a2+3*r2)./((r3+a2*(2+r)).^2);

%=================
% Kerr Hamiltonian
%=================
% Pure orbital part of the Kerr Hamiltonian (at pi/2)
Horb = sqrt(A.*(1+pph2.*uc2) + prstar2);

% Kerr Hamiltonian (orbital plane)
Hso  = Gs.*a.*pph;
Heff = Horb + Hso;

% derivative of the effective Hamiltonian with respect to r
%dHeffdr = dGs_dr.*a.*pph + 1./(2*Horb).*(dA_dr.*(1 + pph2.*uc2) - 2.*A.*uc3.*drc_dr.*pph2);

% +++ TO CHECK +++
x  = X(1);
y  = X(2);
px = P(1);
py = P(2);

dHorbdx = (1/2).*B.^(-2).*r.^(-4).*rc.^(-3).*(A.*(1+rc.^(-2).*(py.*x+(-1).* ...
    px.*y).^2+B.^(-1).*r.^(-2).*(px.*x+py.*y).^2)).^(-1/2).*(B.*dA.* ...
    r.*rc.*x.*(rc.^2.*(px.*x+py.*y).^2+B.*r.^2.*(rc.^2+(py.*x+(-1).* ...
    px.*y).^2))+A.*((-1).*dB.*r.*rc.^3.*x.*(px.*x+py.*y).^2+(-2).*B.* ...
    rc.^3.*(px.*x+py.*y).*(px.*((-1).*r.^2+x.^2)+py.*x.*y)+2.*B.^2.* ...
    r.^2.*(py.*x+(-1).*px.*y).*(drc.*r.*x.*((-1).*py.*x+px.*y)+rc.*( ...
    py.*(r.^2+(-1).*x.^2)+px.*x.*y))));
dHorbdy = (1/2).*B.^(-2).*r.^(-4).*rc.^(-3).*(A.*(1+rc.^(-2).*(py.*x+(-1).* ...
    px.*y).^2+B.^(-1).*r.^(-2).*(px.*x+py.*y).^2)).^(-1/2).*(B.*dA.* ...
    r.*rc.*y.*(rc.^2.*(px.*x+py.*y).^2+B.*r.^2.*(rc.^2+(py.*x+(-1).* ...
    px.*y).^2))+A.*((-1).*dB.*r.*rc.^3.*y.*(px.*x+py.*y).^2+(-2).*B.* ...
    rc.^3.*(px.*x+py.*y).*(px.*x.*y+py.*((-1).*r.^2+y.^2))+2.*B.^2.* ...
    r.^2.*((-1).*py.*x+px.*y).*(drc.*r.*y.*(py.*x+(-1).*px.*y)+rc.*( ...
    py.*x.*y+px.*(r.^2+(-1).*y.^2)))));
dHorbdz = 0;

dHorbdpx = A.*B.^(-1).*r.^(-2).*rc.^(-2).*(B.*r.^2.*y.*((-1).*py.*x+px.*y)+ ...
    rc.^2.*x.*(px.*x+py.*y)).*(A.*(1+rc.^(-2).*(py.*x+(-1).*px.*y).^2+ ...
    B.^(-1).*r.^(-2).*(px.*x+py.*y).^2)).^(-1/2);
dHorbdpy = A.*B.^(-1).*r.^(-2).*rc.^(-2).*(B.*r.^2.*x.*(py.*x+(-1).*px.*    y)+ ...
    rc.^2.*y.*(px.*x+py.*y)).*(A.*(1+rc.^(-2).*(py.*x+(-1).*px.*y).^2+ ...
    B.^(-1).*r.^(-2).*(px.*x+py.*y).^2)).^(-1/2);
dHorbdpz = 0;

dHsodx = a.*r.^(-2).*(dGs.*r.*x.*(py.*x+(-1).*px.*y)+Gs.*(py.*(r.^2+(-1).* ...
    x.^2)+px.*x.*y));
dHsody = (-1).*a.*r.^(-2).*(dGs.*r.*y.*((-1).*py.*x+px.*y)+Gs.*(py.*x.*y+ ...
    px.*(r.^2+(-1).*y.^2)));
dHsodz = 0;

dHsodpx = -a.*Gs.*y;
dHsodpy = +a.*Gs.*x;
dHsodpz = 0;

dHorb.dx = [dHorbdx;  dHorbdy;  dHorbdz];
dHorb.dp = [dHorbdpx; dHorbdpy; dHorbdpz];
dHso.dx  = [dHsodx;   dHsody;   dHsodz];
dHso.dp  = [dHsodpx;  dHsodpy;  dHsodpz];


dHeff.dx = dHorb.dx + dHso.dx;
dHeff.dp = dHorb.dp + dHso.dp;

% do not evolve spin
dHeff.dchi1 = 0*dHeff.dx;
dHeff.dchi2 = dHeff.dchi1; 

my_zeros = Heff*0;
H = my_zeros + 1;

dH.dx    = my_zeros;
dH.dp    = my_zeros;
dH.dchi1 = my_zeros;
dH.dchi2 = my_zeros;

dH.H    = H;
dH.Heff = Heff;
dH.Horb = Horb;
dH.Gs   = Gs;
dH.A    = A;
dH.Bp   = B;
dH.Bnp  = B*0;
dH.dp = dHeff.dp; %%added for analogy with 3D hamiltonians. Needed for Omg

return



