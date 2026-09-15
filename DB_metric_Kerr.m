function [A,Bp,Bnp,Benp,Gs,dA,dBp,dBnp,dBenp,dGs] = DB_metric_Kerr(X,chi1)

%==========================================================================
% Kerr Cartesian metric
%==========================================================================

r = vecnorm(X,2,1);
n = X./r;
a = norm(chi1);

x = X(1,:);
y = X(2,:);
z = X(3,:);

% shorthands
r2 = r.^2;
u  = 1./r;
u3 = u.^3;
u4 = u.^4;
a2 = a.^2;
a4 = a.^4;

% squared centrifugal radius
rc2 = r2 + a2.*(1+2*u);
rc4 = rc2.^2;
rc  = sqrt(rc2);
drc = r./rc.*(1-a2.*u3);

D = r2-2.*r+a2;
dD = 2.*r -2;

% A function
Aeq  = (1-2./rc).*(1+2./rc)./(1+2./r);
dAeq = 2.*r2./rc4.*(1+2.*(r-2).*u3.*a2 + a4.*u4);

% n.chi1 and its square, reused throughout A/B/Gs and their derivatives
% (sum(n.*chi1,1) instead of dot(n,chi1): dot() requires matching sizes,
% but chi1 is a fixed 3x1 broadcast against the 3xN array of points n)
nchi  = sum(n.*chi1,1);
nchi2 = nchi.^2;

% common denominators reused across A, B and Gs (and their derivatives)
QN   = r2 + nchi2;            % r^2 + (n.chi1)^2
QD   = r2.*rc2 + D.*nchi2;    % r^2*rc^2 + D*(n.chi1)^2
Fben = r2 + 2.*r + nchi2;     % numerator shorthand for Benp

A    = Aeq.*QN.*rc2./QD;
Bp   = r2./QN;
Bnp  = (D-r2)./QN;
Benp = -r2.*Fben./(QN.*QD);

% Gs function
Gs = 2.*r./QD;

%derivatives

% dot(d n/dxi, chi1), reused in every derivative below. Written as an
% explicit linear combination of the (fixed, scalar) chi1 components
% instead of building a length-3 vector and calling dot() on it, so this
% vectorizes over N points without any reshaping.
dndx_chi = chi1(1).*r.^(-3).*(r2-x.^2) + chi1(2).*(-1).*r.^(-3).*x.*y + chi1(3).*(-1).*r.^(-3).*x.*z;
dndy_chi = chi1(1).*(-1).*r.^(-3).*x.*y + chi1(2).*r.^(-3).*(r2-y.^2) + chi1(3).*(-1).*r.^(-3).*y.*z;
dndz_chi = chi1(1).*(-1).*r.^(-3).*x.*z + chi1(2).*(-1).*r.^(-3).*y.*z + chi1(3).*r.^(-3).*(r2-z.^2);

QNinv2 = QN.^(-2);
QDinv2 = QD.^(-2);

% ---- dA ----
Ecoef     = (-2).*D.*drc.*r+(-2).*D.*rc+dD.*r.*rc+2.*rc.^3;
Fcoef_dir = (-2).*D.*drc+dD.*rc;
Gcoef     = 2.*r.^3.*rc.*(D-rc.^2);
Kcoef     = r.*Ecoef + Fcoef_dir.*nchi2;

BigTermX = x.*nchi.*Kcoef + Gcoef.*dndx_chi;
BigTermY = y.*nchi.*Kcoef + Gcoef.*dndy_chi;
BigTermZ = z.*nchi.*Kcoef + Gcoef.*dndz_chi;

Apref = rc./r.*QDinv2;

dAdx = Apref.*(dAeq.*rc.*x.*QN.*QD - Aeq.*nchi.*BigTermX);
dAdy = Apref.*(dAeq.*rc.*y.*QN.*QD - Aeq.*nchi.*BigTermY);
dAdz = Apref.*(dAeq.*rc.*z.*QN.*QD - Aeq.*nchi.*BigTermZ);

% ---- dBp ----
dBpdx = 2.*nchi.*QNinv2.*(x.*nchi - r2.*dndx_chi);
dBpdy = 2.*nchi.*QNinv2.*(y.*nchi - r2.*dndy_chi);
dBpdz = 2.*nchi.*QNinv2.*(z.*nchi - r2.*dndz_chi);

% ---- dBnp ----
Bnp_c1 = r.*(dD.*r-2.*D);
Bnp_c2 = dD-2.*r;
Bnp_c3 = 2.*r.*(r2-D);
Kbnp   = Bnp_c1 + Bnp_c2.*nchi2;

dBnpdx = (1./r).*QNinv2.*(x.*Kbnp + Bnp_c3.*nchi.*dndx_chi);
dBnpdy = (1./r).*QNinv2.*(y.*Kbnp + Bnp_c3.*nchi.*dndy_chi);
dBnpdz = (1./r).*QNinv2.*(z.*Kbnp + Bnp_c3.*nchi.*dndz_chi);

% ---- dBenp ----
Benp_c1 = 2.*drc.*r2.*rc + 2.*r.*rc.^2;
Benp_c3 = 2.*D.*r;
Kbenp   = Benp_c1 + dD.*nchi2;

termAx = x.*(1+r) + r.*nchi.*dndx_chi;
termAy = y.*(1+r) + r.*nchi.*dndy_chi;
termAz = z.*(1+r) + r.*nchi.*dndz_chi;

termBx = r2.*dndx_chi - x.*nchi;
termBy = r2.*dndy_chi - y.*nchi;
termBz = r2.*dndz_chi - z.*nchi;

termCx = x.*Kbenp + Benp_c3.*nchi.*dndx_chi;
termCy = y.*Kbenp + Benp_c3.*nchi.*dndy_chi;
termCz = z.*Kbenp + Benp_c3.*nchi.*dndz_chi;

dBenpdx = QNinv2.*QDinv2.*((-2).*r.*QN.*QD.*termAx + 2.*nchi.*QD.*Fben.*termBx + r.*QN.*Fben.*termCx);
dBenpdy = QNinv2.*QDinv2.*((-2).*r.*QN.*QD.*termAy + 2.*nchi.*QD.*Fben.*termBy + r.*QN.*Fben.*termCy);
dBenpdz = QNinv2.*QDinv2.*((-2).*r.*QN.*QD.*termAz + 2.*nchi.*QD.*Fben.*termBz + r.*QN.*Fben.*termCz);

% ---- dGs ----
Gs_c1 = 2.*drc.*r.^3.*rc + r2.*rc.^2;
Gs_c2 = dD.*r-D;
Gs_c3 = 2.*D.*r2;
Kgs   = Gs_c1 + Gs_c2.*nchi2;

dGsdx = (-2./r).*QDinv2.*(x.*Kgs + Gs_c3.*nchi.*dndx_chi);
dGsdy = (-2./r).*QDinv2.*(y.*Kgs + Gs_c3.*nchi.*dndy_chi);
dGsdz = (-2./r).*QDinv2.*(z.*Kgs + Gs_c3.*nchi.*dndz_chi);

dA.dx = [dAdx;dAdy;dAdz];

dBp.dx = [dBpdx;dBpdy;dBpdz];

dBnp.dx = [dBnpdx;dBnpdy;dBnpdz];

dBenp.dx = [dBenpdx;dBenpdy;dBenpdz];

dGs.dx = [dGsdx;dGsdy;dGsdz];

return
