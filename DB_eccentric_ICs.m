function [r,p] = DB_eccentric_ICs(obj)

if abs(obj.chi1(1))>1e-10 || abs(obj.chi1(2))>1e-10
    error('Kerr spin not aligned with z-axis!')
elseif norm(obj.chi2)>1e-10
    error('Spinning particle!')
end

a   = obj.chi1(3);
e0  = obj.e0;
sr0 = obj.sr0;
anomaly0 = obj.anomaly0;

[r0_eq, ph0_eq, ~, pr0_eq, pph0_eq] = Polar_EccentricIn(a, e0, sr0, anomaly0);

th0 = obj.th0;
[x,y,z,px,py,pz] = DB_coords_spherical2cart(r0_eq,ph0_eq,th0,pr0_eq,pph0_eq,0);

r = [ x; y; z];
p = [px;py;pz]; %0.235+px for the dyn Simone is running on cluster

return

function [r0, ph0, prstar0, pr0, pph0, E0, Omg0] = Polar_EccentricIn(a, e0, sr0, anomaly0)
r1  = sr0./(1 - e0);
r2  = sr0./(1 + e0);

[A1, ~, ~, rc1] = kerr_metric_equatorial(r1,a);
[A2, ~, ~, rc2] = kerr_metric_equatorial(r2,a);

B1  = A1./rc1.^2;
B2  = A2./rc2.^2;

G1  = 2./r1./rc1.^2*a;
G2  = 2./r2./rc2.^2*a;

dG  = G1 - G2;
dA  = A1 - A2;
dB  = B1 - B2;

A12 = A1 + A2;
B12 = B1 + B2;

% compute constants of motion
pph2 = (A12.*dG.^2 - dA.*dB + dG.*sqrt(4*A1.*A2.*dG.^2 + 2.*dA.*(B12.*dA - A12.*dB)))./(dB.^2 - 2*B12.*dG.^2 + dG.^4);
pph0 = sqrt(pph2);
E0   = sqrt(A1.*(1 + pph2./rc1.^2)) + G1.*pph0;

% get actual initial radius based on mean anomaly
r0 = sr0/(1-e0*cos(anomaly0));
[A0, ~, B0, rc0] = kerr_metric_equatorial(r0,a);
G0 = 2./r0./rc0.^2*a;

% get radial momentum and frequency
if anomaly0>=0 && anomaly0<pi
    pr_sign = -1;
else
    pr_sign = 1;
end
diff2 = (E0- G0.*pph0).^2 - A0*(1+pph2./rc0.^2);
if abs(diff2)<1e-14
    diff2 = 0;
end

prstar0 = pr_sign*sqrt( diff2 );
pr0  = prstar0*sqrt(B0./A0);
Omg0 = G0 + A0.*pph0./rc0.^2./(E0 - G0.*pph0);
ph0  = anomaly0;

return


function [A, dA, B, rc] = kerr_metric_equatorial(r,a)
r2 = r.^2;
u  = 1./r;
u2 = u.^2;
u3 = u.^3;
u4 = u.^4;
a2 = a.^2;
a4 = a.^4;

rc2 = r2 + a2.*(1+2*u);
rc4 = rc2.^2;
rc  = sqrt(rc2);

A  = (1-2./rc).*(1+2./rc)./(1+2./r);
dA = 2.*r2./rc4.*(1+2.*(r-2).*u3.*a2 + a4.*u4);
B  = 1./(1-2*u + a2.*u2);
return
