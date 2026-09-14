function P_ic = DB_postspherical_ICs(obj)

%==========================================================================
% Calculates initial data (px0,py0) for quasi-spherical orbits of arbitraty inclination
% Takes instance of DB_class as input
% LN: 09/2025
% SA: updated with class logic
% assumes pth = 0 and phi = 0
% ==========================================================================

q = obj.q;
chi1 = obj.chi1;
chi2 = obj.chi2;
r0 = obj.r0;
th0 = obj.th0;
phi0 = obj.phi0;
iota = obj.iota;
sep = DB_LSSO(chi1(3),iota);

delta = min(1.5,r0-sep);
if delta<0.5
    delta = 0.5;
end

r = linspace(r0-delta,r0+delta,101);
% Allocate memory
py0 = r*0;
Fy = r*0;

for i=1:length(r)
    py0(i) = fzero(@(py) dHeff_dr(obj,r(i),py),1);
end

for i=1:length(r)
    x0 = [r(i).*sin(th0).*cos(phi0);0;r(i).*cos(th0)];
    p0 = [0;py0(i);0];

    [~,~,dHeff] = DB_Hamiltonian_Kerr(obj,x0,p0,chi1);

    F = DB_flux(x0,p0,dHeff.dp,q,chi1,chi2);
    Fy(i) = F(2);
end

dpphi_dr = DB_D1(r.*py0,r,4);

pr_ic = fzero(@(pr) first_PA(obj,r(51),py0(51),chi1*0,pr,dpphi_dr(51),Fy(51)),-1e-4);

px_ic = pr_ic.*sin(th0).*cos(phi0);
py_ic = py0(51);
pz_ic = pr_ic.*cos(th0);

P_ic = [px_ic;py_ic;pz_ic];

return

function dHeff_dr = dHeff_dr(obj, r, py)

if obj.iota > pi/2
    pr_flag = -1;
else
    pr_flag = 1;
end

phi0 = obj.phi0;
th0 = obj.th0;

[x,y,z] = DB_coords_spherical2cart(r,phi0,th0,0,0,0);
[~,~,~,~,pphi] = DB_coords_cart2spherical(x,y,z,0,py,0);

R = [x,y,z];
P = [0;py;0];

obj = obj;                 %DEBUG
[~,~,dHeff] = DB_Hamiltonian_Kerr(obj,R,P,obj.chi1);

sth = sin(th0);
cth = cos(th0);
r2 = r.^2;

dHeff_dr = pr_flag.*dHeff.dx(1,:).*sth + dHeff.dx(3,:).*cth - pr_flag.*dHeff.dp(2,:).*pphi./(r2.*sth);

return

function first_PA = first_PA(obj,x,py,chi1,px,dpphi_dr,Fy)

    r = [x;0;0];
    p = [px;py;0];

    obj = obj;                 %DEBUG
    [~,~,dHeff] = DB_Hamiltonian_Kerr(obj,r,p,chi1);

    dHeffdpr  = dHeff.dp(1,:);

    first_PA = dpphi_dr.*dHeffdpr - x.*Fy;     %first_PA = 0 -> dp_phi/dt = F_phi
return
