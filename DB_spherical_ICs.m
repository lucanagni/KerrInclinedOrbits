function pycirc = DB_spherical_ICs(obj)

%==========================================================================
% Calculates initial quasispherical py0 for spherical orbits of arbitrary inclination
% Takes instance of DB_class as input
% PR: 07/02/2019
% LN: 09/2025
%
% SA: updated with class logic
% assumes pth = 0 and phi = 0
%==========================================================================

py0 = 1;
pycirc = fzero(@(py) dHeff_dr(obj,py),py0);

return

function dHeff_dr = dHeff_dr(obj, py)

if obj.iota > pi/2
    pr_flag = -1;
else
    pr_flag = 1;
end

r = obj.r0;
th = obj.th0;
phi = obj.phi0;

[x,y,z] = DB_coords_spherical2cart(r,phi,th,0,0,0);
[~,~,~,~,pphi] = DB_coords_cart2spherical(x,y,z,0,py,0);

R = [x,y,z];
P = [0;py;0];

[~,~,dHeff] = DB_Hamiltonian_Kerr(obj,R,P,obj.chi1);

sth = sin(th);
cth = cos(th);
r2 = r.^2;

dHeff_dr = pr_flag.*dHeff.dx(1,:).*sth + dHeff.dx(3,:).*cth - pr_flag.*dHeff.dp(2,:).*pphi./(r2.*sth);

return
