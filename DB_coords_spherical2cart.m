function [x,y,z,px,py,pz] = DB_coords_spherical2cart(r,ph,th,pr,pph,pth)

%==========================================================================
% Convert spherical coords and momenta to cartesian
% To be tested
%==========================================================================

costh = cos(th);
sinth = sin(th);
cosph = cos(ph);
sinph = sin(ph);

x = r.*sinth.*cosph;
y = r.*sinth.*sinph;
z = r.*costh;

d  = 1./(r.*sinth);
px = pr.*sinth.*cosph - pph.*sinph.*d + (1./r).*costh.*cosph.*pth;
py = pr.*sinth.*sinph + pph.*cosph.*d + (1./r).*costh.*sinph.*pth;
pz = pr.*costh        - (1./r).*sinth.*pth;

return