function [r,ph,th,pr,pph,pth] = DB_coords_cart2spherical(x,y,z,px,py,pz)

%==========================================================================
% Convert cartesian coords and momenta to spherical
% To be tested
%==========================================================================

x2 = x.*x;
y2 = y.*y;
z2 = z.*z;

r  = sqrt(x2 + y2 + z2);
th = atan2( sqrt(x2 + y2), z);
ph = unwrap(atan2(y , x));

costh = cos(th);
sinth = sin(th);
cosph = cos(ph);
sinph = sin(ph);

c   = (px.*cosph + py.*sinph);

pr  = c.*sinth + costh.*pz;
pph = r.*sinth.*(py.*cosph - px.*sinph); 
pth = r.*(-pz.*sinth + costh.*c);

return