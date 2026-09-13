function [A,dA,B,dB] = DB_metric_Kerr_eq(x,a)

%==========================================================================
% Kerr equatorial metric in polar coords
% input in cartesian coordinates
%==========================================================================

r = norm(x);

% shorthands
r2 = r.^2;
u  = 1./r;
u2 = u.^2;
u3 = u.^3;
u4 = u.^4;
a2 = a.^2;
a4 = a.^4;

% squared centrifugal radius
rc2 = r2 + a2.*(1+2*u);
rc4 = rc2.^2;
rc  = sqrt(rc2);

% A function
A  = (1-2./rc).*(1+2./rc)./(1+2./r);

dA = 2.*r2./rc4.*(1+2.*(r-2).*u3.*a2 + a4.*u4);

% B function
B  = 1./(1-2*u + a2.*u2);
dB = 2.*(a2+(-1).*r).*r.*(a2+((-2)+r).*r).^(-2);

return