function [rc,drc] = DB_rc(x,q,chi1,chi2)

%==========================================================================
% Damour-Balmelli centrifugal radius
% Cartesian coordinates
% PR: 01/02/2019
%==========================================================================

% Useful variables
[~, X1, X2] = DB_nuX1X2(q);

a0   = X1.*chi1 + X2.*chi2;

% Centrifugal radius
rc = sqrt(dot(x,x) + dot(a0,a0).*(1+2./sqrt(dot(x,x))));

% Derivatives
drc.dx = x./rc.*(1 - dot(a0,a0)./dot(x,x).^(3/2));

drc.dchi1 = X1.*a0./rc.*(1 + 2./sqrt(dot(x,x)));
drc.dchi2 = X2.*a0./rc.*(1 + 2./sqrt(dot(x,x)));

end