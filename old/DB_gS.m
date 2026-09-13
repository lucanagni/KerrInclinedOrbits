function [gS,gSstar,dgS,dgSstar] = DB_gS(x,p,q,chi1,chi2)

%==========================================================================
% Damour-Balmelli gyro-gravitomagnetic functions
% NON-resummed
% Cartesian coordinates
% PR: 01/02/2019
%==========================================================================

% Useful variables
nu = q/(1+q).^2;

% Radial versor
n = x./sqrt(dot(x,x));

% Centrifugal radius
[rc,drc] = DB_rc(x,q,chi1,chi2);

% Gyro-gravitomagnetic functions
gS     = 2 - 27./8.*nu.*dot(n,p).^2 - 5./8.*nu./rc ...
    + 5./8.*nu.*(1 + 7.*nu).*dot(n,p).^4 - (21./2 - 23./8.*nu).*nu.*dot(n,p).^2./rc...
    - (51./4 + 1./8.*nu).*nu./rc.^2;

gSstar =  3./2 - (15./8 + 9./4.*nu).*dot(n,p).^2 - (9./8+3./4.*nu)./rc ...
    + (35./16 + 5./2.*nu + 45./16.*nu.^2).*dot(n,p).^4 ...
    + (69./16 - 9./4.*nu + 57./16.*nu.^2).*dot(n,p).^2./rc...
    - (27./16 + 39./4.*nu + 3./16.*nu.^2)./rc.^2;

% Derivatives
dnpdx = 1./sqrt(dot(x,x)).*(p - n.*dot(n,p));

dgS.dx    = (-27./4.*nu + 5./2.*nu.*(1 + 7.*nu).*dot(n,p).^2 - (21 - 23./4.*nu).*nu./rc).*dot(n,p).*dnpdx ...
    + (5./8.*nu + (21./2 - 23./8.*nu).*nu.*dot(n,p).^2 + (51./2 + nu./4).*nu./rc)./rc.^2.*drc.dx;
dgS.dp    = (-27./4.*nu + 5./2.*nu.*(1 + 7.*nu).*dot(n,p).^2 - (21 - 23./4.*nu).*nu./rc).*dot(n,p).*n;
dgS.dchi1 = (5./8.*nu + (21./2 - 23./8.*nu).*nu.*dot(n,p).^2 + (51./2 + nu./4).*nu./rc)./rc.^2.*drc.dchi1;
dgS.dchi2 = (5./8.*nu + (21./2 - 23./8.*nu).*nu.*dot(n,p).^2 + (51./2 + nu./4).*nu./rc)./rc.^2.*drc.dchi2;

dgSstar.dx    = (-(15./4 + 9./2.*nu).*nu + (35./4 + 10.*nu + 45./4.*nu.^2).*dot(n,p).^2 ...
    + (69./8 - 9./2.*nu + 57./8.*nu.^2)./rc).*dot(n,p).*dnpdx ...
    + ((9./8. + 3./4.*nu) - (69./8 - 9./2.*nu + 57./8.*nu.^2).*dot(n,p).^2 ...
    + (27./8 + 39./2.*nu + 3./8.*nu.^2)./rc)./rc.^2.*drc.dx;
dgSstar.dp    = (-(15./4 + 9./2.*nu).*nu + (35./4 + 10.*nu + 45./4.*nu.^2).*dot(n,p).^2 ...
    + (69./8 - 9./2.*nu + 57./8.*nu.^2)./rc).*dot(n,p).*n;
dgSstar.dchi1 = ((9./8. + 3./4.*nu) - (69./8 - 9./2.*nu + 57./8.*nu.^2).*dot(n,p).^2 ...
    + (27./8 + 39./2.*nu + 3./8.*nu.^2)./rc)./rc.^2.*drc.dchi1;
dgSstar.dchi2 = ((9./8. + 3./4.*nu) - (69./8 - 9./2.*nu + 57./8.*nu.^2).*dot(n,p).^2 ...
    + (27./8 + 39./2.*nu + 3./8.*nu.^2)./rc)./rc.^2.*drc.dchi2;

end