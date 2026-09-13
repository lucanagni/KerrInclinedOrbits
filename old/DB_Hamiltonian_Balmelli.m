function [Heff,H,dHeff,dH] = DB_Hamiltonian_Balmelli(x,p,q,chi1,chi2)
%function [H,Horb,Gs,A,Bp,Bnp,Heff,dHeff,dH] = DB_Hamiltonian(x,p,q,chi1,chi2)

%==========================================================================
% Damour-Balmelli Hamiltonian
% Cartesian coordinates
% PR: 31/01/2019
%==========================================================================

% Useful variables
[nu, X1, X2] = DB_nuX1X2(q);

a0    = X1.*chi1 + X2.*chi2;
S     = X1.^2.*chi1 + X2.^2.*chi2;
Sstar = nu.*(chi1+chi2);

Delta = dot(x,x) - 2*sqrt(dot(x,x)) + dot(a0,a0);
R4    = dot(x,x).^2 + dot(x,x).*dot(a0,a0) + 2.*sqrt(dot(x,x)).*dot(a0,a0);

% Radial versor
n = x./sqrt(dot(x,x));

% Centrifugal radius
[rc,drc] = DB_rc(x,q,chi1,chi2);

% Gyro-gravitomagnetic functions
[gS,gSstar,dgS,dgSstar] = DB_gS(x,p,q,chi1,chi2);

% Angular momentum
L = cross(x,p);

% Metric functions
[A,Bp,Bnp,dA,dBp,dBnp] = DB_metric_Balmelli(x,q,chi1,chi2);
Q4  = 2.*(4 - 3.*nu)./dot(x,x).*dot(n,p).^4; %FIXME?

% Orbital Hamiltonian
Horb = sqrt(A.*(1 + Bp.*dot(p,p) + Bnp.*dot(n,p).^2 ...
    - 1./(1 + dot(n,a0).^2./dot(x,x)).*((dot(x,x) + 2.*sqrt(dot(x,x)) + dot(n,a0).^2)...
    ./(R4 + Delta.*dot(n,a0).^2)).*(dot(cross(n,p),a0)).^2 + Q4));

% Spin-orbit Hamiltonian
Hso = 1./sqrt(dot(x,x))./rc.^2./(1 + Delta.*dot(n,a0).^2./dot(x,x)./rc.^2).*gS.*dot(L,S)...
    + 1./rc.^3.*gSstar.*dot(L,Sstar);

% Total Hamiltonian
Heff = Horb + Hso;

% Total Hamiltonian
H = 1./nu.*sqrt(1+2.*nu.*(Heff-1));


% Derivatives
dDelta.dx    = 2.*(sqrt(dot(x,x))-1).*n;
dDelta.dchi1 = 2.*a0.*X1;
dDelta.dchi2 = 2.*a0.*X2;

dR4.dx    = 2.*(2.*dot(x,x).^(3/2) + sqrt(dot(x,x)).*dot(a0,a0) + dot(a0,a0)).*n;
dR4.dchi1 = 2.*a0.*X1.*(dot(x,x) + 2.*sqrt(dot(x,x)));
dR4.dchi2 = 2.*a0.*X2.*(dot(x,x) + 2.*sqrt(dot(x,x)));

N1 = dot(x,x) + 2.*sqrt(dot(x,x)) + dot(n,a0).^2;
D1 = 1 + dot(n,a0).^2./dot(x,x);
D2 = R4 + Delta.*dot(n,a0).^2;

dN1dx    = 2.*x + 2.*x./sqrt(dot(x,x)) ...
    + 2.*dot(n,a0)./sqrt(dot(x,x)).*(a0 - 2.*x.*dot(n,a0)./sqrt(dot(x,x)));
dN1dchi1 = 2.*dot(n,a0).*X1.*n;
dN1dchi2 = 2.*dot(n,a0).*X2.*n;

dD1dx    = 2.*dot(n,a0)./dot(x,x).^(3/2).*(a0 - 2.*x.*dot(n,a0)./sqrt(dot(x,x)));
dD1dchi1 = 2.*dot(n,a0).*X1.*n./dot(x,x);
dD1dchi2 = 2.*dot(n,a0).*X2.*n./dot(x,x);

dD2dx    = dR4.dx + dDelta.dx.*dot(n,a0).^2 ...
    + 2.*Delta.*dot(n,a0)./sqrt(dot(x,x)).*(a0 - x.*dot(n,a0)./sqrt(dot(x,x)));
dD2dchi1 = dR4.dchi1 + dDelta.dchi1.*dot(n,a0).^2 + 2.*Delta.*dot(n,a0).*X1.*n;
dD2dchi2 = dR4.dchi2 + dDelta.dchi2.*dot(n,a0).^2 + 2.*Delta.*dot(n,a0).*X2.*n;

fact       = N1./D1./D2;
dfactdx    = dN1dx./D1./D2 - N1./D1.^2./D2.*dD1dx - N1./D1./D2.^2.*dD2dx;
dfactdchi1 = dN1dchi1./D1./D2 - N1./D1.^2./D2.*dD1dchi1 - N1./D1./D2.^2.*dD2dchi1;
dfactdchi2 = dN1dchi2./D1./D2 - N1./D1.^2./D2.*dD1dchi2 - N1./D1./D2.^2.*dD2dchi2;

dQ4.dx    = 8.*(4 - 3.*nu)./dot(x,x).*dot(n,p).^3./sqrt(dot(x,x)).*(p - dot(n,p).*n);
dQ4.dp    = 8.*(4 - 3.*nu)./dot(x,x).*dot(n,p).^3.*n;
dQ4.dchi1 = 0;
dQ4.dchi2 = 0;

dHorb.dx    = 1./2./Horb.*(dA.dx.*(1 + Bp.*dot(p,p) + Bnp.*dot(n,p).^2 ...
    - fact.*(dot(cross(n,p),a0)).^2 + Q4)...
    + A.*(dBp.dx.*dot(p,p) + dBnp.dx.*dot(n,p).^2 ...
    + 2.*Bnp.*dot(n,p)./sqrt(dot(x,x)).*(p - n.*dot(n,p)) - dfactdx.*(dot(cross(n,p),a0)).^2 ...
    - 2.*fact.*dot(cross(n,p),a0)./sqrt(dot(x,x)).*(cross(p,a0) - (dot(cross(n,p),a0)).*n)...
    + dQ4.dx));
dHorb.dp    = 1./2./Horb.*A.*(2.*Bp.*p + 2.*Bnp.*dot(n,p).*n ...
    + 2.*fact.*dot(cross(n,p),a0).*cross(n,a0) + dQ4.dp);
dHorb.dchi1 = 1./2./Horb.*(dA.dchi1.*(1 + Bp.*dot(p,p) + Bnp.*dot(n,p).^2 ...
    - 1./(1 + dot(n,a0).^2./dot(x,x)).*((dot(x,x) + 2.*sqrt(dot(x,x)) + dot(n,a0).^2)...
    ./(R4 + Delta.*dot(n,a0).^2)).*(dot(cross(n,p),a0)).^2 + Q4)...
    + A.*(dBp.dchi1.*dot(p,p) + dBnp.dchi1.*dot(n,p).^2 ...
    - dfactdchi1.*(dot(cross(n,p),a0)).^2 - 2.*fact.*dot(cross(n,p),a0).*X1.*cross(n,p))...
    + dQ4.dchi1);
dHorb.dchi2 = 1./2./Horb.*(dA.dchi2.*(1 + Bp.*dot(p,p) + Bnp.*dot(n,p).^2 ...
    - 1./(1 + dot(n,a0).^2./dot(x,x)).*((dot(x,x) + 2.*sqrt(dot(x,x)) + dot(n,a0).^2)...
    ./(R4 + Delta.*dot(n,a0).^2)).*(dot(cross(n,p),a0)).^2 + Q4)...
    + A.*(dBp.dchi2.*dot(p,p) + dBnp.dchi2.*dot(n,p).^2 ...
    - dfactdchi2.*(dot(cross(n,p),a0)).^2 - 2.*fact.*dot(cross(n,p),a0).*X2.*cross(n,p))...
    + dQ4.dchi2);

factSO       = 1 + Delta.*dot(n,a0).^2./dot(x,x)./rc.^2;
dfactSOdx    = dDelta.dx.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^2.*(1./sqrt(dot(x,x)).*(a0 - dot(n,a0).*n)...
    - x./dot(x,x) - 1./rc.*drc.dx);
dfactSOdchi1 = dDelta.dchi1.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*dot(n,a0).*X1.*n./dot(x,x)./rc.^2 - 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^3.*drc.dchi1;
dfactSOdchi2 = dDelta.dchi2.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*dot(n,a0).*X2.*n./dot(x,x)./rc.^2 - 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^3.*drc.dchi2;

dHso.dx    = 1./rc.^2./sqrt(dot(x,x))./factSO.*(-n./sqrt(dot(x,x)).*gS.*dot(L,S) ...
    - 2./rc.*drc.dx.*gS.*dot(L,S) - 1./factSO.*dfactSOdx.*gS.*dot(L,S)...
    + dgS.dx.*dot(L,S) + gS.*cross(p,S))...
    + 1./rc.^3.*(-3./rc.*drc.dx.*gSstar.*dot(L,Sstar) + dgSstar.dx.*dot(L,Sstar)...
    + gSstar.*cross(p,Sstar));
dHso.dp    = 1./sqrt(dot(x,x))./rc.^2./factSO.*dgS.dp.*dot(L,S) + 1./rc.^3.*dgSstar.dp.*dot(L,Sstar)...
    - 1./sqrt(dot(x,x))./rc.^2./factSO.*gS.*cross(x,S) - 1./rc.^3.*gSstar.*cross(x,Sstar);
dHso.dchi1 = - 2./sqrt(dot(x,x))./rc.^3./factSO.*gS.*dot(L,S).*drc.dchi1 ...
    - 3./rc.^4.*gSstar.*dot(L,Sstar).*drc.dchi1...
    - 1./sqrt(dot(x,x))./rc.^2./factSO.^2.*dfactSOdchi1.*gS.*dot(L,S)...
    + 1./sqrt(dot(x,x))./rc.^2./factSO.*dgS.dchi1.*dot(L,S) + 1./rc.^3.*dgSstar.dchi1.*dot(L,Sstar)...
    + X1.^2./sqrt(dot(x,x))./rc.^2./factSO.*gS.*L + X1.*X2./rc.^3.*gSstar.*L;
dHso.dchi2 = - 2./sqrt(dot(x,x))./rc.^3./factSO.*gS.*dot(L,S).*drc.dchi2...
    - 3./rc.^4.*gSstar.*dot(L,Sstar).*drc.dchi2...
    - 1./sqrt(dot(x,x))./rc.^2./factSO.^2.*dfactSOdchi2.*gS.*dot(L,S)...
    + 1./sqrt(dot(x,x))./rc.^2./factSO.*dgS.dchi2.*dot(L,S) + 1./rc.^3.*dgSstar.dchi2.*dot(L,Sstar)...
    + X1.*X2./sqrt(dot(x,x))./rc.^2./factSO.*gS.*L + X2.^2./rc.^3.*gSstar.*L;

dHeff.dx    = dHorb.dx + dHso.dx;
dHeff.dp    = dHorb.dp + dHso.dp;
dHeff.dchi1 = dHorb.dchi1 + dHso.dchi1;
dHeff.dchi2 = dHorb.dchi2 + dHso.dchi2;

dH.dx    = 1./nu./H.*dHeff.dx;
%dH.dp    = 1./nu./H.*dHeff.dp;
dH.dchi1 = 1./nu./H.*dHeff.dchi1;
dH.dchi2 = 1./nu./H.*dHeff.dchi2;


dH.H    = H;
dH.Heff = Heff;
dH.Horb = Horb;
dH.Gs   = gS;
dH.A    = A;
dH.Bp   = Bp;
dH.Bnp  = Bnp;
dH.dp = 1./nu./H.*dHeff.dp; %added for analogy with other hamiltonians. Needed for Omg

%fprintf('%15.12f \n',H);

end
