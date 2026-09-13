function [A,Bp,Bnp,dA,dBp,dBnp] = DB_metric_Balmelli(x,q,chi1,chi2)

%==========================================================================
% Damour-Balmelli metric functions
% NON-resummed
% Cartesian coordinates
% PR: 01/02/2019
%==========================================================================

% Useful variables
pi2 = pi.^2;
pi4 = pi2.^2;
eulergamma = 0.57721566490153286061;

[nu, X1, X2] = DB_nuX1X2(q);

nu2 = nu.*nu;

a0   = X1.*chi1 + X2.*chi2;

Delta = dot(x,x) - 2*sqrt(dot(x,x)) + dot(a0,a0);

% Radial versor
n = x./sqrt(dot(x,x));

% Centrifugal radius
[rc,drc] = DB_rc(x,q,chi1,chi2);
uc = 1./rc;

% Padé coefficients
a5c0 = -4237./60 + 2275./512.*pi.^2 + 256./5.*log(2) + 128./5.*eulergamma;
a5c1 = -221./6 + 41./32.*pi.^2;
a5   = a5c0 + nu.*a5c1;
a6   = 3097.3.*nu.^2 - 1330.6.*nu + 81.38;

loguc  = log(uc);
a5tot  = a5  + 64/5*loguc;
a6tot  = a6  + (-7004/105 - 144/5*nu).*loguc;
a5tot2 = a5tot.^2;

N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*pi2)))...
    ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))...
    ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))...
    ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*pi2 - 123*nu*pi2))...
    ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*pi2)) ...
    + nu*(2048*(5582 + 9*nu) - 834432*pi2 + 15129*pi4)))./(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*pi2)...
    - (-3776 + 123*pi2)*(-3008 - 96*nu + 123*pi2))))./(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));

% Orbital metric functions
Num  = 1 + N1.*uc;
Den  = 1 + D1.*uc + D2.*uc.^2 + D3.*uc.^3 + D4.*uc.^4 + D5.*uc.^5;
Aorb = Num./Den;

Dorb = (1 + 6.*nu.*uc.^2 + 2.*(26 - 3.*nu).*nu.*uc.^3).^(-1);

% Kerr-like metric function
Aeq  = (1+2.*uc)./(1+2./sqrt(dot(x,x))).*Aorb;

AnuK0   = Aeq.*(1 + dot(n,a0).^2./dot(x,x))./(1+Delta.*dot(n,a0).^2./(dot(x,x).*rc.^2));
BpnuK0  = 1./(1 + dot(n,a0).^2./dot(x,x));
BnpnuK0 = 1./(1 + dot(n,a0).^2./dot(x,x)).*(Aeq./Dorb.*rc.^2./dot(x,x) - 1);

% DB precessing corrections
AchiQ  = (3.*X1 - nu./2).*nu.*dot(chi1,chi1) + (2 - nu).*nu.*dot(chi1,chi2)...
    + (3.*X2 - nu./2).*nu.*dot(chi2,chi2);
BchiQ  = (9.*X1 - 15./4.*nu).*nu.*dot(chi1,chi1) + (6 + 9./2.*nu).*nu.*dot(chi1,chi2)...
    + (9.*X2 - 15./4.*nu).*nu.*dot(chi2,chi2);

AnchiQ = (2.*X1 + 5./2.*nu).*nu.*dot(n,chi1).^2 + (3 - 7.*nu).*nu.*dot(n,chi1).*dot(n,chi2)...
    + (2.*X2 + 5./2.*nu).*nu.*dot(n,chi2).^2;
BnchiQ = (9.*X1 - 15./4.*nu).*nu.*dot(n,chi1).^2 + (6 + 9./2.*nu).*nu.*dot(n,chi1).*dot(n,chi2)...
    + (9.*X2 - 15./4.*nu).*nu.*dot(n,chi2).^2;

% Total metric function
A   = AnuK0 + (AchiQ - AnchiQ)./dot(x,x).^2;
Bp  = BpnuK0 - BnchiQ./sqrt(dot(x,x)).^3;
Bnp = BnpnuK0 + BchiQ./sqrt(dot(x,x)).^3;

% Derivatives
dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
    ./(7.*power(1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*uc);
dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
    ./(7.*power(1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*uc);
dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
    ./(7.*power(1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*uc);
dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
    ./(7.*power(1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*uc);
dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
    ./(7.*power(1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*uc);
dD5 = (nu*(-8400*nu*(-24*(a6 - (4*loguc*(1751 + 756*nu))/105.).*(1536 + nu*(-3776 + 123*pi2)) ...
    + nu*(-2304*power(a5 + (64*loguc)/5.,2) + 96*(a5 + (64*loguc)/5.).*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) ...
    - (1536*loguc*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2))).*(4128768*loguc*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*pi2) - 756*nu*(1536 + nu*(-3776 + 123*pi2))))))...
    ./(2625.*power(-768 + nu*(3584 + 24*(a5 + (64*loguc)/5.) - 123*pi2),2).*uc);

% First derivative
dNum  = dN1.*uc + N1;
dDen  = D1 + uc.*(dD1 + 2*D2) + uc.^2.*(dD2 + 3*D3) + uc.^3.*(dD3 + 4*D4) ...
    + uc.^4.*(dD4 + 5*D5) + dD5.*uc.^5;

% derivative of A function with respect to uc
prefactor = Aorb./(Num.*Den);
dAorb_duc = prefactor.*(dNum.*Den - dDen.*Num);

dDorb_duc = -6.*nu.*uc.*(2 + (26 - 3.*nu).*uc).*Dorb.^2;

% Derivative of A with respect to r
dAorb.dx    = -uc.^2.*dAorb_duc.*drc.dx;
dAorb.dchi1 = -uc.^2.*dAorb_duc.*drc.dchi1;
dAorb.dchi2 = -uc.^2.*dAorb_duc.*drc.dchi2;

dDorb.dx    = -uc.^2.*dDorb_duc.*drc.dx;
dDorb.dchi1 = -uc.^2.*dDorb_duc.*drc.dchi1;
dDorb.dchi2 = -uc.^2.*dDorb_duc.*drc.dchi2;

dAeq.dx    = -2.*uc.^2.*drc.dx./(1 + 2./sqrt(dot(x,x))).*Aorb ...
    + 2.*(1 + 2.*uc)./(1 + 2./sqrt(dot(x,x))).^2.*n./dot(x,x).*Aorb...
    + (1 + 2.*uc)./(1 + 2./sqrt(dot(x,x))).*dAorb.dx;
dAeq.dchi1 = -2.*uc.^2.*drc.dchi1./(1 + 2./sqrt(dot(x,x))).*Aorb...
    + (1+2.*uc)./(1+2./sqrt(dot(x,x))).*dAorb.dchi1;
dAeq.dchi2 = -2.*uc.^2.*drc.dchi2./(1 + 2./sqrt(dot(x,x))).*Aorb...
    + (1+2.*uc)./(1+2./sqrt(dot(x,x))).*dAorb.dchi2;

dDelta.dx    = 2.*x - 2.*x./sqrt(dot(x,x));
dDelta.dchi1 = 2.*a0.*X1;
dDelta.dchi2 = 2.*a0.*X2;

numAnuK0       = 1 + dot(n,a0).^2./dot(x,x);
dnumAnuK0dx    = 2.*dot(n,a0)./dot(x,x).^(3/2).*(a0 - 2.*x.*dot(n,a0)./sqrt(dot(x,x)));
dnumAnuK0dchi1 = 2.*dot(n,a0).*X1.*n./dot(x,x);
dnumAnuK0dchi2 = 2.*dot(n,a0).*X2.*n./dot(x,x);

denAnuK0       = 1 + Delta.*dot(n,a0).^2./dot(x,x)./rc.^2;
ddenAnuK0dx    = dDelta.dx.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^2.*(1./sqrt(dot(x,x)).*(a0 - dot(n,a0).*n)...
    - x./dot(x,x) - 1./rc.*drc.dx);
ddenAnuK0dchi1 = dDelta.dchi1.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*dot(n,a0).*X1.*n./dot(x,x)./rc.^2 - 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^3.*drc.dchi1;
ddenAnuK0dchi2 = dDelta.dchi2.*dot(n,a0).^2./dot(x,x)./rc.^2 ...
    + 2.*dot(n,a0).*X2.*n./dot(x,x)./rc.^2 - 2.*Delta.*dot(n,a0).^2./dot(x,x)./rc.^3.*drc.dchi2;

dAnuK0.dx    = AnuK0.*(dAeq.dx./Aeq + dnumAnuK0dx./numAnuK0 - ddenAnuK0dx./denAnuK0);
dAnuK0.dchi1 = AnuK0.*(dAeq.dchi1./Aeq + dnumAnuK0dchi1./numAnuK0 - ddenAnuK0dchi1./denAnuK0);
dAnuK0.dchi2 = AnuK0.*(dAeq.dchi2./Aeq + dnumAnuK0dchi2./numAnuK0 - ddenAnuK0dchi2./denAnuK0);

dBpnuK0.dx    = -1./numAnuK0.^2.*dnumAnuK0dx;
dBpnuK0.dchi1 = -1./numAnuK0.^2.*dnumAnuK0dchi1;
dBpnuK0.dchi2 = -1./numAnuK0.^2.*dnumAnuK0dchi2;

dBnpnuK0.dx    = dBpnuK0.dx.*(Aeq./Dorb.*rc.^2./dot(x,x) - 1) ...
    + BpnuK0.*(Aeq./Dorb.*rc.^2./dot(x,x)).*(dAeq.dx./Aeq - dDorb.dx./Dorb ...
    + 2.*drc.dx./rc - 2.*x./dot(x,x));
dBnpnuK0.dchi1 = dBpnuK0.dchi1.*(Aeq./Dorb.*rc.^2./dot(x,x) - 1) ...
    + BpnuK0.*(Aeq./Dorb.*rc.^2./dot(x,x)).*(dAeq.dchi1./Aeq - dDorb.dchi1./Dorb ...
    + 2.*drc.dchi1./rc);
dBnpnuK0.dchi2 = dBpnuK0.dchi2.*(Aeq./Dorb.*rc.^2./dot(x,x) - 1) ...
    + BpnuK0.*(Aeq./Dorb.*rc.^2./dot(x,x)).*(dAeq.dchi2./Aeq - dDorb.dchi2./Dorb ...
    + 2.*drc.dchi2./rc);

dAchiQ.dchi1 = 2.*(3.*X1 - nu./2).*nu.*chi1 + (2 - nu).*nu.*chi2;
dAchiQ.dchi2 = 2.*(3.*X1 - nu./2).*nu.*chi2 + (2 - nu).*nu.*chi1;

dBchiQ.dchi1 = 2.*(9.*X1 - 15./4.*nu).*nu.*chi1 + (6 + 9./2.*nu).*nu.*chi2;
dBchiQ.dchi2 = 2.*(9.*X1 - 15./4.*nu).*nu.*chi2 + (6 + 9./2.*nu).*nu.*chi1;

dAnchiQ.dchi1 = 2.*(2.*X1 + 5./2.*nu).*nu.*dot(n,chi1).*n + (3 - 7.*nu).*nu.*n.*dot(n,chi2);
dAnchiQ.dchi2 = 2.*(2.*X1 + 5./2.*nu).*nu.*dot(n,chi2).*n + (3 - 7.*nu).*nu.*n.*dot(n,chi1);
dAnchiQ.dx    = 2.*(2.*X1 + 5./2.*nu).*nu.*dot(n,chi1).*(chi1./sqrt(dot(x,x)) - x.*dot(n,chi1)./dot(x,x)) ...
    + (3 - 7.*nu).*nu.*(dot(n,chi1).*(chi2./sqrt(dot(x,x)) - x.*dot(n,chi2)./dot(x,x))...
    + dot(n,chi2).*(chi1./sqrt(dot(x,x)) - x.*dot(n,chi1)./dot(x,x)))...
    + 2.*(2.*X2 + 5./2.*nu).*nu.*dot(n,chi2).*(chi2./sqrt(dot(x,x)) - x.*dot(n,chi2)./dot(x,x));

dBnchiQ.dchi1 = 2.*(9.*X1 - 15./4.*nu).*nu.*dot(n,chi1).*n + (6 + 9./2.*nu).*nu.*n.*dot(n,chi2);
dBnchiQ.dchi2 = 2.*(9.*X1 - 15./4.*nu).*nu.*dot(n,chi2).*n + (6 + 9./2.*nu).*nu.*n.*dot(n,chi1);
dBnchiQ.dx    = (9.*X1 - 15./4.*nu).*nu.*dot(n,chi1).*(chi1./sqrt(dot(x,x)) - x.*dot(n,chi1)./dot(x,x))...
    + (6 + 9./2.*nu).*nu.*(dot(n,chi1).*(chi2./sqrt(dot(x,x)) - x.*dot(n,chi2)./dot(x,x)) ...
    + dot(n,chi2).*(chi1./sqrt(dot(x,x)) - x.*dot(n,chi1)./dot(x,x)))...
    + (9.*X2 - 15./4.*nu).*nu.*dot(n,chi2).*(chi2./sqrt(dot(x,x)) - x.*dot(n,chi2)./dot(x,x));

dA.dx    = dAnuK0.dx - dAnchiQ.dx./dot(x,x).^2 - 4.*(AchiQ - AnchiQ)./dot(x,x).^3.*x;
dA.dchi1 = dAnuK0.dchi1 + (dAchiQ.dchi1 - dAnchiQ.dchi1)./dot(x,x).^2;
dA.dchi2 = dAnuK0.dchi2 + (dAchiQ.dchi2 - dAnchiQ.dchi2)./dot(x,x).^2;

dBp.dx    = dBpnuK0.dx - dBnchiQ.dx./sqrt(dot(x,x)).^3 + 3.*BnchiQ./dot(x,x).^2.*n;
dBp.dchi1 = dBpnuK0.dchi1 - dBnchiQ.dchi1./sqrt(dot(x,x)).^3;
dBp.dchi2 = dBpnuK0.dchi2 - dBnchiQ.dchi2./sqrt(dot(x,x)).^3;

dBnp.dx    = dBnpnuK0.dx - 3.*BchiQ./dot(x,x).^2.*n;
dBnp.dchi1 = dBnpnuK0.dchi1 + dBchiQ.dchi1./sqrt(dot(x,x)).^3;
dBnp.dchi2 = dBnpnuK0.dchi2 + dBchiQ.dchi2./sqrt(dot(x,x)).^3;

end