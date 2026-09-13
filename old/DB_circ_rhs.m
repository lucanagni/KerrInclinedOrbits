function dydt = DB_circ_rhs(t,yy)

%==========================================================================
% Right-hand side of DB Hamilton's equations
% Circular dynamics
% Cartesian coordinates
% PR: 06/02/2019
%==========================================================================

global q

% allocate memory
dydt = zeros(12,1);

% unwrap yy
x     = yy(1);
y     = yy(2);
z     = yy(3);
px    = yy(4);
py    = yy(5);
pz    = yy(6);
chi1x = yy(7);
chi1y = yy(8);
chi1z = yy(9);
chi2x = yy(10);
chi2y = yy(11);
chi2z = yy(12);

r    = [x;y;z];
p    = [px;py;pz];
chi1 = [chi1x;chi1y;chi1z];
chi2 = [chi2x;chi2y;chi2z];

[~,~,~,dH] = DB_Hamiltonian(r,p,q,chi1,chi2);

dydt(1) = dH.dp(1,:);
dydt(2) = dH.dp(2,:);
dydt(3) = dH.dp(3,:);

dydt(4) = -dH.dx(1,:);
dydt(5) = -dH.dx(2,:);
dydt(6) = -dH.dx(3,:);

dchi1dt = cross(dH.dchi1,chi1);
dydt(7) = dchi1dt(1,:);
dydt(8) = dchi1dt(2,:);
dydt(9) = dchi1dt(3,:);

dchi2dt = cross(dH.dchi2,chi2);
dydt(10) = dchi2dt(1,:);
dydt(11) = dchi2dt(2,:);
dydt(12) = dchi2dt(3,:);

return

