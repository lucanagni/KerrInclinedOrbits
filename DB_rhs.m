function dydt = DB_rhs(obj, ~, Y)

%==========================================================================
% Right-hand side of DB Hamilton's equations
% Cartesian coordinates
% PR: 07/02/2019
%
% SA: 30/09/2024
% Updated with class-logic
%==========================================================================

q = obj.q;

% allocate memory
n    = length(Y(1));
dydt = zeros(12,n);

% unwrap Y
x     = Y(1);
y     = Y(2);
z     = Y(3);
px    = Y(4);
py    = Y(5);
pz    = Y(6);
chi1x = Y(7);
chi1y = Y(8);
chi1z = Y(9);
chi2x = Y(10);
chi2y = Y(11);
chi2z = Y(12);

r    = [x;y;z];
p    = [px;py;pz];
chi1 = [chi1x;chi1y;chi1z];
chi2 = [chi2x;chi2y;chi2z];

[~,~,dHeff] = DB_Hamiltonian(obj,r,p,chi1,chi2);
%if contains(obj.hamiltonian, 'kerr')
%    dH = dHeff;
%end

if obj.geodesics
    F = [0;0;0];
elseif obj.flux_KOS
    F = DB_flux_KOS(obj,r,p,dHeff.dp,q,chi1,chi2);
elseif obj.flux_nucorrections == 1
    F = DB_flux(r,p,dHeff.dp,q,chi1,chi2);
else
    F = DB_flux2(r,p,dHeff.dp,q,chi1,chi2); %3x1 column vector
end

dydt(1) = dHeff.dp(1,:);
dydt(2) = dHeff.dp(2,:);
dydt(3) = dHeff.dp(3,:);

dydt(4) = -dHeff.dx(1,:) + F(1,:);
dydt(5) = -dHeff.dx(2,:) + F(2,:);
dydt(6) = -dHeff.dx(3,:) + F(3,:);

dchi1dt = cross(dHeff.dchi1,chi1);
dydt(7) = dchi1dt(1,:);
dydt(8) = dchi1dt(2,:);
dydt(9) = dchi1dt(3,:);

dchi2dt = cross(dHeff.dchi2,chi2);
dydt(10) = dchi2dt(1,:);
dydt(11) = dchi2dt(2,:);
dydt(12) = dchi2dt(3,:);

return

