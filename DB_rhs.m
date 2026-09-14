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
chi1 = obj.chi1; % fixed parameter, not evolved

% allocate memory
n    = length(Y(1));
dydt = zeros(6,n);

% unwrap Y
x     = Y(1);
y     = Y(2);
z     = Y(3);
px    = Y(4);
py    = Y(5);
pz    = Y(6);

r    = [x;y;z];
p    = [px;py;pz];

[~,~,dHeff] = DB_Hamiltonian_Kerr(obj,r,p,chi1);

if obj.geodesics
    F = [0;0;0];
elseif obj.flux_KOS
    F = DB_flux_KOS(obj,r,p,dHeff.dp,q,chi1,0*chi1);
elseif obj.flux_nucorrections == 1
    F = DB_flux(r,p,dHeff.dp,q,chi1);
else
    F = DB_flux2(r,p,dHeff.dp,q,chi1); %3x1 column vector
end

dydt(1) = dHeff.dp(1,:);
dydt(2) = dHeff.dp(2,:);
dydt(3) = dHeff.dp(3,:);

dydt(4) = -dHeff.dx(1,:) + F(1,:);
dydt(5) = -dHeff.dx(2,:) + F(2,:);
dydt(6) = -dHeff.dx(3,:) + F(3,:);

return

