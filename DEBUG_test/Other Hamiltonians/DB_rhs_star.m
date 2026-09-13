function dydt = DB_rhs_star(obj, ~, Y)

%==========================================================================
% Right-hand side of DB Hamilton's equations
% Cartesian coordinates with tortoise momentum
%
% FIXME (see DB_Hamiltonian_kerr_star)
%==========================================================================

q = obj.q;

% allocate memory
n    = length(Y(1));
dydt = zeros(12,n);

% unwrap Y
x     = Y(1);
y     = Y(2);
z     = Y(3);
pxs   = Y(4);
pys   = Y(5);
pzs   = Y(6);
chi1x = Y(7);
chi1y = Y(8);
chi1z = Y(9);
chi2x = Y(10);
chi2y = Y(11);
chi2z = Y(12);

r    = [x;y;z];
ps   = [pxs;pys;pzs];
chi1 = [chi1x;chi1y;chi1z];
chi2 = [chi2x;chi2y;chi2z];

[T,Tinv,dT,~] = DB_Tmatrix(r,chi1);

[~,~,dHeff,dH] = DB_Hamiltonian(obj,r,ps,chi1,chi2);
if contains(obj.hamiltonian, 'kerr')
    dH = dHeff;
end

if obj.geodesics
    F = [0;0;0];
else
    F = DB_flux2(r,Tinv*ps,dH.dps,q,chi1,chi2);
end

Xdot = dH.dps'*T;
dydt(1) = Xdot(1);
dydt(2) = Xdot(2);
dydt(3) = Xdot(3);

dTijk(:,:,1) = dT.dx;   
dTijk(:,:,2) = dT.dy;
dTijk(:,:,3) = dT.dz;

P = Tinv*ps;

dPsdX = zeros(3);
for i =1:3
    dPsdX = dPsdX+reshape(dTijk(:,i,:),[3,3]).*P(i);
end
%dpsx = dT.dx*P; %attempt to calculate Pdot differently. Same result.
%dpsy = dT.dy*P;
%dpsz = dT.dz*P;
%dPs = [dpsx,dpsy,dpsz];

Pdot = T*(-dH.dx + F) + dPsdX*(dH.dps'*T)';
%Pdot = T*(-dH.dx + F) + dPs*Xdot';
dydt(4) = Pdot(1);
dydt(5) = Pdot(2);
dydt(6) = Pdot(3);

dchi1dt = cross(dH.dchi1,chi1);
dydt(7) = dchi1dt(1,:);
dydt(8) = dchi1dt(2,:);
dydt(9) = dchi1dt(3,:);

dchi2dt = cross(dH.dchi2,chi2);
dydt(10) = dchi2dt(1,:);
dydt(11) = dchi2dt(2,:);
dydt(12) = dchi2dt(3,:);

return