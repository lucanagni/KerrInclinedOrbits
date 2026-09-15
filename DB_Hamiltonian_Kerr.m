function [Heff,Horb,dHeff,dHorb,dHso] = DB_Hamiltonian_Kerr(obj,X,p,chi1)

%==========================================================================
% Kerr Cartesian Hamiltonian
%==========================================================================

%error('DB_Hamiltonian_Kerr not implemented')

[A,Bp,Bnp,Benp,Gs,dA,dBp,dBnp,dBenp,dGs] = DB_metric_Kerr(X,chi1);

r = vecnorm(X,2,1);
L = cross(X,p,1);

x  = X(1,:);
y  = X(2,:);
z = X(3,:);
px = p(1,:);
py = p(2,:);
pz = p(3,:);

r2 = r.^2;

pdotp = dot(p,p,1);
pdotX = dot(X,p,1);      % = r * dot(n,p)
LchiX = sum(L.*chi1,1);  % = r * dot(cross(n,p),chi1); chi1 is a fixed 3x1
                         % broadcast against the 3xN array of points, so
                         % dot() (which needs matching sizes) won't do.
Lchi  = LchiX./r;      % = dot(cross(n,p),chi1)

% S = r^2 * (Horb-argument), so that A.*S == r^2.*Horb.^2. This avoids
% recomputing the whole Horb radicand (and its normalized n/n.p/L terms)
% again in every derivative below.
S  = r2.*(1+Bp.*pdotp+Benp.*Lchi.^2) + Bnp.*pdotX.^2;
AS = A.*S;
Horb = sqrt(AS)./r;
my_zeros = Horb.*0;

if obj.so_coupling
    Hspin = Gs.*LchiX;
else
    Hspin = my_zeros;
end

Heff = Horb + Hspin;

dAdx = dA.dx(1,:);
dAdy = dA.dx(2,:);
dAdz = dA.dx(3,:);
dBpdx = dBp.dx(1,:);
dBpdy = dBp.dx(2,:);
dBpdz = dBp.dx(3,:);
dBnpdx = dBnp.dx(1,:);
dBnpdy = dBnp.dx(2,:);
dBnpdz = dBnp.dx(3,:);
dBenpdx = dBenp.dx(1,:);
dBenpdy = dBenp.dx(2,:);
dBenpdz = dBenp.dx(3,:);
dGsdx = dGs.dx(1,:);
dGsdy = dGs.dx(2,:);
dGsdz = dGs.dx(3,:);

% dot(cross(e_i,p),chi1) for i=x,y,z -- reused both for the L.chi1
% x/y/z-derivatives below and for dHso.dx (they don't depend on X).
Vx = (-1).*pz.*chi1(2)+py.*chi1(3);
Vy = pz.*chi1(1)+(-1).*px.*chi1(3);
Vz = (-1).*py.*chi1(1)+px.*chi1(2);

% dot(cross(n,e_i),chi1) * r, for i=x,y,z -- reused both for the L.chi1
% p-derivatives below and for dHso.dp (they don't depend on p).
Wpx = chi1(2).*z+(-1).*chi1(3).*y;
Wpy = chi1(3).*x+(-1).*chi1(1).*z;
Wpz = chi1(1).*y+(-1).*chi1(2).*x;

InnerX = dBpdx.*pdotp.*r2 + dBnpdx.*pdotX.^2 + dBenpdx.*r2.*Lchi.^2 + 2.*Benp.*Lchi.*(r.*Vx-Lchi.*x);
InnerY = dBpdy.*pdotp.*r2 + dBnpdy.*pdotX.^2 + dBenpdy.*r2.*Lchi.^2 + 2.*Benp.*Lchi.*(r.*Vy-Lchi.*y);
InnerZ = dBpdz.*pdotp.*r2 + dBnpdz.*pdotX.^2 + dBenpdz.*r2.*Lchi.^2 + 2.*Benp.*Lchi.*(r.*Vz-Lchi.*z);

Tx = (-2).*Bnp.*pdotX.*(x.*pdotX-px.*r2) + r2.*InnerX;
Ty = (-2).*Bnp.*pdotX.*(y.*pdotX-py.*r2) + r2.*InnerY;
Tz = (-2).*Bnp.*pdotX.*(z.*pdotX-pz.*r2) + r2.*InnerZ;

dHorbdx = (dAdx.*r2.*S + A.*Tx)./(2.*r2.^2.*Horb);
dHorbdy = (dAdy.*r2.*S + A.*Ty)./(2.*r2.^2.*Horb);
dHorbdz = (dAdz.*r2.*S + A.*Tz)./(2.*r2.^2.*Horb);

dHorbdpx = A.*(Bnp.*x.*pdotX + r2.*Bp.*px + r.*Benp.*Lchi.*Wpx)./(r2.*Horb);
dHorbdpy = A.*(Bnp.*y.*pdotX + r2.*Bp.*py + r.*Benp.*Lchi.*Wpy)./(r2.*Horb);
dHorbdpz = A.*(Bnp.*z.*pdotX + r2.*Bp.*pz + r.*Benp.*Lchi.*Wpz)./(r2.*Horb);

%chi1 = chi1./cos(pi/6); %DEBUG
dHsodx = Gs.*Vx + dGsdx.*LchiX;
dHsody = Gs.*Vy + dGsdy.*LchiX;
dHsodz = Gs.*Vz + dGsdz.*LchiX;

dHsodpx = Gs.*Wpx;
dHsodpy = Gs.*Wpy;
dHsodpz = Gs.*Wpz;
%chi1 = a; %DEBUG

dHorb.dx = [dHorbdx;  dHorbdy;  dHorbdz];
dHorb.dp = [dHorbdpx; dHorbdpy; dHorbdpz];
dHso.dx = [dHsodx;  dHsody;  dHsodz];
dHso.dp = [dHsodpx; dHsodpy; dHsodpz];

if obj.so_coupling
    dHeff.dx = dHorb.dx + dHso.dx;
    dHeff.dp = dHorb.dp + dHso.dp;
else
    dHeff.dx = dHorb.dx;
    dHeff.dp = dHorb.dp;
end

%{
dH.H    = H;
dH.Heff = Heff;
dH.Horb = Horb;
dH.Gs   = Gs;
dH.A    = A;
dH.Bp   = Bp;
dH.Bnp  = Bnp;
dH.dp = dHeff.dp;
%}

end
