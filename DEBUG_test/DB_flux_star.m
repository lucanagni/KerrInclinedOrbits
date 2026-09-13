function F = DB_flux_star(x,ps,dHdp,q,chi1,chi2)

%==========================================================================
% Precessing flux
% NON-resummed
% Cartesian coordinates
% Omega?
% PN corrections do not work
% PR: 07/02/2019
%==========================================================================

% Useful variables
nu         = q/(1+q).^2;
eulergamma = 0.57721566490153286061;
X1         = q./(1+q);
X2         = 1./(1+q);

% Angular velocity
radxy = sqrt(x(1,:).^2 + x(2,:).^2); 
theta = atan(radxy./x(3,:));
phi   = atan(x(2,:)./x(1,:));

vr    = dHdp(1,:).*sin(theta).*cos(phi)...
    + dHdp(2,:).*sin(theta).*sin(phi) + dHdp(3,:).*cos(theta);

omega = 1./sqrt(dot(x,x)).*sqrt(dot(dHdp,dHdp) - vr.^2);
v_omg = omega.^(1/3);

[T,Tinv,~,~] = DB_Tmatrix(x,chi1);

% Angular momentum
L = cross(x,Tinv*ps);
l = L./sqrt(dot(L,L));

% Newtonian flux
dEdtN = -32./5.*nu.*v_omg.^10;

f2    = -1247./336 -35./12.*nu;
f3    = 4.*pi;
f4    = -44711./9072 + 9271./504.*nu + 65./18.*nu.^2;
f5    = -(8191./672 + 583./24.*nu).*pi;
f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma...
    -(134543./7776 - 41./48.*pi.^2).*nu - 94403./3024.*nu.^2 - 775./324.*nu.^3;
fl6   = -1712./105;
f7    = -(16285./504 - 214745./1728.*nu - 193385./3024.*nu.^2).*pi;

f3so  = -0.25.*(11.*X1 + 5.*X2).*X1.*dot(l,chi1)...
    -0.25.*(11.*X2 + 5.*X1).*X2.*dot(l,chi2);
f4ss  = nu./48.*(289.*dot(l,chi1).*dot(l,chi2) - 103.*dot(chi1,chi2));


dEdt  = dEdtN.*(1 + f2.*v_omg.^2 + (f3 + f3so).*v_omg.^3 + (f4 + f4ss).*v_omg.^4 ...
    + f5.*v_omg.^5 + (f6 + fl6.*log(4.*v_omg)).*v_omg.^6 + f7.*v_omg.^7);

%%%% FIXME
%dEdt = 0;
dEdt = dEdtN.*DB_pade(3,4,f2,f3,f4,f5,f6,fl6,f7,f3so,f4ss,v_omg);
%a = sqrt(dot(chi2,chi2));

%dEdt = dEdtN*(DB_flux_spinless(v_omg)+f3so.*v_omg.^3+f4ss.*v_omg.^4);

Fspin = (61.*X1 + 48.*X2).*X1.*dot(Tinv*ps,chi1) + (61.*X2 + 48.*X1).*X2.*dot(Tinv*ps,chi2);

% Total flux
F     = 1./omega./sqrt(dot(L,L)).*dEdt.*Tinv*ps...
    + 8./15.*nu.^2.*v_omg.^8./dot(L,L)./sqrt(dot(x,x)).*Fspin.*L;

end