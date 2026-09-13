function P_ic = DB_ic(obj)

%==========================================================================
% Calculates initial data
% Takes r0,q,chi1,chi2 as input
% PR: 07/02/2019
%==========================================================================
q = obj.q;
chi1 = obj.chi1;
chi2 = obj.chi2;
r0 = obj.r0;

pphi0 = 4;
r = r0-1.5:0.03:r0+1.5;
pphicirc = r*0;
px0 = r*0;

for i=1:length(r)
    pphicirc(i) = fzero(@(pphi) dHeff_dr(obj,r(i),pphi,q,chi1,chi2),pphi0);
end

py0 = pphicirc./r;

for i=1:length(r)
    x0 = [r(i);0;0];
    p0 = [0;py0(i);0];
    
    [~,~,dH] = DB_Hamiltonian(obj,x0,p0,chi1,chi2);
    F = DB_flux(x0,p0,dH.dp,q,chi1,chi2);
    Fy(i) = F(2);
end

dpphi_dr = DB_D1(pphicirc,r,4)';


% DEBUG
%{
pxspan = linspace(-1,1,100);

%py = linspace*0;
for i=1:length(pxspan)
    PX(i) = first_PA(obj,r(51),py0(51),q,chi1,chi2,pxspan(i),dpphi_dr(51),Fy(51));
    R = [r(51);0;0];
    P = [pxspan(i);py0(51);0];
    
    [~,~,dHeff] = DB_Hamiltonian(obj,R,P,chi1,chi2);
    dHdpphi(i) = (1./r(51)).*dHeff.dp(2,:) ;

end

figure
hold on
plot(pxspan,PX,'DisplayName','first PA');
plot(pxspan,dHdpphi,'DisplayName','dHdpphi')
legend
yline(0);
ylim([-0.02,0.08])
grid on
title('first PA')
pause
%}

for i=1:length(r)
    %pr = @(px) first_PA(obj,r(i),py0(i),q,chi1,chi2,px,dpphi_dr(i),Fy(i));
    %px0(i) = fzero(pr,-10^(-4));

    px0(i) = fzero(@(px) first_PA(obj,r(i),py0(i),q,chi1,chi2,px,dpphi_dr(i),Fy(i)),-1e-4);
end

px_ic = px0(51);
py_ic = py0(51);
pz_ic = 0;

P_ic = [px_ic;py_ic;pz_ic];

return

function dHeff_dr = dHeff_dr(obj,x,pphi,q,chi1,chi2)

r = [x;0;0];
p = [0;pphi./x;0];

[~,~,dHeff] = DB_Hamiltonian(obj,r,p,chi1,chi2);

dHeff_dr = dHeff.dx(1,:) - 1./x.^2.*pphi.*dHeff.dp(2,:);

return

function first_PA = first_PA(obj,x,py,q,chi1,chi2,px,dpphi_dr,Fy)

% check if it needs H or Heff 

r = [x;0;0];
p = [px;py;0];

[~,~,dHeff] = DB_Hamiltonian(obj,r,p,chi1,chi2);

dHeffdphi = (1./x).*(dHeff.dp(2,:) - Fy);
dHeffdpr  = dHeff.dp(1,:);

%first_PA = dHeffdphi + dpphi_dr.*dHeffdpr; %first_PA = 0 > dp_phi/dt = F_phi
first_PA = dpphi_dr.*dHeffdpr - x.*Fy;
return

%{
function F = root2d(x)

F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;

fun = @root2d;
x0 = [0,0];
x = fsolve(fun,x0)
%}