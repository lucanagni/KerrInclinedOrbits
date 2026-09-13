function DB_circEJ

%==========================================================================
% E-J curve for DB circular equatorial dynamics
% Cartesian coordinates
% Hard-coded parameters
% PR: 04/02/2019
%==========================================================================

q    = 1;
nu   = q/(1+q)^2;
chi1 = [0.1;0  ;0];
chi2 = [0  ;0.6;0];

pphi0  = 5;
N      = 100;
xmin   = 4;
xmax   = 34;
dx     = (xmax - xmin)/N;
x      = xmin:dx:xmax;

Eb       = zeros(1,N+1);
pphicirc = zeros(1,N+1);
for j=1:N+1
    pphicirc(j) = fzero(@(pphi) dH_dr(x(j),pphi,q,chi1,chi2),pphi0);
    r           = [x(j);0;0];
    p           = [0;pphicirc(j)./x(j);0];
    [~, H]      = DB_Hamiltonian(r,p,q,chi1,chi2);
    Eb(j)       = sqrt(1 + 2.*nu.*(H-1)) - 1;
end

figure
plot(pphicirc,Eb)
xlabel('$p_\varphi$','FontSize',15,'Interpreter','Latex')
ylabel('$E_b$','FontSize',15,'Interpreter','Latex')

return

function dH_dr = dH_dr(x,pphi,q,chi1,chi2)

r = [x;0;0];
p = [0;pphi./x;0];

[~,dH] = DB_Hamiltonian(r,p,q,chi1,chi2);

dH_dr = dH.dx(1,:) - 1./x.^2.*pphi.*dH.dp(2,:);

return
