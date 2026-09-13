function obj = DB_solver(r0,write)
%write è una variabile che se vale 1 mi crea la cartella in cui scrive la
%dinamica altrimenti no. È un po' rozzo ma era un metodo veloce per non
%finire con l'avere mille cartelle

%==========================================================================
% First try at implementing BD dynamics
% Cartesian coordinates
% Hard-coded parameters
% PR: 07/02/2019
%==========================================================================
global q
global obj
global i
i = 0;

% Time grid
Tmin =  0;
Tmax =  1e6;
dt   =  0.1;

reltol           =  1e-8;
abstol           =  1e-10;

% Mass ratio and spins
q  = 1e-3;
X2  = 1/(1+q);
obj.a = 0.2;

chi1 = [0;0;0];
chi2 = [0;0;obj.a/X2^2];


% Initial conditions
py0 = DB_circpy0(r0,q,chi1,chi2);

r = [r0;0;0];
p = [-0.0001;py0;0];

y0(1) = r(1);
y0(2) = r(2);
y0(3) = r(3);
y0(4) = p(1);
y0(5) = p(2);
y0(6) = p(3);
y0(7) = chi1(1);
y0(8) = chi1(2);
y0(9) = chi1(3);
y0(10) = chi2(1);
y0(11) = chi2(2);
y0(12) = chi2(3);
obj.H = [];
obj.Horb = [];
obj.Gs = [];
obj.A = [];
obj.Bp = [];
obj.Bnp = [];
obj.tEOB = [];

rend = 1+sqrt(1-obj.a^2)+1e-4;
%rend =  X2+sqrt(X2.^2-0.2.^2)+1e-4;
%rend = 5;
%-----------------------
% Integrate the dynamics
%-----------------------
options = odeset('events',@(T,Y) DB_ode_stop(T,Y,rend),'RelTol',reltol,'AbsTol',abstol);
[T,Y]   = ode113(@DB_rhs,Tmin:dt:Tmax,y0,options);

%fprintf('HAM: length T = %c')
%size(T)
% rename output for convenience
x     = Y(:,1);
y     = Y(:,2);
z     = Y(:,3);
px    = Y(:,4);
py    = Y(:,5);
pz    = Y(:,6);
chi1x = Y(:,7);
chi1y = Y(:,8);
chi1z = Y(:,9);
chi2x = Y(:,10);
chi2y = Y(:,11);
chi2z = Y(:,12);

for i=1:length(T)
    r    = [x(i);y(i);z(i)];
    p    = [px(i);py(i);pz(i)];
    chi1 = [chi1x(i), chi1y(i), chi1z(i)];
    chi2 = [chi2x(i), chi2y(i), chi2z(i)];
    [~,~,~,dH] = DB_Hamiltonian(r,p,q,chi1,chi2);
    obj.H(i) = dH.H;
    obj.Horb(i) = dH.Horb;
    obj.A(i) = dH.A;
    obj.Bp(i) = dH.Bp;
    obj.Bnp(i) = dH.Bnp;
    obj.Gs(i) = dH.Gs;
    obj.tEOB(i) = T(i);
end

obj.x = x;
obj.y = y;
obj.z = z;
obj.px = px;
obj.py = py;
obj.pz = pz;
obj.time = T;
obj.r0 = r0;

display('creating plots')

figure('Name','orbit')
plot3(x,y,z)
xlabel('$x$','FontSize',15,'Interpreter','Latex')
ylabel('$y$','FontSize',15,'Interpreter','Latex')
zlabel('$z$','FontSize',15,'Interpreter','Latex')

%%{
radius = sqrt(x.^2 + y.^2 + z.^2);
radxy  = sqrt(x.^2 + y.^2);
phi    = unwrap(2.*atan(y./x))./2;
theta  = unwrap(2.*atan(radxy./z))./2;

figure('Name','angle-radius')
polarplot(theta, radius)
%thetaticks([])
thetaticklabels({})

figure('Name','radius')
plot(T,radius)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$r$','FontSize',15,'Interpreter','Latex')

figure('Name','phi')
plot(T,phi)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$\varphi$','FontSize',15,'Interpreter','Latex')

figure('Name','theta')
plot(T,theta)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$\theta$','FontSize',15,'Interpreter','Latex')

DB_write_dynamics(write);
%{
chi1vect = [chi1x chi1y chi1z]';
chi2vect = [chi2x chi2y chi2z]';

figure('Name','spin components')
hold on
plot(T,chi1vect(1,:))
plot(T,chi1vect(2,:))
plot(T,chi1vect(3,:))
xlabel('$t$','FontSize',15,'Interpreter','Latex')
leg = legend('$\chi_{1,x}$','$\chi_{1,y}$','$\chi_{1,z}$','location','northeast');
set(leg,'Interpreter','Latex','FontSize',10);


figure
hold on
plot(T,chi2vect(1,:))
plot(T,chi2vect(2,:))
plot(T,chi2vect(3,:))
xlabel('$t$','FontSize',15,'Interpreter','Latex')
leg = legend('$\chi_{2,x}$','$\chi_{2,y}$','$\chi_{2,z}$','location','northeast');
set(leg,'Interpreter','Latex','FontSize',10);


rvect = [x y z]';
pvect = [px py pz]';
Lvect = cross(rvect,pvect);
lvect = Lvect./norm(Lvect);

chi1_eff = dot(lvect,chi1vect);
chi2_eff = dot(lvect,chi2vect);

figure('Name','eff spins')
hold on
plot(T,chi1_eff)
plot(T,chi2_eff)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
leg = legend('$\chi_{1,eff}$','$\chi_{2,eff}$','location','northeast');
set(leg,'Interpreter','Latex','FontSize',10);



X1 = q./(1+q);
X2 = 1./(1+q);
rvect    = [x y z]';
pvect    = [px py pz]';
Lvect = cross(rvect,pvect);
Jvect = Lvect + X1.^2.*chi1vect + X2.^2.*chi2vect;
J2    = dot(Jvect,Jvect);

figure
plot(T,Jvect(1,:))
%}
%}
%display('Calling DB_write_dynamics \n')

%fprintf('RHS called %d times. \n Length of T = %d. \n Length of tEOB = %d \n ',i, length(obj.time), length(obj.tEOB));
%{
fileID  = fopen(['T_debug.dat']   ,'w');
for n=1:length(obj.tEOB)
    fprintf(fileID, '%15.12f %15.12f \n',[obj.time(n) obj.tEOB(n)]);
end
for n=(length(obj.tEOB)+1):length(obj.time)
    fprintf(fileID,'%15.12f \n',obj.time(n));
end
fclose(fileID);

%}



return

