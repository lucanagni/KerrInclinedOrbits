function DB_circ_solver

%==========================================================================
% First try at implementing BD dynamics
% Circular dynamics
% Cartesian coordinates
% Hard-coded parameters
% PR: 31/01/2019
%==========================================================================
global q

% Time grid
Tmin =  0;
Tmax =  50000;
dt   =  0.5;

reltol           =  1e-8;
abstol           =  1e-10;

% Binary parameters
q  = 0.0001;
chi1 = [0;0; 0.2];
chi2 = [0.4;0;-0.2];

% Initial radius
r0 = 3.2733644609e+01;

% Initial phase space variables
py0 = DB_circpy0(r0,q,chi1,chi2);
r = [r0;0;0];
p = [0;py0;0];

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

rend = 5;

%-----------------------
% Integrate the dynamics
%-----------------------
options = odeset('events',@(T,Y) DB_ode_stop(T,Y,rend),'RelTol',reltol,'AbsTol',abstol);
[T,Y]   = ode113(@DB_circ_rhs,Tmin:dt:Tmax,y0,options);

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




figure
plot3(x,y,z)
xlabel('$x$','FontSize',15,'Interpreter','Latex')
ylabel('$y$','FontSize',15,'Interpreter','Latex')
zlabel('$z$','FontSize',15,'Interpreter','Latex')

%%{
radius = sqrt(x.^2 + y.^2 + z.^2);
radxy  = sqrt(x.^2 + y.^2);
phi    = unwrap(2.*atan(y./x))./2;
theta  = unwrap(2.*atan(radxy./z))./2;

figure
plot(T,radius)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$r$','FontSize',15,'Interpreter','Latex')

figure
plot(T,phi)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$\varphi$','FontSize',15,'Interpreter','Latex')

figure
plot(T,theta)
xlabel('$t$','FontSize',15,'Interpreter','Latex')
ylabel('$\theta$','FontSize',15,'Interpreter','Latex')
%}

chi1vect = [chi1x chi1y chi1z]';
chi2vect = [chi2x chi2y chi2z]';

figure
hold on
plot(T,chi1vect(1,:))
plot(T,chi1vect(2,:))
plot(T,chi1vect(3,:))
leg = legend('$\chi_{1,x}$','$\chi_{1,y}$','$\chi_{1,z}$','location','northeast');
set(leg,'Interpreter','Latex','FontSize',10);

figure
hold on
plot(T,chi2vect(1,:))
plot(T,chi2vect(2,:))
plot(T,chi2vect(3,:))
leg = legend('$\chi_{2,x}$','$\chi_{2,y}$','$\chi_{2,z}$','location','northeast');
set(leg,'Interpreter','Latex','FontSize',10);

%{
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


return

