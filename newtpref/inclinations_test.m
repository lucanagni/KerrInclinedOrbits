function inclinations_test(varargin)

% LN: Script originally use to compute multipole. It is now a mess and has been replaced by DB_waveforms_test. I think IT can be eliminated
m = 1;

radii        = [   6.1,     6.1,     6.1,     6.1];
inclinations = [  pi/2,    pi/3,    pi/4,    pi/6];
hl_pi        = ["\pi/2", "\pi/3", "\pi/4", "\pi/6"];

inputDB.chi1 = [0;0;0];

addpath ../

maxtau = -100000;
hl = "";
% ================ ANALYTICAL WAVEFORMS FOR VARIOUS INCLINATIONS ================
figure('Position', [626 142 818 724])
dirs = ["/home/luca/waveforms/K/plunge/th0_90/a02/wf.mat", ...
        "/home/luca/waveforms/K/plunge/th0_60/a02/wf.mat", ...
        "/home/luca/waveforms/K/plunge/th0_45/a02/wf.mat", ...
        "/home/luca/waveforms/K/plunge/th0_30/a02/wf.mat"];
for i=1:length(inclinations)
    load(dirs(i));
    DB = s.dyn;

    %inputDB.r0  = radii(i);
    %inputDB.th0 = inclinations(i);
    %inputDB.flux_nucorrections = 1;
    %DB = DB_class(inputDB);
    %w = generic_wave(DB.t, DB.r, DB.phi, DB.th);
    w = generic_wave_cart(DB,2);

    tLR = tLR_splined(DB);

    range = 1:length(DB.t)-10;
    tau   = DB.t(range)-tLR;
   
    %if tau(1)>maxtau
    %    maxtau = tau(1);
    %end
    maxtau = -420;

    subplot(2,1,1)
    plot(tau, abs(w(range)), 'LineWidth', 1.0)
    %plot(T, abs(w), 'LineWidth', 1.0)
    ylim([0,1.7])

    hold on

    subplot(2,1,2)
    plot(tau, DB.th(range), 'LineWidth', 1.0)
    hold on

    hl(i) = "$\theta_0  ="+hl_pi(i)+"$";
end

ylabs = {'$h_{22}^{\rm Quad}$', '$\theta$'};

for i=1:2
    subplot(2,1,i)
    legend(hl,'location','NorthWest','Box','off','NumColumns',2, ...
              'FontSize', 12, 'Interpreter', 'latex')
    xlim([maxtau,100])
    xlabel('$t$', 'FontSize', 24, 'Interpreter', 'latex')
    ylabel(ylabs{i}, 'FontSize', 24, 'Interpreter', 'latex')
end

% ================ NUMERICAL WAVEFORMS FOR VARIOUS INCLINATIONS ================

%dirs = ["/home/luca/waveforms/K/plunge/th0_90/a02/wf.mat", ...
%        "/home/luca/waveforms/K/plunge/th0_60/a02/wf.mat", ...
%        "/home/luca/waveforms/K/plunge/th0_45/a02/wf.mat", ...
%        "/home/luca/waveforms/K/plunge/th0_30/a02/wf.mat"];
%
%dirs = ["/home/luca/waveforms/S/plunge/th0_90/more_timesteps/wf.mat", ...
%        "/home/luca/waveforms/S/plunge/th0_60/wf.mat", ...
%        "/home/luca/waveforms/S/plunge/th0_45/wf.mat", ...
%        "/home/luca/waveforms/S/plunge/th0_30/wf.mat"];

figure('Position', [626 142 818 724])
for i=1:length(inclinations)
    load(dirs(i));
    DB = s.dyn;
    tLR = tLR_splined(DB);

    w = s.ell(2).emm(3).hlm;
    T = s.ell(2).emm(3).t-tLR;

    range = 1:length(DB.t);
    tau   = DB.t(range)-tLR;
   
    %if tau(1)>maxtau
    %    maxtau = tau(1);
    %end
    maxtau = -220;

    subplot(2,1,1)
    plot(T, abs(w), 'LineWidth', 1.0)
    ylim([0,1.7])
    hold on

    subplot(2,1,2)
    plot(tau, DB.th(range), 'LineWidth', 1.0)
    hold on

    hl(i) = "$\theta_0  ="+hl_pi(i)+"$";
end

ylabs = {'$h_{22}^{\rm Num}$', '$\theta$'};

for i=1:2
    subplot(2,1,i)
    legend(hl,'location','NorthWest','Box','off','NumColumns',2, ...
              'FontSize', 12, 'Interpreter', 'latex')
    xlim([maxtau,100])
    xlabel('$t$', 'FontSize', 24, 'Interpreter', 'latex')
    ylabel(ylabs{i}, 'FontSize', 24, 'Interpreter', 'latex')
    if i==1
      ylim([0,1.5])
    end
end
  
% ================ COMPARE ANALYTICAL AND NUMERICAL ================
if isempty(varargin)
  N = 1;
else
  N = varargin{1};
end
load(dirs(N))

%inputDB.r0  = radii(N);
%inputDB.th0 = inclinations(N);
%inputDB.flux_nucorrections = 1;
%DB = DB_class(inputDB);
DB = s.dyn;

w_anal = generic_wave_cart(DB,m);
tau = DB.t-tLR_splined(DB);

load(dirs(N))
w_num = s.ell(2).emm(m+1).hlm;
T = s.ell(2).emm(m+1).t-tLR_splined(DB);

splined_num = spline(T,w_num,tau);
DeltaA = (abs(w_anal)-abs(splined_num))./abs(splined_num);
DeltaPhi = (-unwrap(angle(w_anal)) + unwrap(angle(-splined_num)));

my_linewidth = 1;
axes_fontsize = 12;
legend_fontsize = axes_fontsize+2;
labels_fontsize = 20;

f = figure;
f.Position(3:4) = f.Position(3:4)*1.2;
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
t_in = -420;
t_end = 0;
label = KerrLabel(DB);

nexttile(tl,[2,1])
set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
xlim([t_in,t_end])
hold on
plot(T,-real(w_num),'LineWidth',my_linewidth,'DisplayName','$\Re[h_{22}]$')
plot(tau,real(w_anal),'LineWidth',my_linewidth,'DisplayName','$\Re[h_{22}^{\rm EOB}]$','LineStyle','--') 
plot(tau,abs(w_anal),'LineWidth',my_linewidth,'DisplayName','$|h_{22}^{\rm EOB}|$','LineStyle','--','Color',[1 0 0 .2])
plot(T,abs(w_num),'LineWidth',my_linewidth,'DisplayName','$|h_{22}|$','Color',[1 0 0 .2])
text(-150,1.5,label,'Interpreter','latex','FontSize',16)

nexttile
set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
xlim([t_in,t_end])
hold on
grid on
plot(tau,DeltaPhi,'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
plot(tau,DeltaA,'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
%yscale log
%yticks(-0.5:0.25:0.5)

xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')

%hold on
%plot(tau,real(w_anal),'DisplayName','Analytical')
%plot(T,-real(w_num),'DisplayName','Numerical')
%legend
%xlim([maxtau,100])

if abs(tLR_splined(DB)-tLR_splined(s.dyn))>1e-5
  fprintf('WARNING: analytical tLR is %.5f while numerical tLR is %.5f',tLR_splined(DB),tLR_splined(s.dyn))
end 

% DEBUG
%r0    = 1.213061319425267e+00;   % 2/sqrt(e);
%Omega = DB.Omg;
%Heff = DB.Heff;
%m=2;
%Tail = DB_Tail(2,m.*Omega,m.*Omega,r0);
%
%
%figure
%plot(DB.t,(DB.Omg))

rmpath ../
return



function out = generic_wave(t, r, phi, th)

phidot  = DB_D1(phi, t, 4);
phi2dot = DB_D1(phidot, t, 4);

rdot   = DB_D1(r, t, 4);
r2dot  = DB_D1(rdot, t, 4);
thdot  = DB_D1(th, t, 4);
th2dot = DB_D1(thdot, t, 4);

% circ term here is (-2).*phidot.^2.*r.^2
out = (1/2).*exp((1i*(-2)).*phi).* ...
  (r.^2.*thdot.^2.*cos(th).^2+r.*(4.*rdot.*thdot+r.*(th2dot+(sqrt( ...
  -1)*(-4)).*phidot.*thdot)).*cos(th).*sin(th)+((1i*(-1)).* ...
  phi2dot.*r.^2+(-2).*phidot.^2.*r.^2+r.*r2dot+(1i*(-4)).* ...
  phidot.*r.*rdot+rdot.^2+(-1).*r.^2.*thdot.^2).*sin(th).^2);

% fix norm 
out = out.*8.*sqrt(pi/5);


return

function out = generic_wave_cart(DB,m)

fprintf('CARTESIAN FORMULA\n')

x = DB.x;
y = DB.y;
z = DB.z;
t = DB.t;
r = DB.r;
rOmg = r.*(1+DB.chi1(3).*DB.r.^(-3/2)).^(2/3);

[x,y,z] = DB_coords_spherical2cart(rOmg,DB.phi,DB.th,0,0,0);

r0    = 1.213061319425267e+00;   % 2/sqrt(e);
Omega = DB.Omg;
Heff = DB.Heff;
Tail = DB_Tail(2,m.*Omega,m.*Omega,r0);

a = DB.chi1(3);
rho = DB_rholm(Omega,a,r);
delta = DB_deltalm(Omega,a);

xdot = DB_D1(x,t,4);
x2dot = DB_D1(xdot,t,4);
x3dot = DB_D1(x2dot,t,4);
ydot = DB_D1(y,t,4);
y2dot = DB_D1(ydot,t,4);
y3dot = DB_D1(y2dot,t,4);
zdot = DB_D1(z,t,4);
z2dot = DB_D1(zdot,t,4);
z3dot = DB_D1(z2dot,t,4);

if m==2
  MassQuadrupole = 2.*(x2dot.*x + xdot.^2) - 2.*(y2dot.*y + ydot.^2) - 2.*1i.*(x2dot.*y + 2.*xdot.*ydot + x.*y2dot);
elseif m==1
  MassQuadrupole = (-2).*x2dot.*z+(-2).*x.*z2dot+(-4).*xdot.*zdot+(sqrt(-1)*2).*(y2dot.*z+y.*z2dot+2.*ydot.*zdot);
elseif m==-2
  MassQuadrupole = (-2).*x2dot.*z+(-2).*x.*z2dot+(-4).*xdot.*zdot+(sqrt(-1)*2).*(y2dot.*z+y.*z2dot+2.*ydot.*zdot);
end

if m==2
  CurrQuadrupole = (sqrt(-1).*x+y).*((-1).*x3dot.*z+sqrt(-1).*y3dot.*z+(-1).*xdot.* ...
  z2dot+sqrt(-1).*ydot.*z2dot+x.*z3dot+(sqrt(-1)*(-1)).*y.*z3dot+( ...
  x2dot+(sqrt(-1)*(-1)).*y2dot).*zdot);
elseif m==1
  CurrQuadrupole = x.*x3dot.*y+(sqrt(-1)*(-1)).*x3dot.*y.^2+x.*xdot.*y2dot+(sqrt(-1)* ...
  (-1)).*xdot.*y.*y2dot+(-1).*x.^2.*y3dot+sqrt(-1).*x.*y.*y3dot+ ...
  sqrt(-1).*x3dot.*z.^2+y3dot.*z.^2+sqrt(-1).*xdot.*z.*z2dot+ydot.*( ...
  (-1).*x2dot.*(x+(sqrt(-1)*(-1)).*y)+z.*z2dot)+(sqrt(-1)*(-1)).*x.* ...
  z.*z3dot+(-1).*y.*z.*z3dot+((sqrt(-1)*(-1)).*x2dot.*z+(-1).* ...
  y2dot.*z).*zdot;
end

Ulm = sqrt(8.*pi./5).*MassQuadrupole;
Vlm = (-8/3).*((2/5).*pi).^(1/2).*CurrQuadrupole;

out = -1/sqrt(2).*(Ulm - 1i*Vlm);

% add tail and energy contribution
%rOmg = r.*(1+a.*r.^(-3/2)).^(2/3);
%rOmg = r;
%v_phi = rOmg.*Omega;
%phi = DB.phi;
%out = -8.*sqrt(pi/5)* v_phi.^2 .*exp(-2*1i.*phi).*Heff.*Tail.*(rho.^2).*exp(1i.*delta);
%out = out.*Heff.*Tail.*(rho.^2).*exp(1i.*delta);

return


function out_odd = generic_wave_odd(t, r, phi, th,DotProd)
    
phidot  = DB_D1(phi, t, 4);
phi2dot = DB_D1(phidot, t, 4);

rdot   = DB_D1(r, t, 4);
r2dot  = DB_D1(rdot, t, 4);
thdot  = DB_D1(th, t, 4);
th2dot = DB_D1(thdot, t, 4);

% circ term here is (-2).*phidot.^2.*r.^2
out = (1/8).*Conjugate.*DotProd.*exp(1).^((sqrt(-1)*2).*phi).*(5.*pi.^( ...
  -1)).^(1/2).*((-1).*phi2dot.*r.^2+(sqrt(-1)*(-2)).*phidot.^2.* ...
  r.^2+sqrt(-1).*r.*r2dot+(-4).*phidot.*r.*rdot+sqrt(-1).*rdot.^2+( ...
  phi2dot.*r.^2+sqrt(-1).*(2.*phidot.^2.*r.^2+(-1).*r.*r2dot+(sqrt( ...
  -1)*(-4)).*phidot.*r.*rdot+(-1).*rdot.^2+2.*r.^2.*thdot.^2)).*cos( ...
  2.*th)+sqrt(-1).*r.*(r.*th2dot+(sqrt(-1)*4).*phidot.*r.*thdot+4.* ...
  rdot.*thdot).*sin(2.*th));

% fix norm 
out = out.*32.*sqrt(2)*pi/5;

return