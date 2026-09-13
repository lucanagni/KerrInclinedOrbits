function DB = wave_test

X = load("h_Yl2m2_x10.0000.dat");
inputDB.r0   = 6;
inputDB.th0  = pi/3;
inputDB.chi1 = [0;0;0.2];

addpath ../

DB = DB_class(inputDB);

t   = DB.t;
r   = DB.r;
th  = DB.th;
phi = DB.phi;

h22new  = generic_wave(t, r, phi, th, 0);
h22circ = generic_wave(t, r, phi, pi/2+0.*r, 1);

figure('position', [476 139 730 727])

subplot(2,2,1)
plot(t, r)

subplot(2,2,2)
plot(t, h22circ)
legend('circ')

subplot(2,2,3) 
plot(t, h22new, '-')
hold on
plot(t, abs(h22new))

subplot(2,2,4)
tt = X(:,1);
h = X(:,2)+1j*X(:,3);
plot(tt, real(h))
hold on
plot(tt, imag(h))
plot(tt, abs(h))

for i=1:3
    subplot(2,2,i)
    xlim([0, t(end)-50])
end

rmpath ../

return 

function out = generic_wave(t, r, phi, th, kill_phi2dot)

phidot  = DB_D1(phi, t, 4);
if kill_phi2dot
    phi2dot = r.*0;
else
    phi2dot = DB_D1(phidot, t, 4);
end

rdot   = DB_D1(r, t, 4);
r2dot  = DB_D1(rdot, t, 4);
thdot  = DB_D1(th, t, 4);
th2dot = DB_D1(thdot, t, 4);

out = (1/2).*exp((1i*(-2)).*phi).*((15/2).*pi.^(-1)).^(1/2).* ...
  (r.^2.*thdot.^2.*cos(th).^2+r.*(4.*rdot.*thdot+r.*(th2dot+(sqrt( ...
  -1)*(-4)).*phidot.*thdot)).*cos(th).*sin(th)+((1i*(-1)).* ...
  phi2dot.*r.^2+(-2).*phidot.^2.*r.^2+r.*r2dot+(1i*(-4)).* ...
  phidot.*r.*rdot+rdot.^2+(-1).*r.^2.*thdot.^2).*sin(th).^2);
return


