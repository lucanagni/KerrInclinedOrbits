function [Fphi, Fr, r2dot, r3dot, Omgdot, Omg2dot, prstar2dot] = DB_Fphi_Iterate(a, r, prstar, pph, Omg, hatF, mu, iter,A,B,dA,dH)

% Iterative procedure for the computation of the angular radiation reaction
% with full Newtonian prefactor. Time derivatives of the fluxes are
% neglected.

%[~,~,~,dH]    = DB_Hamiltonian(obj,X,P,chi1,chi2);
%[A, Bp,Bnp,~,~,dA] = DB_metric_Kerr(x,chi1);

%B = 1./(Bp+Bnp);
A2           = A.*A;
B2           = B.*B;

squAB        = sqrt(A./B);
%H            = dH.Heff;
Horb         = dH.Horb;
dHdr         = dH.dr;
onebyHorb    = 1./Horb;
rdot         = squAB.*prstar./Horb;
rdot2        = rdot.*rdot;

u            = 1./r;
u2           = u.*u;
u3           = u2.*u;
u4           = u3.*u;

r2           = r.*r;
r3           = r2.*r;

a2           = a*a;
a4           = a2*a2;

rc2 = r2 + a2.*(1+2*u);
rc = sqrt(rc2);

uc           = 1./rc;
uc2          = uc.*uc;
uc3          = uc2.*uc;
uc4          = uc3.*uc;

pph2         = pph.*pph;
%pph2dot      = 0;
%Radial flux; not necessary to iterate here, but maybe someday

[~, ~, Fr]   = KOS_Fr_BD(r, 0, prstar, pph, 0);
Fr           = mu*Fr;
Frdot        = 0;

drc_dr       = r.*uc.*(1 - a2.*u3);
d2rc_dr2     = uc.*(1 + 2*a2.*u3) - drc_dr.^2.*uc;

Gs           =  2*u.*uc2;
dGs_dr       = -2*(a2+3*r2)./((r3+a2*(2+r)).^2);
d2Gs_dr2     = 4.*(a4+3.*a2.*((-2)+r).*r+6.*r.^4).*(r3+a2.*(2+r)).^(-3);

Gtilde       = Gs*a;
dGtilde_dr   = dGs_dr*a;
d2Gtilde_dr2 = d2Gs_dr2*a;

dB       = -2*B2.*u2.*(1-a2*u);
d2B      = 2*dB.*dB./B - 2*a2*u4.*B2 - 2*u.*dB;

d2A          = 4.*(2+r).^(-3).*rc.^(-4).*((-6).*drc_dr.^2.*r.*(2+r).^2+2.*(2+r).*( ...
               4.*drc_dr+d2rc_dr2.*r.*(2+r)).*rc+4.*rc.^2+(-1).*rc.^4);

factor       = (dA.*B - A.*dB)./B2;
dsquAB_dr    = 0.5./squAB.*factor;
d2squAB_d2r  = -1/4.*squAB.*(dA.*B - A.*dB).^2./(A2.*B2)...
             + 0.5./squAB.*(d2A.*B2 - A.*B.*d2B - 2*dA.*dB.*B + 2.*A.*dB.^2)./(B.*B2);
         
dAbyrc2      = (dA.*uc2 - 2*A.*uc3.*drc_dr);      
d2Abyrc2     = d2A.*uc2 - 4*dA.*uc3.*drc_dr  + 6*A.*uc4.*(drc_dr).^2 - 2*A.*uc3.*d2rc_dr2;

Omgdot_0     = onebyHorb.*pph.*rdot.*dAbyrc2...
                         + dGtilde_dr.*rdot;
prstardot    = -squAB.*(-Fr + dHdr);
Fphi         = -32.d0/5.d0*mu*r.^4.*Omg.^5.*hatF;


% Kerr Kepler law (check Straumann's book, Eq.(7.256)
r_Omg = r.*(1+a.*r.^(-3/2)).^(2/3);

Fphi   = Fphi.*r_Omg.^4./(r.^4);

% Iterations
for k = 0 : iter
    
    Hdot       = rdot.*Fr + Omg.*Fphi;
    HSdot      = dGtilde_dr.*rdot.*pph + Gtilde.*Fphi;
    Horbdot    = Hdot - HSdot;
   
    Omgdot     = Omgdot_0 + onebyHorb.*A.*uc2.*(Fphi-pph.*onebyHorb.*Horbdot);
    
    r2dot      = dsquAB_dr./squAB.*rdot2 + squAB.*(...
                 onebyHorb.*(prstardot - prstar.*onebyHorb.*Horbdot));
             
    pph2dot    = Fphi.*(4.*rdot./r + 5.*Omgdot./Omg);
             
    H2dot      = r2dot.*Fr + Omgdot.*Fphi + rdot.*Frdot + Omg.*pph2dot;
    HS2dot     = d2Gtilde_dr2.*rdot.^2.*pph + dGtilde_dr.*r2dot.*pph + 2*dGtilde_dr.*rdot.*Fphi + Gtilde.*pph2dot;
    Horb2dot   = H2dot - HS2dot;
             
    prstar2dot = dsquAB_dr./squAB.*rdot.*prstardot + squAB.*(Frdot ...
                 -0.5*onebyHorb.*(rdot.*(d2A + d2Abyrc2.*pph2) - dA.*Horbdot.*onebyHorb ...
                 + dAbyrc2.*(pph.*(2*Fphi - pph.*Horbdot.*onebyHorb)))...
                 -(dGtilde_dr.*Fphi + pph.*(d2Gtilde_dr2.*rdot)));
             
    Omg2dot    = onebyHorb.*(d2Abyrc2.*pph.*rdot2 + dAbyrc2.*pph.*r2dot + 2*dAbyrc2.*rdot.*(Fphi-pph.*Horbdot.*onebyHorb)...
                                     + A.*uc2.*(pph2dot - 2*Fphi.*Horbdot.*onebyHorb ...
                                     + pph.*(2*(Horbdot.*onebyHorb).^2 - Horb2dot.*onebyHorb)))...
               + (d2Gtilde_dr2.*rdot2 + dGtilde_dr.*(r2dot));
           
    r3dot = rdot./squAB.*(dsquAB_dr.*(3*r2dot - 2*dsquAB_dr./squAB.*rdot2) + d2squAB_d2r.*rdot2)...
                  + onebyHorb.*squAB.*(prstar2dot - 2*prstardot.*Horbdot.*onebyHorb ...
                  + prstar.*(2*(Horbdot.*onebyHorb).^2 - Horb2dot.*onebyHorb));
%               r3dot = 0;
              
%     r3dot = -B./A.*dsquAB_dr.^2.*rdot2 + d2squAB_d2r./squAB.*rdot.^3 + 2./squAB.*dsquAB_dr.*rdot.*r2dot + ...
%             dsquAB_dr.*(prstardot.*onebyHorb - 2.*prstar.*Horbdot.*onebyHorb.^2) +...
%             squAB.*onebyHorb.*(prstar2dot - 2*prstardot.*Horbdot.*onebyHorb - prstar.*H2dot.*onebyHorb + 2*prstar.*Hdot.^2.*onebyHorb.^2);
    
    Fphi_NPref = KOS_Fphi_NewtPref(r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot);
    Fphi       = -32.d0/5.d0*mu*Fphi_NPref.*hatF; % here rOmg is missing
    Fphi       = Fphi.*r_Omg.^4./(r.^4);
    
end


end

function [Fr, Fr_s, FrC] = KOS_Fr_BD(r, pr, prstar, pph, nu)

    % Computes the radial radiation reaction using the complete formula for
    % general orbits to 2PN, from Bini-Damour (2012), as a function of pr and prstar; as well as its circular
    % approximation, obtained by setting pr = dpr/dt = 0.
    
    % Shorthands
    
    u = 1./r;
    u2 = u.*u;
    u3 = u2.*u;
    u4 = u2.*u2;
    u5 = u3.*u2;
    %u6 = u3.*u3;
    %u7 = u4.*u3;
    
    %pr2 = pr.*pr;
    %pr4 = pr2.*pr2;
    %pr6 = pr4.*pr2;
    
    prs2 = prstar.*prstar;
    prs4 = prs2.*prs2;
    
    pph2 = pph.*pph;
    pph4 = pph2.*pph2;
    pph6 = pph2.*pph4;
    %pph8 = pph6.*pph2;
    
    nu2 = nu*nu;
    %nu3 = nu2*nu;
    
    % 0PN
    
    f0pn   = -8/15 + 56/5*pph2.*u;
    
    f0pn_s = f0pn;
    
    % 1pn
    %{
    f1pn   = pr2*(-1228/105 + 556/105*nu) + u*(-624/35 + 16/21*nu) + ...
             pr2.*pph2.*u*(-124/105 - 436/105*nu) + pph2.*u2*(-496/7 - 1268/105*nu) + ...
             pph4.*u3*(1252/105 - 2588/105*nu);  
    %}
    f1pn_s = prs2*(-1228/105 + 556/105*nu) + u*(-1984/105 + 16/21*nu) + ...
             prs2.*pph2.*u*(-124/105 - 436/105*nu) + pph2.*u2*(-1696/35 - 1268/105*nu) + ...
             pph4.*u3*(1252/105 - 2588/105*nu);  
    
    % 2pn
    %{
    f2pn = pr4*(323/315 + 1061/315*nu - 1273/315*nu2) + ...
           pph6.*u5*(-3229/315 - 718/63*nu + 3277/105*nu2) + ...
           pph4.*u4*(-6103/45 + 16418/105*nu + 25217*nu2) + ...
           u2*(33338/567 + 9526/105*nu - 3548/315*nu2) + ...
           pr2.*pph2.*u2*(-8126/105 + 3628/63*nu - 8804/315*nu2) + ...
           pr2.*u*(1222/9 + 57926/945*nu - 218/189*nu2) + ...
           pr4.*pph2.*u*(-461/315 - 983/315*nu + 131/63*nu2) + ...
           pph2.*u3*(218/63 + 66032/315*nu - 1752/35*nu2) + ...
           pr2.*u3.*pph4*(-628/105 - 1052/105*nu + 194/7*nu2);
    %}
    f2pn_s = prs4*(323/315 + 1061/315*nu - 1273/315*nu2) + ...
             pph6.*u5*(-3229/315 - 718/63*nu + 3277/105*nu2) + ...
             pph4.*u4*(-35209/315 + 1606/15*nu + 25217/315*nu2) + ...
             u2*(59554/2835 + 9686/105*nu - 3548/315*nu2) + ...
             prs2.*pph2.*u2*(-1774/21 + 10292/315*nu - 8804/315*nu2) + ...
             prs2.*u*(20666/315 + 17590/189*nu - 218/189*nu2) + ...
             prs4.*pph2.*u*(-461/315 - 983/315*nu + 131/63*nu2) + ...
             pph2.*u3*(-29438/315 + 58424/315*nu - 1752/35*nu2) + ...
             prs2.*u3.*pph4*(-628/105 - 1052/105*nu + 194/7*nu2);
    
    %{
    % 3pn
    
    % f3pn = (1/6930).*(5809+7543.*nu+(-7946).*nu2+23539.*nu3).*pr6+( ...
    %     -1/6930).*((-120973)+232359.*nu+(-442512).*nu2+115839.*nu3).* ...
    %     pph8.*u7+(1/3118500).*(749760300.*pph6+(-324411750).*nu.* ...
    %     pph6+(-811430550).*nu2.*pph6+(-889123050).*nu3.*pph6).* ...
    %     u6+((1/198).*((-14883)+119740.*nu+(-267138).*nu2+64684.* ...
    %     nu3).*pph2.*pr2+(-2/51975).*(753498+1797370.*nu+1357515.* ...
    %     nu2+186440.*nu3).*pr4+(-1/6930).*((-626951)+4206571.*nu+( ...
    %     -9387844).*nu2+2270073.*nu3).*pph2.*pr6).*u+u2 ...
    %     .*((1/34650).*(7494231+6214655.*nu+225330.*nu2+1570705.*nu3).* ...
    %     pph2.*pr4+(1/3118500).*pr2.*((-2).*(528415558+375572700.* ...
    %     nu2+39513700.*nu3+nu.*(733823225+(-6392925).*pi.^2))+( ...
    %     -40677120).*log(r)))+u3.*((-1/6930).*((-770588)+4479136.*nu+( ...
    %     -9578599).*nu2+2449281.*nu3).*pph4.*pr4+(1/3118500).*(4.*( ...
    %     185963428+(-140175).*nu2+(-36225300).*nu3+25.*nu.*(2398799+ ...
    %     340956.*pi.^2))+15750.*((-14883)+119740.*nu+(-267138).*nu2+ ...
    %     64684.*nu3).*pph4+(-54236160).*log(r))+(1/3118500).*pr2.*( ...
    %     5008242222.*pph2+(-754415550).*nu2.*pph2+263985600.*nu3.* ...
    %     pph2+(-3300).*nu.*((-510112)+23247.*pi.^2).*pph2+(-122031360) ...
    %     .*pph2.*log(r)))+u4.*((1/13860).*(1354169+8798367.*nu+( ...
    %     -5520051).*nu2+858666.*nu3).*pph4.*pr2+(1/3118500).*((-40) ...
    %     .*(20246351+(-21424500).*nu2+2653080.*nu3+(-5).*nu.*(( ...
    %     -42402709)+767151.*pi.^2)).*pph2+254232000.*pph2.*log(r)))+ ...
    %     u5.*((-1/4620).*((-294697)+456833.*nu+(-6287).*nu2+472826.* ...
    %     nu3).*pph6.*pr2+(1/3118500).*(1674944820.*pph4+ ...
    %     1917472050.*nu.*pph4+(-3462805350).*nu2.*pph4+1154248200.* ...
    %     nu3.*pph4+(-152539200).*pph4.*log(r)));
    %}
             
    % Totes
    
    %a1 = f1pn./f0pn;
    %a2 = f2pn./f0pn;
    %a3 = f3pn./f0pn;
    
    % Fr = u4.*pr.*f0pn./(1-a1+a1.*a1-a2);%-a1.^3+2*a1.*a2-a3);
    
    %Fr   = u4.*pr.*(f0pn + f1pn + f2pn + 0*f3pn);
    
    % Hijacking Fr to get the raw version, as function of prstar
    
    Fr = u4.*prstar.*(f0pn_s + f1pn_s + f2pn_s);
    
    % Resummed version, with prstar
    
    a1 = f1pn_s./f0pn_s;
    a2 = f2pn_s./f0pn_s;
    
    Fr_s = u4.*prstar.*f0pn_s./(1 - a1 + (a1.*a1 - a2));
    
    %Fr_s = u4.*prstar.*f0pn_s.*(1 + a1 - a2./a1)./(1 - a2./a1);
    
    % Circular part, very useless, but whatevs
    
    frc = 1 + u.*(-(1133/280 + 118/35*nu) + (-175549/15120 + 1707/80*nu + 1311/280*nu*nu)*u);
    
    FrC = 32/3*nu*u.*pr.*frc.*u.*u.*u;
    
end

function Fphi_NPref = KOS_Fphi_NewtPref(r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot)

    % Complete Newtonian prefactor for the angular momentum flux/angular
    % radiation reaction.
    
    %r3dot = 0;
    
    Omg2 = Omg.*Omg;
    Omg3 = Omg2.*Omg;
    Omg4 = Omg3.*Omg;
    
    %Omgdot2 = Omgdot.*Omgdot;
    
    r2 = r.*r;
    r3 = r.*r2;
    r4 = r3.*r;
    
    rdot2 = rdot.*rdot;
    rdot3 = rdot2.*rdot;
    rdot4 = rdot3.*rdot;
    
    r2dot2 = r2dot.*r2dot;
    
    Fphi_NPref = (1/8).*(2.*Omg.*(4.*Omg4+(-1).*Omg.*Omg2dot+3.*Omgdot.^2).*r4+ ...
        6.*Omgdot.*r.*rdot3+6.*Omg.*rdot4+r3.*((-16).*Omg3.*r2dot+ ...
        Omg2dot.*r2dot+(-1).*Omgdot.*r3dot+24.*Omg2.*Omgdot.*rdot)+ ...
        r2.*(6.*Omg.*r2dot2+3.*Omgdot.*r2dot.*rdot+rdot.*((-4).*Omg.* ...
        r3dot+32.*Omg3.*rdot+Omg2dot.*rdot)));
    
    % Fphi_NPref = Fphi_NPref./(r.^4.*Omg.^5);
    
    end