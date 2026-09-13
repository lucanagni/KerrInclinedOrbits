function F = DB_flux_KOS(obj,X,P,dHdp,q,chi1,chi2)
%=============================================================================
% Attemot to implement KOS RR (Newt) to compare non-geodesic dynamics
% Not working
%=============================================================================
    mu = q/(1+q)^2;
    a = norm(chi1);
    [r,phi,th,pr,pph] = DB_coords_cart2spherical(X(1),X(2),X(3),P(1),P(2),P(3));
    u = 1./r;

    r_Omg  = r.*(1+a.*r.^(-3/2)).^(2/3);
    Omg    = -dHdp(1).*sin(phi).*u + dHdp(2).*cos(phi).*u;
    v_phi  = r_Omg.*Omg;
    jNewt  = r_Omg.*v_phi;   % "Formal" Newtonian angular momentum
    jhat   = pph./jNewt;     %  Newton-normalized angular momentum
    x      = v_phi.^2;

    [H,Horb,dHeff,dHorb]      = DB_Hamiltonian(obj,X,P,chi1,chi2); 
    [A, Bp,Bnp,~,~,dAC] = DB_metric_Kerr(X,chi1);
    dA                  = dAC.dx(1).*cos(phi)+dAC.dx(2).*sin(phi);
    dHdr                = dHeff.dx(1).*cos(phi) + dHeff.dx(2).*sin(phi) + ...
                          dHeff.dp(1).*pph.*sin(phi)./r.^2 - dHeff.dp(2).*pph.*cos(phi)./r.^2;
    dH.dr               = dHdr;
    dH.Horb = Horb;

    B                   = 1./(Bp+Bnp); % crf. Balmelli-Damour 2015
    squAB               = sqrt(A./B);
    prstar              = squAB.*pr;

    [hatF, ~] = KOS_Flux(x,Omg,1,H,jhat,8,a);
    Fphi = DB_Fphi_Iterate(a, r, prstar, pph, Omg, hatF, mu, 2,A,B,dA,dH);
    [~,~,~,Fx,Fy,Fz] = DB_coords_spherical2cart(r,phi,th,0,Fphi,0); % convert Fphi to cartesian

    F = [Fx;Fy;Fz];
return

function [hatF_resum, F_resum] = KOS_Flux(x,Omega,E,Heff,jhat,lmax,a) % needed in DB_flux_KOS

    F    = struct('F',[1 8]);
    % checks
    if lmax>8
        error('lmax must not be larger than 8!');
    end

    %=========================================
    % Resummed PN-correction to the waveform
    %-----------------------------------------
    % modulus of hhatlm: - optimized
    Mod_hhat = KOS_Modhhatlm(x,Omega,E,Heff,jhat,a);

    FNewt = KOS_FlmNewt(x);

    %===============
    % The total flux
    %---------------

    Flm.l2m2 = Mod_hhat.l2m2.^2.*FNewt.l2m2; 
    Flm.l2m1 = Mod_hhat.l2m1.^2.*FNewt.l2m1;
    Flm.l3m3 = Mod_hhat.l3m3.^2.*FNewt.l3m3;
    Flm.l3m2 = Mod_hhat.l3m2.^2.*FNewt.l3m2;
    Flm.l3m1 = Mod_hhat.l3m1.^2.*FNewt.l3m1;
    Flm.l4m4 = Mod_hhat.l4m4.^2.*FNewt.l4m4;
    Flm.l4m3 = Mod_hhat.l4m3.^2.*FNewt.l4m3;
    Flm.l4m2 = Mod_hhat.l4m2.^2.*FNewt.l4m2;
    Flm.l4m1 = Mod_hhat.l4m1.^2.*FNewt.l4m1;
    Flm.l5m5 = Mod_hhat.l5m5.^2.*FNewt.l5m5;
    Flm.l5m4 = Mod_hhat.l5m4.^2.*FNewt.l5m4;
    Flm.l5m3 = Mod_hhat.l5m3.^2.*FNewt.l5m3;
    Flm.l5m2 = Mod_hhat.l5m2.^2.*FNewt.l5m2;
    Flm.l5m1 = Mod_hhat.l5m1.^2.*FNewt.l5m1;
    Flm.l6m6 = Mod_hhat.l6m6.^2.*FNewt.l6m6;
    Flm.l6m5 = Mod_hhat.l6m5.^2.*FNewt.l6m5;
    Flm.l6m4 = Mod_hhat.l6m4.^2.*FNewt.l6m4;
    Flm.l6m3 = Mod_hhat.l6m3.^2.*FNewt.l6m3;
    Flm.l6m2 = Mod_hhat.l6m2.^2.*FNewt.l6m2;
    Flm.l6m1 = Mod_hhat.l6m1.^2.*FNewt.l6m1;
    Flm.l7m7 = Mod_hhat.l7m7.^2.*FNewt.l7m7;
    Flm.l7m6 = Mod_hhat.l7m6.^2.*FNewt.l7m6;
    Flm.l7m5 = Mod_hhat.l7m5.^2.*FNewt.l7m5;
    Flm.l7m4 = Mod_hhat.l7m4.^2.*FNewt.l7m4;
    Flm.l7m3 = Mod_hhat.l7m3.^2.*FNewt.l7m3;
    Flm.l7m2 = Mod_hhat.l7m2.^2.*FNewt.l7m2;
    Flm.l7m1 = Mod_hhat.l7m1.^2.*FNewt.l7m1;
    Flm.l8m8 = Mod_hhat.l8m8.^2.*FNewt.l8m8;
    Flm.l8m7 = Mod_hhat.l8m7.^2.*FNewt.l8m7;
    Flm.l8m6 = Mod_hhat.l8m6.^2.*FNewt.l8m6;
    Flm.l8m5 = Mod_hhat.l8m5.^2.*FNewt.l8m5;
    Flm.l8m4 = Mod_hhat.l8m4.^2.*FNewt.l8m4;
    Flm.l8m3 = Mod_hhat.l8m3.^2.*FNewt.l8m3;
    Flm.l8m2 = Mod_hhat.l8m2.^2.*FNewt.l8m2;
    Flm.l8m1 = Mod_hhat.l8m1.^2.*FNewt.l8m1;
    %-
    %=================================================
    % partial, Newton-normalized, multipoles up to l=8
    %-------------------------------------------------
    F(2).l =  Flm.l2m2 + Flm.l2m1;
    F(3).l =  Flm.l3m3 + Flm.l3m2 + Flm.l3m1;
    F(4).l =  Flm.l4m4 + Flm.l4m3 + Flm.l4m2 + Flm.l4m1;
    F(5).l =  Flm.l5m5 + Flm.l5m4 + Flm.l5m3 + Flm.l5m2 + Flm.l5m1;
    F(6).l =  Flm.l6m6 + Flm.l6m5 + Flm.l6m4 + Flm.l6m3 + Flm.l6m2 + Flm.l6m1;
    F(7).l =  Flm.l7m7 + Flm.l7m6 + Flm.l7m5 + Flm.l7m4 + Flm.l7m3 + Flm.l7m2 + Flm.l7m1;
    F(8).l =  Flm.l8m8 + Flm.l8m7 + Flm.l8m6 + Flm.l8m5 + Flm.l8m4 + Flm.l8m3 + Flm.l8m2 + Flm.l8m1;
    %

    %==================================================
    % Summing over lmax multipoles and Newton-normalize
    % to the 22 multipole. This is the f(v) function
    %--------------------------------------------------
    % sum over all multipoles
    hatF_resum = 0.; % initializing
    for n=2:lmax
        hatF_resum = hatF_resum + F(n).l;
    end

    F_resum    = hatF_resum;                % not normalized to Newtonian
    hatF_resum = hatF_resum./FNewt.l2m2;     % normalized to Newtonian
return

function s = KOS_FlmNewt(x)
    % DINFlmNewt.m 
    %            This function computes the GW energy flux from a circularized
    %            binary in its multipolar decoposition.
    %            The output is given by Flm/nu^2, where nu=m1m2/(m1+m2)^2 is
    %            the symmetric mass ratio. 
    %            USAGE:
    %            
    %            F=DINFlmNewt(x,nu)
    %
    %            where F is a structure of the form s.lm. Values up to l=8 are
    %            explicitly included.
    %
    %            Author: Alessandro Nagar
    %            When:   @IHES, November 2008 & November 2009
    %
    
    % shorthands: variable
    x5  = x.^5;
    x6  = x.*x5;
    x7  = x.*x6;
    x8  = x.*x7;
    x9  = x.*x8;
    x10 = x.*x9;
    x11 = x.*x10;
    x12 = x.*x11;
    
    
    % Newtonian partial fluxes
    s.l2m2    = x5* 32.d0/5.d0;                       %l2m2
    s.l2m1    = x6* 8.d0/45.d0;                       %l2m1
    s.l3m3    = x6* 243.d0/28.d0;                     %l3m3
    s.l3m2    = x7* 32.d0/63.d0;                      %l3m2
    s.l3m1    = x6/ 1260.d0;                          %l3m1
    s.l4m4    = x7* 8192.d0/567.d0;                   %l4m4
    s.l4m3    = x8* 729.d0/700.d0;                    %l4m3
    s.l4m2    = x7* 32.d0/3969.d0;                    %l4m2
    s.l4m1    = x8/ 44100.d0;                         %l4m1
    s.l5m5    = x8* 1953125.d0/76032.d0;              %l5m5
    s.l5m4    = x9* 131072.d0/66825.d0;               %l5m4
    s.l5m3    = x8* 2187.d0/70400.d0;                 %l5m3
    s.l5m2    = x9* 256.d0/400950.d0;                 %l5m2
    s.l5m1    = x8/ 19958400.d0;                      %l5m1
    s.l6m6    = x9* 839808.d0/17875.d0;               %l6m6
    s.l6m5    = x10* 48828125.d0/13621608.d0;         %l6m5
    s.l6m4    = x9* 4194304.d0/47779875.d0;           %l6m4
    s.l6m3    = x10* 59049.d0/15415400.d0;            %l6m3
    s.l6m2    = x9* 128.d0/28667925.d0;               %l6m2
    s.l6m1    = x10/ 1123782660.d0;                   %l6m1
    s.l7m7    = x10* 96889010407.d0/1111968000.d0;    %l7m7
    s.l7m6    = x11* 5668704.d0/875875.d0;            %l7m6
    s.l7m5    = x10* 1220703125.d0/5.666588928d9;     %l7m5
    s.l7m4    = x11* 4194304.d0/3.07432125d8;         %l7m4
    s.l7m3    = x10* 1594323.d0/3.2064032d10;         %l7m3
    s.l7m2    = x11* 32.d0/1.35270135d8;              %l7m2
    s.l7m1    = x10/ 9.3498717312d11;                 %l7m1
    s.l8m8    = x11* 274877906944.d0/1688511825.d0;   %l8m8
    s.l8m7    = x12* 4747561509943.d0/4.083146496d11; %l8m7
    s.l8m6    = x11* 51018336.d0/1.04229125d8;        %l8m6
    s.l8m5    = x12* 30517578125.d0/8.00296713216d11; %l8m5
    s.l8m4    = x11* 4194304.d0/1.5679038375d10;      %l8m4
    s.l8m3    = x12* 177147.d0/3.96428032d10;         %l8m3
    s.l8m2    = x11* 32.d0/3.4493884425d10;           %l8m2
    s.l8m1    = x12/ 8.174459284992d13;               %l8m0
    
return

function Mod_hhat = KOS_Modhhatlm(xin,Omega,E,Heff,jhat,a) %needed in KOS_Flux

    x = xin;   
    %x = Omega.^(2/3);

    %-----------------------------------------------------------
    % compute the flm's. Use the factored rhof of Pan et al 2011
    %-----------------------------------------------------------
    %s=kerr_flm_rho(x,a,sqrt(x));
    s = DB_KOS_rholm_iResum(x, a);
    
    % Rename structures (awful but it works!)
    f22 = s.ell(2).emm(2).flm;
    f21 = s.ell(2).emm(1).flm;
    f33 = s.ell(3).emm(3).flm;
    f32 = s.ell(3).emm(2).flm;
    f31 = s.ell(3).emm(1).flm;
    f44 = s.ell(4).emm(4).flm;
    f43 = s.ell(4).emm(3).flm;
    f42 = s.ell(4).emm(2).flm;
    f41 = s.ell(4).emm(1).flm;
    f55 = s.ell(5).emm(5).flm;
    f54 = s.ell(5).emm(4).flm;
    f53 = s.ell(5).emm(3).flm;
    f52 = s.ell(5).emm(2).flm;
    f51 = s.ell(5).emm(1).flm;
    f66 = s.ell(6).emm(6).flm;
    f65 = s.ell(6).emm(5).flm;
    f64 = s.ell(6).emm(4).flm;
    f63 = s.ell(6).emm(3).flm;
    f62 = s.ell(6).emm(2).flm;
    f61 = s.ell(6).emm(1).flm;
    f77 = s.ell(7).emm(7).flm;
    f76 = s.ell(7).emm(6).flm;
    f75 = s.ell(7).emm(5).flm;
    f74 = s.ell(7).emm(4).flm;
    f73 = s.ell(7).emm(3).flm;
    f72 = s.ell(7).emm(2).flm;
    f71 = s.ell(7).emm(1).flm;
    f88 = s.ell(8).emm(8).flm;
    f87 = s.ell(8).emm(7).flm;
    f86 = s.ell(8).emm(6).flm;
    f85 = s.ell(8).emm(5).flm;
    f84 = s.ell(8).emm(4).flm;
    f83 = s.ell(8).emm(3).flm;
    f82 = s.ell(8).emm(2).flm;
    f81 = s.ell(8).emm(1).flm;
    
    %========================================================
    % The PN correction to the waveform.
    % Note that jhat must be put to one when one has factored
    % the full angular momentum J Omega in the Newtonian part 
    % Beware this is the absolute value of the \hat{h}_lm
    %========================================================
    EOmg = E.*Omega;
        
    % --------------
    % using the rhof
    % --------------
    Mod_hhat.l2m2 = ModTail(2,2,EOmg).*Heff.*f22;
    Mod_hhat.l2m1 = ModTail(2,1,EOmg).*jhat.*f21;
    Mod_hhat.l3m3 = ModTail(3,3,EOmg).*Heff.*f33;
    Mod_hhat.l3m2 = ModTail(3,2,EOmg).*jhat.*f32;
    Mod_hhat.l3m1 = ModTail(3,1,EOmg).*Heff.*f31;
    Mod_hhat.l4m4 = ModTail(4,4,EOmg).*Heff.*f44;
    Mod_hhat.l4m3 = ModTail(4,3,EOmg).*jhat.*f43;
    Mod_hhat.l4m2 = ModTail(4,2,EOmg).*Heff.*f42;
    Mod_hhat.l4m1 = ModTail(4,1,EOmg).*jhat.*f41;
    
    %== ell=5
    Mod_hhat.l5m5 = Heff.*ModTail(5,5,EOmg).*f55;
    Mod_hhat.l5m4 = jhat.*ModTail(5,4,EOmg).*f54;
    Mod_hhat.l5m3 = Heff.*ModTail(5,3,EOmg).*f53;
    Mod_hhat.l5m2 = jhat.*ModTail(5,2,EOmg).*f52;
    Mod_hhat.l5m1 = Heff.*ModTail(5,1,EOmg).*f51;
    %== ell=6
    Mod_hhat.l6m6 = Heff.*ModTail(6,6,EOmg).*f66;
    Mod_hhat.l6m5 = jhat.*ModTail(6,5,EOmg).*f65;
    Mod_hhat.l6m4 = Heff.*ModTail(6,4,EOmg).*f64;
    Mod_hhat.l6m3 = jhat.*ModTail(6,3,EOmg).*f63;
    Mod_hhat.l6m2 = Heff.*ModTail(6,2,EOmg).*f62;
    Mod_hhat.l6m1 = jhat.*ModTail(6,1,EOmg).*f61;
    %== ell=7
    Mod_hhat.l7m7 = Heff.*ModTail(7,7,EOmg).*f77;
    Mod_hhat.l7m6 = jhat.*ModTail(7,6,EOmg).*f76;
    Mod_hhat.l7m5 = Heff.*ModTail(7,5,EOmg).*f75;
    Mod_hhat.l7m4 = jhat.*ModTail(7,4,EOmg).*f74;
    Mod_hhat.l7m3 = Heff.*ModTail(7,3,EOmg).*f73;
    Mod_hhat.l7m2 = jhat.*ModTail(7,2,EOmg).*f72;
    Mod_hhat.l7m1 = Heff.*ModTail(7,1,EOmg).*f71;
    %== ell=8
    Mod_hhat.l8m8 = Heff.*ModTail(8,8,EOmg).*f88;
    Mod_hhat.l8m7 = jhat.*ModTail(8,7,EOmg).*f87;
    Mod_hhat.l8m6 = Heff.*ModTail(8,6,EOmg).*f86;
    Mod_hhat.l8m5 = jhat.*ModTail(8,5,EOmg).*f85;
    Mod_hhat.l8m4 = Heff.*ModTail(8,4,EOmg).*f84;
    Mod_hhat.l8m3 = jhat.*ModTail(8,3,EOmg).*f83;
    Mod_hhat.l8m2 = Heff.*ModTail(8,2,EOmg).*f82;
    Mod_hhat.l8m1 = jhat.*ModTail(8,1,EOmg).*f81;
    
return
    
function mTlm = ModTail( ell, mm, MOmega)
    
    hatk = mm*MOmega;  
    one  = 1.d0;
    %- ======================================================
    %- Compute the MODULUS of the complex tail factor: |T_lm|
    %- ======================================================
    x    = -2.d0.*hatk;
    x2   = x.^2;

    prod = one;
    for s=1:ell
        s2   = s.^2;
        prod = prod.*(s2+x2);
    end
    
    mTlm2 = 1/(factorial(ell).^2) * 4*pi*hatk./(one - exp(-4*pi*hatk)).*prod;
    mTlm  = sqrt(mTlm2);
return
    
