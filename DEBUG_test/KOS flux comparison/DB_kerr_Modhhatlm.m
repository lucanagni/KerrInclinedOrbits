function Mod_hhat = kerr_Modhhatlm(props,xin,Omega,E,Heff,jhat,a)

    rr_flag = props.rr_flag; 
    
    % choose the argument in the flm (separately from the argument of
    % the leading-order, Newtonian part, of the flux)
    
    switch rr_flag
            
        case 'DN'
            %---------------------------------------------
            % default input is x=v_phi^2 = (r_Omg.*Omg)^2;
            % Keep it in the rholm's
            %---------------------------------------------
            x = xin;   
            %x = Omega.^(2/3);
        case 'BB'
            %------------------------------------------------------------------
            % Barausse et al. Use v_Omg = Omg.^(1/3) as argument of the rholm's
            % Explicitly imposing Kepler's constraint
            %------------------------------------------------------------------
            x = Omega.^(2/3);   
    end
    
    %-----------------------------------------------------------
    % compute the flm's. Use the factored rhof of Pan et al 2011
    %-----------------------------------------------------------
    %s=kerr_flm_rho(x,a,sqrt(x));
    s = kerr_flm_rho(props,x,a);
    
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
    
    
    function mTlm=ModTail( ell, mm, MOmega)
    
    hatk = mm*MOmega;  
    one  = 1.d0;
    %- ======================================================
    %- Compute the MODULUS of the complex tail factor: |T_lm|
    %- ======================================================
    x    = -2.d0.*hatk;
    x2   = x.^2;
    %{
    prod_sq_num = one;
    prod_sq_den = one;
    for s=1:ell
        s2 = s.^2;
        prod_sq_num = prod_sq_num.*(s2+x2);
        prod_sq_den = prod_sq_den.*s2;
    end
        
    mTlm = abs(sqrt(pi*x./sinh(pi*x).*prod_sq_num./prod_sq_den).*exp(pi*hatk) );
    %}
    %::::::::::::::::::::::::::::::::::::
    % Simpler (equivalent) implementation
    %::::::::::::::::::::::::::::::::::::
    %%{
    prod = one;
    for s=1:ell
        s2   = s.^2;
        prod = prod.*(s2+x2);
    end
    
    mTlm2 = 1/(factorial(ell).^2) * 4*pi*hatk./(one - exp(-4*pi*hatk)).*prod;
    mTlm  = sqrt(mTlm2);
    %}
    return
    
    