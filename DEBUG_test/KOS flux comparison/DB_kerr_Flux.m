function [hatF_resum, F_resum]=DB_kerr_Flux(props,x,Omega,E,Heff,jhat,lmax,a,r_Omg,v_phi)

    % DINFLUX This function computes the Newton.Normalized energy flux according to
    %         the DIN resummation procedure. It is also designed so to add non-QC 
    %         and non-K corrections to  (2,2) partial flux.
    %
    %         USAGE:
    %
    %         [Flm F hatF hatF_resum]=DINFlux(x,Omega,E,Heff,jhat,nu,lmax,r,pr_star,ddotr)
    %
    %         where:
    %
    %         x       :: PN argument
    %         Omega   :: Orbital frequency
    %         E       :: Energy
    %         Heff    :: Effective energy
    %         jhat    :: Newton-Normalized angular momentum
    %         nu      :: symmetric mass ratio
    %         lmax    :: maximum l
    %         r       :: EOB radius
    %         pr_star :: radial momentum
    %         ddotr   :: \ddot{r}
    %
    % Author: Alessandro Nagar
    % W&W   : @IHES, December 2008
    
    m2flux   = props.m2flux;
    l2m2flux = props.l2m2flux;
    rr_flag  = props.rr_flag;
    
    F    = struct('F',[1 8]);
    % checks
    if lmax>8
        error('lmax must not be larger than 8!');
    end
    
    %=========================================
    % Resummed PN-correction to the waveform
    %-----------------------------------------
    % modulus of hhatlm: - optimized
    Mod_hhat = kerr_Modhhatlm(props,x,Omega,E,Heff,jhat,a);
    
    %===============
    % Newtonian flux
    %---------------
    switch rr_flag
        
        case 'BB'                
            % Newtonian prefactor with particular treatment
            % of modes 21 and 44
            FNewt = kerr_FlmNewtBB(x,v_phi,r_Omg,Omega);
            
        otherwise
            FNewt = kerr_FlmNewt(x);
            
    end
    
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
    hatF_resum = 0.; % initializing
    
    if m2flux
        % computing only the m=2 flux
        hatF_resum = Flm.l2m2 + Flm.l3m2 + Flm.l4m2 + Flm.l5m2 ...
                   + Flm.l6m2 + Flm.l7m2 + Flm.l8m2;  
    elseif l2m2flux
        hatF_resum = Flm.l2m2;
    else
        % sum over all multipoles
        for n=2:lmax
            hatF_resum = hatF_resum + F(n).l;
        end
    end
    
    F_resum    = hatF_resum;                % not normalized to Newtonian
    hatF_resum = hatF_resum./FNewt.l2m2;     % normalized to Newtonian
    
    
    return
    
    
    
    
    