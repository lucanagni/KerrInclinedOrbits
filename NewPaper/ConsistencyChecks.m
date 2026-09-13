function ConsistencyChecks(l,M)
    % =============================================================================================================================
    % Check consistency of rotation procedure and computation of the rholm's
    % 1) Check that the even (odd) parity equatorial multipoles rotate into the even (odd) parity component of the inclined ones
    % This works! note that this is consistent with Mathematica notebook: 'backward' rotation applies d(-iota) with the way angles are computed 
    % =============================================================================================================================

    %{
    load('~/waveforms/S/plunge/th0_90/more_timesteps/wf.mat')
    seq = s;
    load('~/waveforms/S/plunge/th0_45/wf.mat')
    sin = s;
    p = struct;

    p = DB_mode_rotate_timedep(seq,sin.dyn,2,2,p,'direction','backward');
    %}

    mystruct = struct;
    inputDB.th0 = pi/2;
    inputDB.ICs = 'spherical';
    inputDB.chi1(3) = 0;
    inputDB.geodesics = 1;
    inputDB.Tmax = 1000;
    inputDB.r0 = 7;
    inputDB.verbose = 0;

    DBeq = DB_class(inputDB);
    
    inputDB.th0 = pi/6;
    DBin = DB_class(inputDB);

    Ueq = struct;
    Uin = struct;
    Veq = struct;
    Vin = struct;

    Ueq.dyn = DBeq;
    Veq.dyn = DBeq;
    Uin.dyn = DBin;
    Vin.dyn = DBin;

    flags.r = 1;


    for m=-l:l
        n = abs(m);
        [Ulm_eq,Vlm_eq] = DB_Multipoles(DBeq,l,m,flags);
        [Ulm_in,Vlm_in] = DB_Multipoles(DBin,l,m,flags);

        if m<0
            Ueq.ell(l).emminus(n+1).hlm = Ulm_eq;
            Veq.ell(l).emminus(n+1).hlm = Vlm_eq;
            Uin.ell(l).emminus(n+1).hlm = Ulm_in;
            Vin.ell(l).emminus(n+1).hlm = Vlm_in; 

            Ueq.ell(l).emminus(n+1).t = DBeq.t;
            Veq.ell(l).emminus(n+1).t = DBeq.t;
            Uin.ell(l).emminus(n+1).t = DBin.t;
            Vin.ell(l).emminus(n+1).t = DBin.t; 
        else
            Ueq.ell(l).emm(n+1).hlm = Ulm_eq;
            Veq.ell(l).emm(n+1).hlm = Vlm_eq;
            Uin.ell(l).emm(n+1).hlm = Ulm_in;
            Vin.ell(l).emm(n+1).hlm = Vlm_in; 

            Ueq.ell(l).emm(n+1).t =  DBeq.t;
            Veq.ell(l).emm(n+1).t =  DBeq.t;
            Uin.ell(l).emm(n+1).t =  DBin.t;
            Vin.ell(l).emm(n+1).t =  DBin.t; 
        end
       
    end
    Ueq.ell(l).emminus(1).hlm = Ueq.ell(l).emm(1).hlm;
    Veq.ell(l).emminus(1).hlm = Veq.ell(l).emm(1).hlm;
    Uin.ell(l).emminus(1).hlm = Uin.ell(l).emm(1).hlm;
    Vin.ell(l).emminus(1).hlm = Vin.ell(l).emm(1).hlm; 

    Ueq.ell(l).emminus(1).t = Ueq.ell(l).emm(1).t;
    Veq.ell(l).emminus(1).t = Veq.ell(l).emm(1).t;
    Uin.ell(l).emminus(1).t = Uin.ell(l).emm(1).t;
    Vin.ell(l).emminus(1).t = Vin.ell(l).emm(1).t;  

    Utilde_eq = struct;
    Utilde_eq = DB_mode_rotate_timedep(Ueq,Uin.dyn,l,M,Utilde_eq,'direction','backward'); 
    Vtilde_eq = struct;
    Vtilde_eq = DB_mode_rotate_timedep(Veq,Vin.dyn,l,M,Vtilde_eq,'direction','backward'); 

    my_linewidth = 1.5;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;
    
    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','southwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',2)%,'Orientation','horizontal')
    hold on
    plot(Utilde_eq.ell(l).emm(M+1).t,real(Utilde_eq.ell(l).emm(M+1).hlm),'LineWidth',my_linewidth,'Color',[1,0,0],'DisplayName',sprintf('$\\Re\\left[\\tilde{U}^{\\rm eq}_{%d%d}\\right]$',l,M))
    plot(Uin.ell(l).emm(M+1).t,real(Uin.ell(l).emm(M+1).hlm),'LineStyle','--','LineWidth',my_linewidth,'Color',[0,1,0],'DisplayName',sprintf('$\\Re\\left[{U}^{\\rm incl}_{%d%d}\\right]$',l,M))
    plot(Utilde_eq.ell(l).emm(M+1).t,abs(Utilde_eq.ell(l).emm(M+1).hlm),'LineWidth',my_linewidth,'Color',[1,0,0,0.2],'DisplayName',sprintf('$|\\tilde{U}^{\\rm eq}_{%d%d}|$',l,M))
    plot(Uin.ell(l).emm(M+1).t,abs(Uin.ell(l).emm(M+1).hlm),'LineStyle','--','LineWidth',my_linewidth,'Color',[0,1,0,0.2],'DisplayName',sprintf('$|{U}^{\\rm incl}_{%d%d}|$',l,M))
    xlim([500,750])

    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','southwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',2)%,'Orientation','horizontal')
    hold on
    plot(Vtilde_eq.ell(l).emm(M+1).t,real(Vtilde_eq.ell(l).emm(M+1).hlm),'LineWidth',my_linewidth,'Color',[1,0,0],'DisplayName',sprintf('$\\Re\\left[\\tilde{V}^{\\rm eq}_{%d%d}\\right]$',l,M))
    plot(Vin.ell(l).emm(M+1).t,real(Vin.ell(l).emm(M+1).hlm),'LineStyle','--','LineWidth',my_linewidth,'Color',[0,1,0],'DisplayName',sprintf('$\\Re\\left[{V}^{\\rm incl}_{%d%d}\\right]$',l,M))
    plot(Vtilde_eq.ell(l).emm(M+1).t,abs(Vtilde_eq.ell(l).emm(M+1).hlm),'LineWidth',my_linewidth,'Color',[1,0,0,0.2],'DisplayName',sprintf('$|\\tilde{V}^{\\rm eq}_{%d%d}|$',l,M))
    plot(Vin.ell(l).emm(M+1).t,abs(Vin.ell(l).emm(M+1).hlm),'LineStyle','--','LineWidth',my_linewidth,'Color',[0,1,0,0.2],'DisplayName',sprintf('$|{V}^{\\rm incl}_{%d%d}|$',l,M))
    xlim([500,750])
return