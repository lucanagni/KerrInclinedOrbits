function DB_flux_spinless(a,n,d,varargin)
    %==========================================================================
    % Compare DB/KOS flux (fhat) with Hughes data (I should probably change script name)
    % a = central BH spin
    % n = Padé numerator
    % d = Padé denominator
    % 
    % if varargin is left empty computes DB flux
    % otherwise varargin must be the output of KoerrOrbitSolverClassy (for all the flags) and computes KOS flux
    % must satisfy n + d = 7
    %==========================================================================

    %case 'log_flux_7p5and8_Fit'
    my_linewidth = 2;
    axes_fontsize = 11;
    legend_fontsize = axes_fontsize+4;
    labels_fontsize = 20;

    ellmax_num = 8;
    mmin = 1;

    [hatF_num, x] = kerr_Flux_num(a,ellmax_num,mmin); 
    v=sqrt(x);

    if isempty(varargin)
        hatF  = hatF_DB(x,n,d);
        spin = SpinFlux(a,x);
        code_name = 'BCD Pade-resummed';
    else
        KOS = varargin{1};
        Omg         = x.^(3/2);
        r_Omg       = 1./x;
        r           = (r_Omg.^(3/2)-a).^(2/3);
        pph         = zeros_dH(1./r,a);
        Dyn.a       = a;
        Dyn.r       = r;
        Dyn.prstar  = 0;
        Dyn.pph     = pph;

        [~,out] = kerr_Hamiltonian(Dyn);

        v_phi  = Omg.*r_Omg;
        jNewt  = r_Omg.*v_phi;   % "Formal" Newtonian angular momentum
        jhat   = pph./jNewt;     %  Newton-normalized angular momentum

        [hatF, ~]     = kerr_Flux(KOS,x,Omg,1,out.H,jhat,8,a,r_Omg,v_phi);
        spin = 0;
        code_name = 'KOS factorized';
    end

    r_LSO = DB_separatrix(0, a);
    x_LSO = ((a+r_LSO.^(3/2)).^(-1)).^(2/3);
    v_LSO = sqrt(x_LSO);

    f=figure;
    size = f.Position(3:4);
    f.Position(3:4) = 1.1*size;
    clear size

    subplot(3,1,[1,2])
    ax1 = gca;
    set(ax1,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','box','on');
    hold on
    plot(v,hatF+spin,'DisplayName',code_name,'LineWidth',my_linewidth)
    plot(v,hatF_num,'DisplayName','Numerical (Hughes)','LineWidth',my_linewidth)
    xlim([0,v_LSO])
    legend('Location','northwest','Interpreter','latex','BackgroundAlpha',0.7,'FontSize',legend_fontsize)

    ylabel('$\hat{f}$','Interpreter','latex','FontSize',labels_fontsize)

    vlims = ylim;
    hlims = xlim;
    yl = vlims(2) - vlims(1);
    xl = hlims(2) - hlims(1);
    xpos = hlims(1) + 0.75*xl;
    ypos = vlims(1) + 0.85*yl;

    text(xpos,ypos,sprintf('a = %.2f',a),'FontSize',legend_fontsize,'Parent',ax1,'Interpreter','latex')
    hold off

    subplot(3,1,3)   
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','box','on');
    hold on
    plot(v,(hatF_num-hatF-spin)./(hatF_num),'LineWidth',my_linewidth)

    xlim([0,v_LSO])
    grid on 
    xlabel('$v_\omega$','Interpreter','latex','FontSize',labels_fontsize)
    ylabel('$\Delta\hat{f}/\hat{f}$','Interpreter','latex','FontSize',labels_fontsize)

return

function [hatf,x] = kerr_Flux_num(a,ellmax,mmin)
    fluxdir = DB_GetBaseDir('fluxesFD'); %copied from eob_ecceeob_eccentric_bis/TeukodeTests/MatlabScripts/GetBasedir.m
    if a==0
        astr = 'a0.0';
    else
        astr = ['a', strrep(num2str(a),'-','m')];
    end
    if contains(fluxdir,'\')
        trail = '\';
    else
        trail = '/';
    end
    mydir   = [fluxdir,  astr, trail];
    if ~isfolder(mydir)
        error('No data found for a=%f\n', a)
    end
    hatf = 0;
    for l=2:ellmax
        lstr = num2str(l);
        for m=mmin:l
            mstr    = num2str(m);
            data_I  = load([mydir, astr, '_eps1e-14_I.Yflux_', lstr, mstr]);
            data_O  = load([mydir, astr, '_eps1e-14_O.Yflux_', lstr, mstr]);
            hatf_lm = [data_O(1:end-800,4); data_I(:,4)];
            hatf    = hatf + hatf_lm;
        end
    end
    v = [data_O(1:end-800,1); data_I(:,1)];
    x = v.^2;
return

function y = SpinFlux(a,x)

    q=1000;
    nu = q/(1+q).^2;
    X1 = q/(1+q);
    X2 = 1/(1+q);

    chi1 = [0 0 a];
    chi2 = [0 0 0];
    l = [0 0 1];

    f3so  = -0.25.*(11.*X1 + 5.*X2).*X1.*dot(l,chi1)...
        -0.25.*(11.*X2 + 5.*X1).*X2.*dot(l,chi2);
    f4ss  = nu./48.*(289.*dot(l,chi1).*dot(l,chi2) - 103.*dot(chi1,chi2));

    y = f3so.*x.^(3./2)+f4ss.*x.^2;
return 

function y = hatF_DB(x,n,d)
    q=1000;
    nu         = q/(1+q).^2;
    eulergamma = 0.57721566490153286061;
    f2    = -1247./336 -35./12.*nu;
    f3    = 4.*pi;
    f4    = -44711./9072 + 9271./504.*nu + 65./18.*nu.^2;
    f5    = -(8191./672 + 583./24.*nu).*pi;
    f6    = 6643739519./69854400 + 16./3.*pi.^2 - 1702./105.*eulergamma...
        -(134543./7776 - 41./48.*pi.^2).*nu - 94403./3024.*nu.^2 - 775./324.*nu.^3;
    fl6   = -1712./105;
    f7    = -(16285./504 - 214745./1728.*nu - 193385./3024.*nu.^2).*pi;

    %y = (1 + f2.*x + (f3).*x.^(3/2) + (f4).*x.^2 ...
        %+ f5.*x.^(5/2) + (f6 + fl6.*log(4.*x.^(1/2))).*x.^3 + f7.*x.^(7/2));
    y = DB_pade(n,d ,f2,f3,f4,f5,f6,fl6,f7,0,0,x.^(1/2));
return 

function pph = zeros_dH(u,a)
    if a==0
        u2 = u.*u;
        pph  = ((u+(-3).*u2).^(-1)).^(1/2);
    else
        u2  = u.*u;
        u3  = u.*u2;
        u4  = u.*u3;
        u5  = u.*u4;
        u6  = u.*u5;
        u7  = u.*u6;
        u8  = u.*u7;
        u9  = u.*u8;
        u10 = u.*u9;
        u11 = u.*u10;
        u12 = u.*u11;
        u13 = u.*u12;
        a2  =  a*a;
        a4  = a2*a2;
        if a>0
            pm = -1;
        else
            pm =  1;
        end
        pph = ((a2.*(2+u.^(-1))+u.^(-3)).^4.*((-4).*a2+(1+(-3).*u).^2.*u.^(-3)) ...
            .^(-1).*u3.*(2.*a2.*(1+2.*u).^(-1).*((-2)+(u.^(-2)+a2.*(1+2.*u)) ...
            .^(1/2)).*(2+(u.^(-2)+a2.*(1+2.*u)).^(1/2)).*u10.*(3+a2.*u2).^2.*( ...
            1+a2.*(1+2.*u).*u2).^(-5)+2.*pm.*(a2.*u13.*(1+a2.*(1+2.*u).*u2).^( ...
            -8).*(3+(-6).*u+4.*a2.*u2+(-2).*a2.*u3+a4.*u4).^2).^(1/2)+(-1).*( ...
            1+a2.*(1+2.*u).*u2).^(-5).*(1+2.*a2.*(1+(-2).*u).*u2+a4.*u4).^2.* ...
            u6+(-1).*(1+2.*u).^(-1).*((-2)+(u.^(-2)+a2.*(1+2.*u)).^(1/2)).*(2+ ...
            (u.^(-2)+a2.*(1+2.*u)).^(1/2)).*(1+a2.*(1+2.*u).*u2).^(-5).*((-1)+ ...
            a2.*u3).*(1+2.*a2.*(1+(-2).*u).*u2+a4.*u4).*u7)).^(1/2);
    end
return