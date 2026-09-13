function [AF,DYN,HPC] = hlm_plots(s,l,m,varargin)
    Dyn = 1;
    AF = 1;
    R = 0;
    I = 0;
    rot = 0;
    hpc = 0;
    emminus = 0;

    norm = sqrt((l+2).*(l+1).*l.*(l-1));

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'real'
                    R = 1;
                case 'imag'
                    I = 1;
                case 'rot'
                    rot = 1;
                case 'hpc'
                    hpc = 1;
                case 'emminus'
                    emminus = 1;
                case 'nodyn'
                    Dyn = 0;
                case 'noAF'
                    AF = 0;
                case 'parity'
                    parity = 1;
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    a = s.dyn.chi1(3);
    %if a<0
    %    s = invert_parity(s);
    %end
    A = abs(a)*10;

    I = rad2deg(pi/2 - sign(a).*s.dyn.th0);
    config_label = sprintf('a0%1.0fi%1.f',A,I);

    % =======================
    % COMPUTE ROTATED MODES IF NECESSARY
    % =======================
    if rot==1
        if hpc==0
            ellmax = l;
        else
            ellmax = 4;
        end
        sr = struct;
        for L=2:ellmax
            for M=-L:L
                sr = DB_mode_rotate_timedep(s,s.dyn,L,M,sr);
            end
        end
        clear s
        s = sr;
    end
    
    if m<0
        t = s.ell(l).emminus(-m+1).t;
        hlm = s.ell(l).emminus(-m+1).hlm;
        if emminus
            hlmminus = s.ell(l).emm(-m+1).hlm;
        end
    else
        t = s.ell(l).emm(m+1).t;
        hlm = s.ell(l).emm(m+1).hlm;
        if emminus
            hlmminus = s.ell(l).emminus(m+1).hlm;
        end
    end

    % =======================
    % QNM frequency
    % =======================
    addpath '/home/luca/repos/kerrorbitsolver/KerrDynamics' 
    a = s.dyn.chi1(3);
    if a>0
        sigma = kerr_FitKerrQNMs(l,m,a,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
        prefactor = 1;
    else
        sigma = -kerr_FitKerrQNMs(l,-m,abs(a),1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/');
        prefactor = 1;
    end
    omg_QNM = imag(sigma);

    
    % =======================
    % Dynamical quantities
    % =======================
    dyn = s.dyn;
    T = dyn.t;
    tLR = tLR_splined(dyn);
    tLSSO = tLSSO_splined(dyn);

    t_in = tLR - 200;
    t_end = tLR + 100;

    th_min = rad2deg(s.dyn.th0);
    th_max = rad2deg(pi-s.dyn.th0);

    % =======================
    % Amplitude-Frequency plots
    % =======================
    if AF==1
        f1 = figure;
        tiledlayout(2, 1,'TileSpacing', 'compact', 'Padding', 'compact');
        f1.Position(1) = f1.Position(1) + 300;
        nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,abs(hlm)./norm,'LineWidth',1,'Color','k','DisplayName',sprintf('$|\\Psi_{%d%d}|/\\nu$',l,m))
        if emminus==1
            plot(t,abs(hlmminus)./norm,'LineWidth',1.5,'LineStyle',':','Color','r','DisplayName',sprintf('$|\\Psi_{%d-%d}|/\\nu$',l,m))
        end
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(tLSSO,'Color',[.7 .7 .7],'LineStyle','-.','HandleVisibility','off')
        legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7)
        xlim([t_in,t_end]) 
        lims = ylim;
        ylim(lims*1.1)
        ylabel(' ') %Add phantom label to prevent cropping values on y axis when exporting

        nexttile
        set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        lnk = plot(t,freq(hlm,t),'LineWidth',1,'Color','k','DisplayName',sprintf('$\\omega_{%d%d}$',l,m));
        yline(omg_QNM,'LineStyle','--','HandleVisibility','off')
        if emminus==1
            plot(t,freq((-1).^l.*conj(hlmminus),t),'LineWidth',1.5,'LineStyle',':','Color','r','DisplayName',sprintf('$-\\omega_{%d-%d}$',l,m))
        end
        lnh = plot(T,prefactor*2*dyn.Omg_orb,'Color',[MyColors('g1'),0.5],'LineWidth',2,'DisplayName','$2\Omega_{\rm orb}$');
        %lnh1 = plot(T,2*(sin(dyn.th)).*DB_D1(dyn.th,dyn.t,4),'Color',[.8 .8 .8],'DisplayName','$2\Omega$');
        uistack(lnh,"bottom");
        uistack(lnk,'top');
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(tLSSO,'Color',[.7 .7 .7],'LineStyle','-.','HandleVisibility','off')
        xlabel('$u$','FontSize',labels_fontsize,'Interpreter','Latex');
        yline(0,'Color',[.9 .9 .9],'LineWidth',.25,'HandleVisibility','off')
        L = legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7);
        L.Direction = 'reverse';
        xlim([t_in,t_end])
        xlength = t_end-t_in;
        lims = ylim;
        ylabel(' ') %Add phantom label to prevent cropping values on y axis when exporting
        if abs(lims(2)-omg_QNM)<0.07
            lims = ylim;
            ylim([lims(1),lims(2)*1.1])
        end
        ylength = abs(lims(2)-lims(1));

        %{
        if strcmp(config_label,'a02i120')
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a05i120')
            ylim([-1,1])
            ylength = 2;
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.5,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(-omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + 0.31*xlength,-new_omg_QNM ,sprintf('$-\\omega_{%d-%d0} = -%.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a05i135')
            ylim([0,1])
            ylength = 1;
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a05i150')
            ylim([0,1])
            ylength = 1;
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif contains(config_label,'a09i1') 
            ylim([-1,1])
            ylength = 2;
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.9,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(-new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + 0.31*xlength,-new_omg_QNM + ylength/15,sprintf('$-\\omega_{%d-%d0} = -%.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        else
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        end
        %}

        if a>0
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a02i180')
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a02i150')
            ylim(1.1.*lims)
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a05i120')
            ylim([-1,1])
            ylength = 2;
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.5,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + 0.31*xlength,new_omg_QNM + ylength/15 ,sprintf('$\\omega_{%d-%d0} = %.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        elseif strcmp(config_label,'a05i180')
            text(t_in + xlength/3,omg_QNM + ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        elseif contains(config_label,'a09i120') 
            ylim([-1.5,1.5])
            ylength = 3;
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.9,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + xlength/3,new_omg_QNM + ylength/15,sprintf('$\\omega_{%d-%d0} = %.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        elseif contains(config_label,'a09i135') 
            ylim([-1,1.2])
            ylength = 2.2;
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.9,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + xlength/3,new_omg_QNM + ylength/15,sprintf('$\\omega_{%d-%d0} = %.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        elseif contains(config_label,'a09i150') 
            ylim([-1,1.5])
            ylength = 2.5;
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.9,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + xlength/3,new_omg_QNM + ylength/15,sprintf('$\\omega_{%d-%d0} = %.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        elseif contains(config_label,'a09i180') 
            ylim([-1.5,1.5])
            ylength = 3;
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$-\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
            new_omg_QNM = imag(kerr_FitKerrQNMs(2,2,0.9,1,'/home/luca/repos/kerrorbitsolver/KerrDynamics/'));
            yline(new_omg_QNM,'LineStyle','--','HandleVisibility','off')
            text(t_in + xlength/3,new_omg_QNM + ylength/15,sprintf('$\\omega_{%d-%d0} = %.5f$',l,m,new_omg_QNM),'Interpreter','latex','FontSize',12)
        else
            text(t_in + xlength/3,omg_QNM - ylength/15,sprintf('$\\omega_{%d%d0} = %.5f$',l,m,omg_QNM),'Interpreter','latex','FontSize',12)
        end
        AF = f1;
    else
        AF = figure;
    end

    axes_fontsize = 11;
    labels_fontsize = 16;
    
    % =======================
    % Trajectory and theta plots
    % =======================
    if Dyn==1
        idx_plunge = find(dyn.t>tLSSO_splined(dyn),1);
        f3 = figure;
        f3.Position(1) = 100;
        tl = tiledlayout(3, 5,'TileSpacing', 'compact', 'Padding', 'compact');

        ax1 = nexttile(tl,4,[2 2]);
        hold on
        plot(ax1,dyn.x(1:idx_plunge),dyn.y(1:idx_plunge),'Color','b','LineWidth',1)
        plot(ax1,dyn.x(idx_plunge:end),dyn.y(idx_plunge:end),'Color',MyColors('r1'),'LineWidth',1)
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        xlabel('$x$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel('$y$','Interpreter','latex','FontSize',labels_fontsize);
        xlim([-s.dyn.r0-.1,s.dyn.r0+.1])
        ylim([-s.dyn.r0-.1,s.dyn.r0+.1])

        ax2 = nexttile(tl,14,[1 2]);
        plot(ax2,T,rad2deg(dyn.th),'Color','r','LineWidth',1)
        yline(90,'Color',[.7 .7 .7],'LineStyle',':')
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        xlabel('$t$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel('$\theta\ (^\circ)$','Interpreter','latex','FontSize',labels_fontsize);
        ytickangle(45)
        if abs(s.dyn.th0-pi/2)/2>1e-5
            yticks([th_min,90,th_max])
            xlim([T(end)-150,T(end)+10])
        end
        temp = ylim;
        ylim([temp(1)-10,temp(2)+10])

        ax3 = nexttile(tl, 1,[3 3]);
        plot3(ax3,dyn.x(1:idx_plunge),dyn.y(1:idx_plunge),dyn.z(1:idx_plunge),'LineWidth',2,'Color','b')
        hold on
        plot3(ax3,dyn.x(idx_plunge:end),dyn.y(idx_plunge:end),dyn.z(idx_plunge:end),'LineWidth',2,'Color',MyColors('r1'))
        xlabel('$x$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
        ylabel('$y$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
        zlabel('$z$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
        if abs(max(dyn.z)-min(dyn.z))<1e-6
            zlim([-1,1])
        else 
            axis equal
        end
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5);
        grid(ax3,'on')

        if strcmp(config_label,'a09i0')
            view(ax3,-38,26)
        end

        xl = xlim;
        yl = ylim;
        zl = zlim;
        text(ax3,0.9*xl(1),0.9*yl(end),0.9*zl(end),sprintf('{\\tt %s}',config_label),'Interpreter','latex','FontSize',14)
        A4Width(f3)

        DYN = f3;
    else
        DYN = figure;
        close
    end

    % =======================
    % Real part plot
    % =======================
    if R==1
        f4 = figure;
        f4.Position(1:2) = [1300,500];
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5);
        plot(t,real(hlm),'LineWidth',1,'Color',MyColors('r1'))
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--')
        xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel(sprintf('$\\Re[h_{%d%d}]$',l,m),'Interpreter','latex','FontSize',labels_fontsize);
        xlim([t_in,t_end])
    end
    % =======================
    % Imaginary part plot
    % =======================
    if I==1
        f5 = figure;
        f5.Position(1:2) = [1300,100];
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5);
        plot(t,imag(hlm),'LineWidth',1,'Color',MyColors('b1'))
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--')
        xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel(sprintf('$\\Im[h_{%d%d}]$',l,m),'Interpreter','latex','FontSize',labels_fontsize);
        xlim([t_in,t_end])
    end

    % =======================
    % plus and cross polarizations
    % =======================
    if hpc==1
        ellmax = 4;
        [h,t] = DB_hpc(s,0*pi/4,0,'ellmax',ellmax);

        f6 = figure;
        f6.Position(1:2) = [1300,500];
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,real(h),'LineWidth',1,'Color',MyColors('r1'))
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','HandleVisibility','off')
        xline(tLSSO,'Color',[.7 .7 .7],'LineStyle','-.','HandleVisibility','off')
        xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel(sprintf('$D_L h_+/\\mu$'),'Interpreter','latex','FontSize',labels_fontsize);
        xlim([t_in-200,t_end])
        xlims = xlim;
        ylim([-1,1.2])
        ylims = ylim;
        xlength = xlims(2)-xlims(1);
        ylength = ylims(2)-ylims(1);
        text(gca,xlims(1) + 0.1*xlength, ylims(2) - ylength/15,'$(\Theta = 0,\Phi = 0)$','Interpreter','latex','FontSize',14) 
        text(gca,tLR - 45, ylims(1) + ylength/15,'${\rm LR} \rightarrow$','Interpreter','latex','FontSize',10)         
        text(gca,tLSSO + 3, ylims(1) + ylength/15,'$\leftarrow {\rm LSSO}$','Interpreter','latex','FontSize',10)         

        f7 = figure;
        f7.Position(1:2) = [1300,100];
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
        hold on
        plot(t,imag(h),'LineWidth',1,'Color',MyColors('b1'))
        xline(tLR,'Color',[.7 .7 .7],'LineStyle','--')
        xlabel('$u$','Interpreter','latex','FontSize',labels_fontsize);
        ylabel(sprintf('$D_L h_{\\times}/\\mu$'),'Interpreter','latex','FontSize',labels_fontsize);
        xlim([t_in-200,t_end])

        HPC = f6;
    else
        HPC = figure;
        close
    end


return

function s = invert_parity(s)
    s.dyn.x = -s.dyn.x;
    s.dyn.px = -s.dyn.px;
    [r,phi,theta,pr,pphi,ptheta] = DB_coords_cart2spherical(s.dyn.x,s.dyn.y,s.dyn.z,s.dyn.px,s.dyn.py,s.dyn.pz);    
    s.dyn.r = r;
    s.dyn.th = theta;
    s.dyn.phi = phi;
    s.dyn.pr = pr;
    s.dyn.pth = ptheta;
    s.dyn.pphi = pphi;

    for l=2:4
        for m=-l:l
            if m<0
                n = abs(m);
                s.ell(l).emminus(n+1).hlm = (-1).^m.*conj(s.ell(l).emminus(n+1).hlm);
            else
                s.ell(l).emm(m+1).hlm = (-1).^m.*conj(s.ell(l).emm(m+1).hlm);
            end
                s.ell(l).emminus(1).hlm = (-1).^m.*conj(s.ell(l).emminus(1).hlm);
        end
    end

return