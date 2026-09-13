function hlm_EOB = DB_anal_wave(DB,l,m,varargin)

% ===================================================================================================
% Compute analytical waveforms using Nretonian multipoles and PN corrections
% ===================================================================================================

flags.logscale_switch = 0; % logscale in relative errors. I think it can be removed
flags.EOB_switch = 1; % old. Can be removed
flags.Newt_switch = 0; % turn on to have Newtonian waveform (should be used with flags.r = 1 too)
flags.analytical_label = 'EOB'; % old.
flags.r = 0; % turn on to use r instead of r_Omg
flags.source = 1; % turn off to compute hlm without effective source
flags.spin = 0; % was used in the rholm, probably useless now
flags.eccentric = 0; % turn on in eccentric orbit
flags.geod = 0; % turn on if geodesic motion
flags.newResumPN = 0; % PN order of the "non-equatorial" rholm. Not really working right now
flags.recap = 0; %only needed for the scipt RecapPlots which exports plots for Figs. 3-5 of PaperII. Sources the correct Schwarzschild DBamics
flags.rhoresummation = 'iotaexp'; %resummation strategy for the rholm. Choices are 'iotaexp','fulliota','equatorial','flm'
flags.deltaresummation = 'iotaexp'; %resummation strategy for the deltalm. Choices are 'iotaexp','fulliota','equatorial'
flags.pade = 0; %apply pade resummation to flm
flags.polarmultipole=0; %turn on to use polar expression on h21 odd multipole
flags.plot = 1; % turn off to avoid plotting wave
flags.newf21incl = 0; %Uses new form of f21odd_incl (attempts to factor out singularity)
flags.f21odd1p5PN = 0; % Include up to 1.5PN order in inclined f21odd
flags.f21odd3PN = 0; % Include up to 3PN order in inclined f21odd
flags.f21odd4p5PN = 0; % Include up to 4.5 PN oder in inclined f21odd

% read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'log'
                    flags.logscale_switch = 1;
                case 'newt'
                    flags.Newt_switch = 1;
                    flags.analytical_label = 'N+S';
                case 'nopn'
                    flags.EOB_switch = 0;
                case 'output'
                    i = i+1;
                    mystruct = varargin{i};
                case 'rhotate'
                    flags.rotatePN = 1;
                case 'r'
                    flags.r = 1;
                case 'nosource'
                    flags.source = 0;
                    flags.analytical_label = erase(flags.analytical_label,"+S");
                case 'spin'
                    flags.spin = 1;
                case 'eccentric'
                    i = i + 1;
                    flags.eccentric = 1;
                    s = varargin{i};
                case 'geod'
                    flags.geod = 1;
                case 'newresumpn'
                    i = i + 1;
                    flags.newResumPN = varargin{i};
                case 'recap'
                    flags.recap = 1;
                case 'rhoresum'
                    i = i + 1;
                    flags.rhoresummation = varargin{i};
                case 'deltaresum'
                    i = i + 1;
                    flags.deltaresummation = varargin{i};
                case 'pade'
                    flags.pade = 1;
                case 'polar'
                    flags.polarmultipole = 1;
                case 'noplot'
                    flags.plot = 0;
                case 'sing'
                    flags.newf21incl = 1;
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

my_linewidth = 2.0;
axes_fontsize = 12;
legend_fontsize = axes_fontsize+2;
labels_fontsize = 20;


%function [full_wave,Ulm,Vlm] = generic_wave_cart(DB,l,m,flags)

r0    = 1.213061319425267e+00;   % 2/sqrt(e);

Heff = DB.Heff;
jhat = DB_jhat(DB,flags);
%jhat = DB.pphi./(rOmg.^2.*Omega);

Omega = DB.Omg;
if m~=0
    Tail = DB_Tail(l,m.*Omega,m.*Omega,r0);
else
    Tail = 1;
end

rho = DB_rholm(DB,flags);
rho_incl = DB_rholm_iotaexp(DB,flags);
rho_full = DB_rholm_fulliota(DB,flags);
rho_hybr = DB_rholm_hybrid(DB,flags);
delta = DB_deltalm(DB);
delta_incl = DB_deltalm_iotaexp(DB,flags);
delta_full = DB_deltalm_fulliota(DB,flags);
delta_hybr = DB_deltalm_hybrid(DB,flags);

[Ulm,Vlm] = DB_Multipoles(DB,l,m,flags);
%[Ulm,Vlm] = DB_MultipolesPolar(DB,l,m,flags);

if strcmp(flags.rhoresummation,'equatorial')
    rho21e = rho.('l2m1');
    rho21o = rho21e;
    rho22e = rho.('l2m2');
    rho22o = rho22e;
    myell = l;
elseif strcmp(flags.rhoresummation,'iotaexp')
    rho21e = rho_incl.('l2m1e');
    rho21o = rho_incl.('l2m1o');
    rho22e = rho_incl.('l2m2e');
    rho22o = rho_incl.('l2m2o');
    myell = l;
elseif strcmp(flags.rhoresummation,'fulliota')
    rho21e = rho_full.('l2m1e');
    rho21o = rho_full.('l2m1o');
    rho22e = rho_full.('l2m2e');
    rho22o = rho_full.('l2m2o');
    myell = l;
elseif strcmp(flags.rhoresummation,'flm')
    rho21e = flm.('l2m1e');
    rho21o = flm.('l2m1o');
    rho22e = flm.('l2m2e');
    rho22o = flm.('l2m2o');
    myell = 1;
elseif strcmp(flags.rhoresummation,'hybrid')
    rho21e = rho_hybr.('l2m1e');
    rho21o = rho_hybr.('l2m1o');
    rho22e = rho_hybr.('l2m2e');
    rho22o = rho_hybr.('l2m2o');
    myell = 1;
elseif strcmp(flags.rhoresummation,'none')
    rho21e = 1;
    rho21o = 1;
    rho22e = 1;
    rho22o = 1;
    myell = 1;
else
    error('Invalid resummation flag')
end

if strcmp(flags.deltaresummation,'equatorial')
    delta21e = delta.('l2m1');
    delta21o = delta21e;
    delta22e = delta.('l2m2');
    delta22o = delta22e;
elseif strcmp(flags.deltaresummation,'iotaexp')
    delta21e = delta_incl.('l2m1e');
    delta21o = delta_incl.('l2m1o');
    delta22e = delta_incl.('l2m2e');
    delta22o = delta_incl.('l2m2o');
elseif strcmp(flags.deltaresummation,'fulliota')
    delta21e = delta_full.('l2m1e');
    delta21o = delta_full.('l2m1o');
    delta22e = delta_full.('l2m2e');
    delta22o = delta_full.('l2m2o');
elseif strcmp(flags.deltaresummation,'hybrid')
    delta21e = delta_hybr.('l2m1e');
    delta21o = delta_hybr.('l2m1o');
    delta22e = delta_hybr.('l2m2e');
    delta22o = delta_hybr.('l2m2o');
elseif strcmp(flags.deltaresummation,'none')
    delta21e = 0;
    delta21o = 0;
    delta22e = 0;
    delta22o = 0;
else
    error('Invalid resummation flag')
end

if flags.source==1
    Ulm = Ulm.*Heff;
    Vlm = Vlm.*jhat;
end

if flags.Newt_switch
    out = 1/sqrt(2).*(Ulm - 1i*Vlm);
else
    %out = 1/sqrt(2).*(Ulm - 1i*Vlm).*(rholm).^(l).*Tail.*exp(1i.*deltalm);
    if m==0
        out = 1/sqrt(2).*(Ulm.*(rho_incl.('l2m0e')).^(l) - 1i*Vlm.*(rho_incl.('l2m0o')).^(l)).*Tail; %for m=0 modes%
    elseif mod(l+m,2)==1
        if flags.polarmultipole
            %[~,Vlm_p] = DB_MultipolesPolar(DB,l,m,flags);
            %Vlm = Vlm_p.*jhat;
        end
        out = 1/sqrt(2).*(Ulm.*(rho21e).^(myell).*exp(1i.*delta21e) - 1i*Vlm.*(rho21o).^(myell).*exp(1i.*delta21o)).*(Tail); %for l+m=odd modes%
        %figure
        %plot(DB.t,angle(Vlm))
        %pause
    else
       out = 1/sqrt(2).*(Ulm.*(rho22e).^(myell).*exp(1i.*delta22e) - 1i.*Vlm.*(rho22o).^(myell).*exp(1i.*delta22o)).*(Tail); %for l+m=even modes%
        %out = 1/sqrt(2).*(Ulm.*(rho_incl.('l2m2e')).^(l).*exp(1i.*deltalm)- 1i.*Vlm.*(rho_incl.('l2m2o')).^(l).*exp(1i.*delta_incl.('l2m2o'))).*Tail;
    end

end

hlm_EOB = out;

if flags.geod
    tLR = DB.t(end);
else
    tLR = tLR_splined(DB);
end
T = DB.t;

% ==================================== PRODUCE PLOT =====================================
iota = (pi/2-min(DB.th));
phi = atan2(DB.y.*cos(iota),DB.x);

if flags.plot
    f = figure;
    f.Position(3:4) = f.Position(3:4)*1.2;
    %f.Position(3) = f.Position(3).*2;
    t_in = tLR-350;
    t_end = tLR+15;
    label = KerrLabel(DB);
    wf_label = sprintf('h_{%d%d}',l,m);

    ax1 = gca;
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',1)%,'Orientation','horizontal')
    xlim([t_in,t_end])
    ylabel(sprintf('$\\Re[%s]/\\nu$',wf_label),'Interpreter','latex','Fontsize',labels_fontsize)
    hold on

    plot(T,real(hlm_EOB),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Re\\left[%s\\right]$',wf_label),'Color','r')
    %plot(T,abs(w_nu),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s|$',wf_label),'Color',[1 0 0 .2])
    plot(T,imag(hlm_EOB),'LineWidth',my_linewidth,'DisplayName',sprintf('$\\Im\\left[%s\\right]$',wf_label),'LineStyle','-','Color','b')
    %plot(tau,abs(w_anal),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s^{\\rm %s}|$',wf_label,flags.analytical_label),'LineStyle','--','Color',[1 0 0 .2])
    plot(T,mod(phi,1.*pi)./10);
    %xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',1.5*my_linewidth,'HandleVisibility','off')

    xs1 = xlim;
    ys1 = ylim;
    if abs(ys1(2)+ys1(1))>1e-3
        ys1(2) = min(abs(ys1));
        ys1(1) = -ys1(2);
    end
    ys1=[-m,m];

    set(ax1,'XTickLabel',[])
    text(xs1(2)-abs(xs1(2)-xs1(1))./3,ys1(2) - abs(ys1(2)-ys1(1))./10,label,'Interpreter','latex','FontSize',16)

    f.Position(3) = f.Position(3).*2;
    %xlim(ax1,[2200,3011])
    %xlim(ax2,[2200,3011])
end
return
