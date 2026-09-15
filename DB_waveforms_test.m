function mystruct = DB_waveforms_test(a,iota,l,m,varargin)

% ===================================================================================================
% Compute analytical waveforms using Nretonian multipoles and PN corrections
% varargin arguments:
% log -> use logscale in relative difference plots%
% newt: plot newtonian waveform (no corrections)
% maggiore: use the Definition of vurrent quadrupole founf in Maggiore's book (old flag, should remove it now. Also, I think Maggiore's formula is lacking a factor 2)
% output: return structure containing the (L,m) mode of waveform. Useful in DB_TEMP_build_analytical_structure
% ===================================================================================================


mystruct = struct;
s = struct;

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
flags.recap = 0; % only needed for the scipt RecapPlots which exports plots for Figs. 3-5 of PaperII. Sources the correct Schwarzschild dynamics
flags.rhoresummation = 'iotaexp'; % resummation strategy for the rholm. Choices are 'iotaexp','fulliota','equatorial','flm'
flags.deltaresummation = 'iotaexp'; % resummation strategy for the deltalm. Choices are 'iotaexp','fulliota','equatorial'
flags.pade = 0; % apply pade resummation to flm
flags.polarmultipole=0; % turn on to use polar expression on h21 odd multipole
flags.newf21incl = 0; % Uses new form of f21odd_incl (attempts to factor out singularity)
flags.f21odd1p5PN = 0; % Include up to 1.5PN order in inclined f21odd
flags.f21odd3PN = 0; % Include up to 3PN order in inclined f21odd
flags.f21odd4p5PN = 0; % Include up to 4.5 PN oder in inclined f21odd
flags.NewAlphaTail = 0; % Use different tails for even and odd components (adds an alpha parameter in the real exponential)


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
                case 'sing'
                    flags.newf21incl = 1;
                case '1p5'
                    flags.f21odd1p5PN = 1;
                    flags.f21odd3PN = 0;
                    flags.f21odd4p5PN = 0;
                case '3'
                    flags.f21odd1p5PN = 1;
                    flags.f21odd3PN = 1;
                    flags.f21odd4p5PN = 0;
                case '4p5'
                    flags.f21odd1p5PN = 1;
                    flags.f21odd3PN = 1;
                    flags.f21odd4p5PN = 1;
                case 'alpha'
                    flags.NewAlphaTail = 1;
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

% ================ COMPARE ANALYTICAL AND NUMERICAL ================
waveformdir = '/home/luca/waveforms';

if a~=0
    mydir = finddir(a,iota);
    dir = sprintf('%s/K/plunge/%s/wf.mat',waveformdir,mydir);
else
    th = abs(90-iota);
    th0 = num2str(th);
    if flags.geod
        base = 'geod';
    else
         base = 'plunge';
    end
    dir = sprintf('%s/S/%s/th0_%s/wf.mat',waveformdir,base,th0);
end

if ~flags.eccentric
    clear s;
    load(dir)
end

DB = s.dyn;
mystruct.dyn = s.dyn;

[w_anal,Ulm, Vlm] = generic_wave_cart(DB,l,m,flags);

if flags.geod
    tLR = DB.t(end);
else
    tLR = tLR_splined(DB);
end
%tau = DB.t-tLR;
tau = DB.t;

if m>=0
    w_num = s.ell(l).emm(m+1).hlm;
    %T = s.ell(l).emm(m+1).t-tLR;
    T = s.ell(l).emm(m+1).t;

    mystruct.ell(l).emm(m+1).hlm = w_anal;
    mystruct.ell(l).emm(m+1).t = DB.t;
else
    w_num = s.ell(l).emminus(-m+1).hlm;
    %T = s.ell(l).emminus(-m+1).t-tLR;
    T = s.ell(l).emminus(-m+1).t;

    mystruct.ell(l).emminus(-m+1).hlm = w_anal;
    mystruct.ell(l).emminus(-m+1).t = DB.t;
end

splined_num = spline(T,w_num,tau);
DeltaA = (abs(w_anal)-abs(splined_num))./abs(splined_num);
DeltaPhi = (-unwrap(angle(w_anal))) - (-unwrap(angle(splined_num)));

M = mean(DeltaPhi(200:end-200));
while abs(M)>pi/2
    if M>0
        DeltaPhi = DeltaPhi - pi;
    else
        DeltaPhi = DeltaPhi + pi;
    end
    M = mean(DeltaPhi(200:end-200));
end

% set appropriate label for the analytical, showing the additional rho contribution and the PN order
%{
% This was useful when the different-parity rholm's were computed differently, mixing equatorial and inclined results. Now it sis probably obsolete
if iota==0
    anlabel = 'Analytical';
else
    if mod(l+m,2)==0
        anlabel = sprintf('Analytical (with $\\rho_{%d%d}^{(1)}$@%dPN',l,m,flags.newResumPN); %even modes need odd correction
    else
        anlabel = sprintf('Analytical (with $\\rho_{%d%d}^{(0)}$@%dPN',l,m,flags.newResumPN); %odd modes need even correction
    end

end
%}
%anlabel = sprintf('Analytical $(\\rho_{%d%d} \\texttt{ %s},\\ \\delta_{%d%d} \\texttt{ %s})$',l,m,flags.rhoresummation,l,m,flags.deltaresummation);
anlabel = sprintf('Analytical');
if flags.eccentric
    legend_pos = 'best';
    delta_x0 = 750;
else
    legend_pos = 'northwest';
    delta_x0 = 350;
end


f = figure;
f.Position(3:4) = f.Position(3:4)*1.2;
f.Position(3) = f.Position(3).*1.2;
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','tight');
t_in = tLR-delta_x0;
t_end = tLR+15;
label = KerrLabel(DB);
%label = sprintf('$a = 0.5, \\iota = 45^\\circ$');
wf_label = sprintf('h_{%d%d}',l,m);

ax1 = nexttile(tl,[2,1]);
set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',1)%,'Orientation','horizontal')
xlim([t_in,t_end])
ylabel(sprintf('$\\Re[%s]/\\nu$',wf_label),'Interpreter','latex','Fontsize',labels_fontsize)
hold on

plot(T,real(w_num),'LineWidth',my_linewidth,'DisplayName','Numerical','Color','k')
%plot(T,abs(w_num),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s|$',wf_label),'Color',[1 0 0 .2],'HandleVisibility','off')
plot(tau,real(w_anal),'LineWidth',my_linewidth,'DisplayName',anlabel,'LineStyle','--','Color',MyColors('r1'))
%plot(tau,abs(w_anal),'LineWidth',my_linewidth,'DisplayName',sprintf('$|%s^{\\rm %s}|$',wf_label,flags.analytical_label),'LineStyle','--','Color',[1 0 0 .2],'HandleVisibility','off')
xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',1.5*my_linewidth,'HandleVisibility','off')

xs1 = xlim;
ys1 = ylim;
if abs(ys1(2)+ys1(1))>1e-3
    ys1(2) = min(abs(ys1));
    ys1(1) = -ys1(2);
end
ys1=[-m,m];

set(ax1,'XTickLabel',[])
text(xs1(2)-abs(xs1(2)-xs1(1))./3,ys1(2) - abs(ys1(2)-ys1(1))./10,label,'Interpreter','latex','FontSize',16)

ax2 = nexttile;
set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
legend('Location','northwest','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.9,'NumColumns',2)
xlim([t_in,t_end])
hold on

if flags.logscale_switch==1
    plot(tau,abs(DeltaPhi),'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
    plot(tau,abs(DeltaA),'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',my_linewidth,'HandleVisibility','off')

    yscale log
    yticks([1e-4,1e-3,1e-2,1e-1,1])
    grid off
    grid on
else
    plot(tau,DeltaPhi,'LineWidth',my_linewidth,'LineStyle','-','Color',[0.25, 0.88, 0.82],'DisplayName','$\Delta \phi$')
    plot(tau,DeltaA,'LineWidth',my_linewidth,'LineStyle','--','Color',[0.98, 0.70, 0.38],'DisplayName','$\Delta A/A$')
    xline(tLR,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',my_linewidth,'HandleVisibility','off')

    xlim([t_in+50,t_end-100])
    ylim auto
    ys = ylim;
    xlim([t_in,t_end])
    ylim(ys*1.5)
    grid on
end
xlabel('$u$','FontSize',labels_fontsize,'Interpreter','latex')

if ~flags.geod
    if abs(tLR-tLR_splined(s.dyn))>1e-5
    fprintf('WARNING: analytical tLR is %.5f while numerical tLR is %.5f',tLR,tLR_splined(s.dyn))
    end
end

ys2 = ylim;
xs2 = xlim;
%ys2
if ~flags.Newt_switch && abs(ys2(2)+ys2(1))/ys2(2) > 1e-3
    ys2(2) = 1.5*min(abs(ys2));
    ys2(1) = -ys2(2);
end
%ys2
xlim(ax1,[t_in,tau(end)-5])
xlim(ax2,[t_in,tau(end)-5])
ylim(ax1,ys1)
%ylim(ax2,ys2)

if flags.eccentric==1 && flags.geod==0
    %xlim(ax1,[xs1(1),xs1(2)+100])
    %ylim(ax1,ys1)

    %xlim(ax2,[xs2(1),xs2(2)+100])
    %ylim(ax2,ys2)
    xlim(ax1,[7200,tLR+100])
    xlim(ax2,[7200,tLR+100])
    xs1 = xlim(ax1);
    ys1 = xlim(ax1);
    text(xs1(2)-abs(xs1(2)-xs1(1))./3,ys1(2) - abs(ys1(2)-ys1(1))./10,label,'Interpreter','latex','FontSize',16)

end

return

function [full_wave,Ulm,Vlm] = generic_wave_cart(DB,l,m,flags)

fprintf('CARTESIAN FORMULA\n')

r0    = 1.213061319425267e+00;   % 2/sqrt(e);

Heff = DB.Heff;
jhat = DB_jhat(DB,flags);
%jhat = DB.pphi./(rOmg.^2.*Omega);

Omega = DB.Omg;
if m~=0
    if flags.NewAlphaTail
        if m==1
            Tail0 = DB_Tail(l,m.*Omega,m.*Omega,r0,2);
            Tail1 = DB_Tail(l,m.*Omega,m.*Omega,r0,1);
        elseif m==2
            Tail0 = DB_Tail(l,m.*Omega,m.*Omega,r0,1);
            Tail1 = DB_Tail(l,m.*Omega,m.*Omega,r0,1/2);
        else
            disp('m value invalid. No Tail computed. Setting Tail = 1')
            Tail1 = 1;
            Tail0 = 1;
        end
    else
        Tail0 = DB_Tail(l,m.*Omega,m.*Omega,r0);
        Tail1 = Tail0;
    end
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

flm = DB_flm(DB,flags);

element = sprintf('l%dm%d',l,abs(m));
if mod(m,2)==0
    parity = 'o';
else
    parity = 'e';
end

element_new = sprintf('l%dm%d%s',l,m,parity);

if m>0
    rholm = rho.(element);
    deltalm = delta.(element);
elseif m<0
    rholm = rho.(element);
    deltalm = -delta.(element);
else
    rholm = 1;
    deltalm = 0;
end

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
    %[Ulm,Vlm] = V21(DB);
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
            [~,Vlm_p] = DB_MultipolesPolar(DB,l,m,flags);
            Vlm = Vlm_p.*jhat;
        end
        out = 1/sqrt(2).*(Ulm.*(rho21e).^(myell).*exp(1i.*delta21e).*Tail0 - 1i*Vlm.*(rho21o).^(myell).*exp(1i.*delta21o).*Tail1); %for l+m=odd modes%
    else
       out = 1/sqrt(2).*(Ulm.*(rho22e).^(myell).*exp(1i.*delta22e).*Tail0 - 1i.*Vlm.*(rho22o).^(myell).*exp(1i.*delta22o).*Tail1); %for l+m=even modes%
    end

%{
    F = figure;
    F.Position(3) = F.Position(3).*2;
    hold on
    T = DB.t;
    plot(T,abs(rho_hybr.debug.l2m1o_orb_incl),'DisplayName','|f21_{incl}|')
    plot(T,abs(Vlm),'DisplayName','|V21|')
    plot(T,abs(Vlm.*rho_hybr.debug.l2m1o_orb_incl),'DisplayName','|f21_{incl}*V21|')
    %plot(DB.t,abs(rho_hybr.debug.l2m2e_orb_incl),'DisplayName','|f22_{e}|','LineWidth',2)
    %plot(DB.t,abs(rho_hybr.debug.l2m2o_orb_incl),'DisplayName','|f22_{o}|','LineWidth',2)
    %plot(DB.t,abs(rho_hybr.debug.l2m1e_orb_incl),'DisplayName','|f21_{e}|','LineWidth',2)
    %plot(DB.CurrTestt,abs(rho_hybr.debug.l2m1o_orb_incl),'DisplayName','|f21_{o}|','LineWidth',2)
    legend()
    yscale log
    xlim([800,1100])
    grid on
    close
%}
end

full_wave = out;


return

function dir = finddir(a,iota)
    if iota>90
        sign = '-';
    else
        sign = '';
    end
    th = abs(90-iota);
    th0 = num2str(th);

    dir = sprintf('th0_%s/a%s0%d',th0,sign,a);
return

function out = wigner_D_function(l,m1,m,alpha,beta,gamma)  %as written, this computes D^l_{m1,m} (order of indices is important)
    cth  = cos(beta*0.5);
    sth  = sin(beta*0.5);
    norm = sqrt( (factorial(l+m1) * factorial(l-m1) * factorial(l+m) * factorial(l-m)) );
    ki   = max(0,m-m1);
    kf   = min(l+m,l-m1);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+m-k) * factorial(l-m1-k) * factorial(k-m+m1) );
        dWig = dWig + div*( (-1).^(k) * cth.^(2*l+m-m1-2*k) * sth.^(2*k+m1-m) );
    end
    dWig_normalized = dWig*norm;
    out = exp(1i.*alpha.*m1).*dWig_normalized.*exp(1i.*gamma.*m);
return

function [U21,V21] = V21(DB)
    % This was needed when attempting to remove the singularity via a different factorization of the Newtonian part, but since we now know that the singularity comes also from the series
    % expansion when computing flm this should not be needed anymore
    r = DB.r;
    t = DB.t;

    x = (DB.r.*DB.Omg).^2;
    %phi = DB.phi;
    iota = (pi/2-min(DB.th));
    phi = atan2(DB.y.*cos(iota),DB.x);
    theta = DB.th;
    %{
    phidot = DB_D1(phi,t,4);
    phi2dot = DB_D1(phidot,t,4);
    phi3dot = DB_D1(phi2dot,t,4);
    thetadot = DB_D1(theta,t,4);
    theta2dot = DB_D1(thetadot,t,4);
    theta3dot = DB_D1(theta2dot,t,4);
%}
    U21 = -8.*exp(-2i.*phi).*sqrt(pi/5).*x;
    V21 = -(4/3).*1i.*sqrt(pi/5).*x.^(3/2).*exp(-1i.*phi);
return
