function fig = DB_conv_test(ID,ny,emm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dir: directory containing Teuk output in the form '/home/luca/waveforms/teukode_convergence_study/K/teuk_HH10_geod_a0.2000_r06.500_th030.000_q1e+03_3601x81_m'
% n_whatever = number of angular points of the 3 resolutions (has to coincide with the number appearing in the directory)
% scaling, rate -> n^p
% Title should be either 'Amplitude' or 'Phase'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    autosave = 0;
    linewidth = 2;
    rate = 2;
    ell = 4;
    %emm = 4;

    basedir = '/home/luca/waveforms/convergence/'

    dt = 700;
    dir = sprintf('%s%s',basedir,ID);

    scaling = 1.5;

    n_low = ny(1);
    n_mid = ny(2);
    n_high = ny(3);

    dir_low = dir;
    dir_mid = strrep(dir_low,num2str(n_low),num2str(n_mid));
    dir_high = strrep(dir_low,num2str(n_low),num2str(n_high));

    wlow = DB_get_wf(dir_low,ell,emm);
    wmid = DB_get_wf(dir_mid,ell,emm);
    whigh = DB_get_wf(dir_high,ell,emm);

    h_low = wlow.h;
    t_low = wlow.t;
    h_mid = wmid.h;
    t_mid = wmid.t;
    h_high = whigh.h;
    t_high = whigh.t;

    Ampl_high = abs(h_high);
    Ampl_mid = abs(h_mid);
    Ampl_low = abs(h_low);

    Phase_high = -unwrap(angle(h_high));
    Phase_mid = -unwrap(angle(h_mid));
    Phase_low = -unwrap(angle(h_low));

    Ampl_hm = Ampl_mid - spline(t_high,Ampl_high,t_mid);
    Ampl_lm = Ampl_low - spline(t_mid,Ampl_mid,t_low);

    Phase_hm = Phase_mid - spline(t_high,Phase_high,t_mid);
    Phase_lm = Phase_low - spline(t_mid,Phase_mid,t_low);

    Phase_hm = Phase_hm - 2.*pi.*round((mean(Phase_hm)./(2*pi)));
    Phase_lm = Phase_lm - 2.*pi.*round((mean(Phase_lm)./(2*pi)));

    x0 = 850;
    x1 = 1250;

    th = '{\theta}';

    fig = figure;
    tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    nexttile
    ax1 = gca;
    set(ax1,'XMinorTick','on','YMinorTick','on','FontSize',10,'FontName','Times','box','on');
    hold on
    plot(t_low+dt,abs(Ampl_lm),'DisplayName',sprintf('$|N_%s=%d - N_%s=%d|$',th,n_low,th,n_mid),'LineWidth',linewidth,'Color',MyColors('b1'))
    plot(t_mid+dt,abs(Ampl_hm),'DisplayName',sprintf('$|N_%s=%d - N_%s=%d|$',th,n_mid,th,n_high),'LineWidth',linewidth,'Color',MyColors('r1'))
    plot(t_mid+dt,scaling^rate.*abs(Ampl_hm),'DisplayName',sprintf('$%g^{%d}|N_%s=%d - N_%s=%d|$',scaling,rate,th,n_mid,th,n_high),'LineStyle','--','LineWidth',linewidth,'Color',MyColors('g1'))
    xlim([x0,x1])
    ylabel('$|h_{22}|/\nu$',"Interpreter","latex",'FontSize',18)
    legend('Interpreter','latex','BackgroundAlpha',.7,'FontSize',14)
    ax1.YAxis.Exponent = -3;
    xlims = xlim;
    ylims = ylim;
    xlength = abs(xlims(2)-xlims(1));
    ylength = abs(ylims(2)-ylims(1));
    %text(xlims(1) + 0.1*xlength,ylims(2)-0.2*ylength,sprintf('{\\tt %s}',label),'Interpreter','latex','FontSize',16);

    nexttile
    ax2 = gca;
    set(ax2,'XMinorTick','on','YMinorTick','on','FontSize',10,'FontName','Times','box','on');
    hold on
    plot(t_low+dt,abs(Phase_lm),'DisplayName',sprintf('$|N_%s=%d - N_%s=%d|$',th,n_low,th,n_mid),'LineWidth',linewidth,'Color',MyColors('b1'))
    plot(t_mid+dt,abs(Phase_hm),'DisplayName',sprintf('$|N_%s=%d - N_%s=%d|$',th,n_mid,th,n_high),'LineWidth',linewidth,'Color',MyColors('r1'))
    plot(t_mid+dt,scaling^rate.*abs(Phase_hm),'DisplayName',sprintf('$%g^{%d}|N_%s=%d - N_%s=%d|$',scaling,rate,th,n_mid,th,n_high),'LineStyle','--','LineWidth',linewidth,'Color',MyColors('g1'))
    xlim([x0,x1])
    ylabel('$\phi_{22}/\nu$',"Interpreter","latex",'FontSize',18)
    xlabel('$u$',"Interpreter","latex",'FontSize',18)
    ax2.YAxis.Exponent = -2;

    A4Width(fig);
    fig.Position(3:4) = [21 12];

    if autosave==1
        exportgraphics(fig, sprintf('~/repos/teobiresumsprecessing/Latex/Thesis/Figures/%s_conv.pdf',label), 'ContentType', 'vector')
    end

return

function wf = DB_get_wf(dir,l,m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ouput is structure contaiing waveform mode (l,m), time and dyn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = abs(m);
    if m<0
        fulldir = sprintf('%s-%d/out0d/h_Yl%dm-%d_x10.0000.dat',dir,n,l,n);
    else
        fulldir = sprintf('%s+%d/out0d/h_Yl%dm%d_x10.0000.dat',dir,n,l,n);
    end

    data = readmatrix(fulldir);
    t = data(:,1);
    h = data(:,2) + 1i.*data(:,3);

    wf.h = h;
    wf.t = t;
    wf.l = l;
    wf.m = m;
return
