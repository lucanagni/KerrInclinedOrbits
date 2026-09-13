function DB_plot_omegas(eq,tilt,twist)
    a_eq = DB_testmass_checks(eq);
    a_tilt = DB_testmass_checks(tilt);
    a_twist = DB_testmass_checks(twist);

    max_eq_index = find(eq.Omg == max(eq.Omg));
    max_tilt_index = find(tilt.Omg == max(tilt.Omg));
    max_twist_index = find(twist.Omg == max(twist.Omg));

    tmax_eq = eq.t(max_eq_index);
    tmax_tilt = tilt.t(max_tilt_index);
    tmax_twist = twist.t(max_twist_index);

    delta_tilt = tmax_eq - tmax_tilt;
    delta_twist = tmax_eq - tmax_twist;

    xmin = tmax_eq - 40;
    xmax = tmax_eq + 40;
    ymin = eq.Omg(max_eq_index)*(1-0.25);
    ymax = eq.Omg(max_eq_index)*(1+0.05);

    Omg_tilt_interp = spline(tilt.t+delta_tilt,tilt.Omg,eq.t);
    Omg_twist_interp = spline(twist.t+delta_twist,twist.Omg,eq.t);

    eq_name = sprintf('%s_{EQ}^{%.1f}','\Omega',a_eq);
    tilt_name = sprintf('%s_{NEQ}^{%.1f}','\Omega',a_tilt);
    twist_name = sprintf('%s_{NEQ}^{%.4f}','\Omega',a_twist);
    
    figure
    tiledlayout(2,1)

    nexttile 
    hold on
    xline(tmax_eq,'HandleVisibility','off',"Alpha",0.4,"LineStyle","--")
    plot(eq.t,eq.Omg,'DisplayName',eq_name)
    plot(tilt.t + delta_tilt,tilt.Omg,'DisplayName',tilt_name)
    plot(twist.t + delta_twist,twist.Omg,'DisplayName',twist_name)
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    %text(tmax_eq,0.5*(ymax-ymin),'t_{max}')
    title('Orbital Frequencies')
    legend("Location","northwest")
    fontsize(legend,12,'points')

    nexttile
    hold on 
    xline(tmax_eq,'HandleVisibility','off',"Alpha",0.4,"LineStyle","--")
    plot(eq.t,abs(eq.Omg-Omg_tilt_interp)./eq.Omg,"DisplayName",sprintf('%s - %s',eq_name,tilt_name))
    plot(eq.t,abs(eq.Omg-Omg_twist_interp)./eq.Omg,"DisplayName",sprintf('%s - %s',eq_name,twist_name))
    xlim([xmin,xmax])
    grid on
    title(sprintf('Relative Differences (normalized to %s)',eq_name))
    legend("Location","northwest")
    fontsize(legend,12,'points')


return