function orbits(Q,P)
    f = figure;
    tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact')

    %figure
    nexttile
    plot3(Q.x,Q.y,Q.z,'LineWidth',0.8)
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',8,'FontName','Times');
    zlim([-5,5])
    %title('Conservative Dynamics')

    %figure
    nexttile
    plot3(P.x,P.y,P.z,'LineWidth',0.8)%,'Color',[0.8500 0.3250 0.0980])
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',8,'FontName','Times');
    %title('Plunging Dynamics')

    hold on
    %[~,Max_locs] = findpeaks(P.Omg);
    %start = Max_locs(end);

    tLSSO = tLSSO_splined(P);
    idx = find(P.t>tLSSO,1);
    plot3(P.x(idx:end),P.y(idx:end),P.z(idx:end),'LineWidth',2)%,'Color',[0.9290 0.6940 0.1250])

    A4Width(f);
return
