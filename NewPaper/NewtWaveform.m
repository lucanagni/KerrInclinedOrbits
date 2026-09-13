function NewtWaveform(varargin)
    if isempty(varargin)
        load("~/waveforms/K/plunge/th0_30/a02/wf.mat")
    else
        s = varargin{1};
    end
    dyn = s.dyn;

    %idx_plunge = find(dyn.t>tLSSO_splined(dyn),1);
    idx_plunge = length(dyn.t);
    f3 = figure;
    f3.Position(1) = 100;

    % =======================
    % Plot parameters
    % =======================
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    ax = gca;
    plot3(ax,dyn.x(1:idx_plunge),dyn.y(1:idx_plunge),dyn.z(1:idx_plunge),'LineWidth',2,'Color','b')
    hold on
    %plot3(ax,dyn.x(idx_plunge:end),dyn.y(idx_plunge:end),dyn.z(idx_plunge:end),'LineWidth',2,'Color',MyColors('r1'))
    xlabel('$x$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
    ylabel('$y$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
    zlabel('$z$','FontSize',14,'Interpreter','Latex','FontSize',labels_fontsize)
    if abs(max(dyn.z)-min(dyn.z))<1e-6
        zlim([-1,1])
    else 
        axis equal
    end
    set(ax,'XMinorTick','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5);
    grid(ax,'on')
    view(ax,-15,19)

    xl = xlim;
    yl = ylim;
    zl = zlim;
    text(ax,0.9*xl(1),0.9*yl(end),0.9*zl(end),sprintf('{\\tt %s}',KerrLabel(dyn)),'Interpreter','latex','FontSize',14)

return