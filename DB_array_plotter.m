function DB_array_plotter(obj, plots_list)
% You are a wizard, Plotter 

T = obj.t;
x = obj.x;
y = obj.y;
z = obj.z;
r = obj.r;
th = obj.th;
phi = obj.phi;

Heff = obj.Heff;
Omg = obj.Omg;
C = obj.C;
Lz = obj.pphi;

bool_plot = @(plot_name) return_bool_plot(plots_list, plot_name);

t_LSO = find_t_LSO(obj);

if bool_plot('orbit')
    figure('Name','orbit')
    plot3(x,y,z)
    xlabel('$x$','FontSize',15,'Interpreter','Latex')
    ylabel('$y$','FontSize',15,'Interpreter','Latex')
    zlabel('$z$','FontSize',15,'Interpreter','Latex')
    if abs(max(z)-min(z))<1e-6
        zlim([-1,1])
    end
    title('Orbit')

    drawnow
end

if bool_plot('angle-radius')
    figure('Name','angle-radius')
    polarplot(th, r)
    thetaticklabels({})
    drawnow
end

if bool_plot('radius')
    figure('Name','radius')
    plot(T,r)
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$r$','FontSize',15,'Interpreter','Latex')
    %xline(t_LSO,'Color',[.7 .7 .7],'LineStyle','--')
    drawnow
end

if bool_plot('phi')
    figure('Name','phi')
    plot(T,mod(phi,2*pi))
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$\varphi$','FontSize',15,'Interpreter','Latex')
    xline(t_LSO,'Color',[.7 .7 .7],'LineStyle','--')
    drawnow
end

if bool_plot('theta')
    figure('Name','theta')
    plot(T,rad2deg(th))
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$\theta$','FontSize',15,'Interpreter','Latex')
    drawnow
end

if bool_plot('energy')
    figure('Name','Heff')
    plot(T,Heff)
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$H_{eff}$','FontSize',15,'Interpreter','Latex')
    title('Energy')
    drawnow
end

if bool_plot('frequency')
    figure('Name','Omg')
    plot(T,Omg)
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$\Omega$','FontSize',15,'Interpreter','Latex')
    title('Orbital Frequency')
    drawnow
end

if bool_plot('Carter')
    figure('Name','C')
    plot(T,C)
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$C$','FontSize',15,'Interpreter','Latex')
    title('Carters Constant')
    drawnow
end

if bool_plot('Lz')
    figure('Name','Axial Angular Momentum')
    plot(T,Lz)
    xlabel('$t$','FontSize',15,'Interpreter','Latex')
    ylabel('$L_z$','FontSize',15,'Interpreter','Latex')
    title('Axial Angular Momentum')
    drawnow
end

return

function bool = return_bool_plot(plots_list, plot_name)
bool = any(ismember(plots_list, plot_name)) || any(ismember(plots_list, 'all'));
return

function t_LSO = find_t_LSO(obj)
    a = obj.chi1(3);
    t = obj.t;
    r = obj.r;
    iota = rad2deg(pi/2 - obj.th0);

    r_LSO = DB_LSSO(a,iota);

    idx = find(r<r_LSO,1);
    t_LSO = t(idx);
return