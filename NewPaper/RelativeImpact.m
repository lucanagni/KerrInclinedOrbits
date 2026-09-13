function out = RelativeImpact(a,varargin)
    %iotas = [0,pi/6,pi/4,pi/3];
    I = 0:5:89;
    iotas = deg2rad(I);
    if isempty(varargin)
        flags.r=0;
        flags.Newt_switch = 0;
        flags.source = 1;
        %flags.analytical_label = 'EOB';

        

        inputDB.chi1(3) = a;
        inputDB.r0 = 7;
        inputDB.Tmax = 2000;
        inputDB.geodesics = 1;
        inputDB.verbose = 0;

        U1 = zeros(length(iotas),1);
        V1 = zeros(length(iotas),1);
        U2 = zeros(length(iotas),1);
        V2 = zeros(length(iotas),1);
        
        for i=1:length(iotas)
            fprintf('Computing configuration %d/%d\n',i,length(iotas))
            inputDB.th0 = pi/2-iotas(i);
            DB = DB_class(inputDB);
            [Ulm,Vlm] = DB_Multipoles(DB,2,1,flags);
            U1(i,:) = mean(Ulm);
            V1(i,:) = mean(Vlm);

            [Ulm,Vlm] = DB_Multipoles(DB,2,2,flags);
            U2(i,:) = mean(Ulm);
            V2(i,:) = mean(Vlm);

            out = [U1,V1,U2,V2];
        end
    else 
        in = varargin{1};
        U1 = in(:,1);
        V1 = in(:,2);
        U2 = in(:,3);
        V2 = in(:,4);
    end
    % ==============
    % Plot Parameters
    % ==============
    my_linewidth = 1;
    axes_fontsize = 12;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 20;

    %t_in = 0;U1 = in(:,1);
    %t_end = 2000;

    f = figure;
    f.Position(3) = 2*f.Position(3);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    ax1 = nexttile;
    set(ax1,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    %legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
    %xlim([t_in,t_end])
    hold on

    scatter(rad2deg(iotas),abs(U1)./abs(V1),'Marker','x')
    ylabel('$U_{21}/V_{21}$','FontSize',labels_fontsize,'FontName','Times','Interpreter','Latex')
    xlabel('$\iota\, (^\circ)$','FontSize',labels_fontsize,'FontName','Times','Interpreter','Latex')
    grid on
    yline(1,'Color',[1 1 1])

    ax2 = nexttile;
    set(ax2,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    %legend('Location','northeast','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7,'NumColumns',2)%,'Orientation','horizontal')
    %xlim([t_in,t_end])
    hold on

    scatter(rad2deg(iotas),abs(V2)./abs(U2),'Marker','v')
    ylabel('$V_{22}/U_{22}$','FontSize',labels_fontsize,'FontName','Times','Interpreter','Latex')
    xlabel('$\iota\, (^\circ)$','FontSize',labels_fontsize,'FontName','Times','Interpreter','Latex')
    grid on
    yline(1,'Color',[1 1 1])
return