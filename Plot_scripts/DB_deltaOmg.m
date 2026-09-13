function [Deltas,iotas]  = DB_deltaOmg
    my_linewidth = 1.5;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    inputDB.verbose = 0;

    angles = flip(1:1:90);
    %angles = angles(angles~=0);
    iotas1 = 90-angles;
    iotas2 = 90+flip(angles);
    iotas = [iotas1 iotas2];

    %spins = [0.2 0.4 0.5 0.7 0.9];
    spins = [0.9 0.9 0.9 0.9 0.9];


    n = length(iotas);
    m = length(spins);

    Deltas = zeros(n,m);

    for j = 1:m
        a = spins(j);
        for i = 1:n
            %inputDB.verbose = 1;
            fprintf('%1.f/%1.f   (a = %.1f, iota = %1.f)\n',i+(j-1)*n,n*m,a,iotas(i))
            inputDB.chi1(3) = a.*sign(90-iotas(i));
            inputDB.th0 = deg2rad(abs(90 - iotas(i)));
            inputDB.r0 = DB_LSSO(a,iotas(i)) + 0.5 + 0.05*(j-1);
            DB = DB_class(inputDB);
            Deltas(i,j) = DB_compare_freq(DB,'noplot');   
        end
        clear DB
    end

    C = autumn(m);
    figure
    set(gca,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    hold on
    for i=1:m
        plot(iotas,Deltas(:,i),'Marker','.','LineWidth',my_linewidth,'LineStyle','-','Color',C(i,:),'DisplayName',sprintf('$a = %.1f$ {\\tt n%d}',spins(i),i));
    end
    legend('Location','best','Interpreter','latex','FontSize',legend_fontsize,'BackgroundAlpha',0.7  )
    xlabel('$\iota$','FontSize',labels_fontsize,'Interpreter','Latex');
    ylabel('$t_{LR} - t_{\Omega_{\rm orb}^{\rm max}}$','FontSize',labels_fontsize,'Interpreter','Latex');


return
