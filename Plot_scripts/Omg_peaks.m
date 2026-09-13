function Omg_peaks
    my_markersize = 4;
    axes_fontsize = 10;
    legend_fontsize = axes_fontsize+2;
    labels_fontsize = 18;

    inputDB.verbose = 0;

    %angles = flip(30:20:90);
    angles = flip(5:5:90);
    iotas1 = 90-angles;
    iotas2 = 90+flip(angles);
    iotas = [iotas1 iotas2];

    spins = 0:0.1:0.9;

    n = length(iotas);
    m = length(spins);

    fprintf('Total Configurations: %d\n',n*m)
    pause

    figure
    ax = gca;
    set(ax,'XMinorTick','on','box','on','YMinorTick','on','FontSize',axes_fontsize,'FontName','Times','LineWidth',0.5,'box','on');
    xlim([-5,185])
    ylim([0-.1,1])
    xlabel('$\iota \ (^\circ)$','Interpreter','latex','FontSize',20)
    ylabel('$a$','Interpreter','latex','FontSize',20)
    hold on

    Deltas = zeros(n,m);

    for j = 1:m
        a = spins(j);
        for i = 1:n
            iota = iotas(i);
            %inputDB.verbose = 1;
            fprintf('%1.f/%1.f   (a = %.1f, iota = %1.f)\n',i+(j-1)*n,n*m,a,iota)
            inputDB.chi1(3) = a.*sign(90-iota);
            inputDB.th0 = deg2rad(abs(90 - iota));
            inputDB.r0 = DB_LSSO(a,iota) + 0.48;
            DB = DB_class(inputDB);
            %Deltas(i,j) = DB_compare_freq(DB,'noplot');   

            point = binary_peak(DB);
            %DB_compare_freq(DB)
            if point
                plot(ax,iota,a,'MarkerSize',my_markersize,'Marker','.','Color','k')
            end
        end
        clear DB
    end

return

function p = binary_peak(dyn)
    tLR = tLR_splined(dyn);
    range = [tLR-10,tLR+10];
    
    [~,t_peak] = findpeaks(dyn.Omg,dyn.t);
    t_peak = t_peak(t_peak>range(1));
    t_peak = t_peak(t_peak<range(2));
    if length(t_peak)>1
        error('More than one peak within chosen interval: n = %d',length(t_peak))
    end
    %while length(t_peak)>1
    %    range = [tLR-10,tLR+10];
    %    t_peak = t_peak(t_peak>range(1));
    %    t_peak = t_peak(t_peak<range(2));
    %end

    p = length(t_peak);
return