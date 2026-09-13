function Deltas = Table_Deltat
    % =====================================================================
    % Produces Table I of Paper I
    % =====================================================================
    spins = 0:0.1:0.9;
    iotas = [0, 30, 60, 85, 95, 120, 150, 180];

    inputDB.verbose = 0;

    m = length(spins);
    n = length(iotas);

    Deltas = zeros(m,n);

    for i=1:m
        a = spins(i);
        for j=1:n
            iota = iotas(j);
            fprintf('%1.f/%1.f   (a = %.1f, iota = %1.f)\n',j+(i-1)*n,n*m,a,iota)
            inputDB.th0 = deg2rad(abs(90 - iotas(j)));
            inputDB.r0 = DB_LSSO(a,iota) + 1;
            if iota<90
                inputDB.chi1(3) = a;
            else
                inputDB.chi1(3) = -a;
            end

            DB = DB_class(inputDB);
            Deltas(i,j) = DB_compare_freq(DB,'noplot');
        end
    end


    fprintf('$a$ & $\\iota = 0$ & $\\iota = 30^\\circ$ & $\\iota = 60^\\circ$ & $\\iota = 85^\\circ$ & $\\iota = 95$ & $\\iota = 120^\\circ$ & $\\iota = 150^\\circ$ & $\\iota = 180^\\circ$ \\\\ \\hline\n')


    for i=1:m
        a = spins(i);
        fprintf('$%.1f$ ',a)
        for j=1:n
            iota = iotas(j);
            fprintf('& $%.2f$ ',Deltas(i,j))
        end
        fprintf('%s \n','\\')
    end


return