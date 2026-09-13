function mode_contribution(s,p,varargin)
    type = 'E';
    start = 'LR';
    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        switch flag
            case 'E'
                type = 'E';
            case 'J'
                type = 'J';
            case 'LR'
                start = 'LR';
            case 'LSSO'
                start = 'LSSO';
            otherwise
                error("'%s' is not a valid flag!", flag)
        end
        i = i + 1;
    end

    [D1,D1lm] = Integrate_flux_new(s,type,start,'multipoles',1);
    [D2,D2lm] = Integrate_flux_new(p,type,start,'multipoles',1);

    fprintf('------------------------------------\n')
    fprintf('( l, m)         old           new \n')
    fprintf('------------------------------------\n')
    for l = 2:4
        for m = flip(-l:l)
            if m < 0
                f1 = D1lm.ell(l).emminus(-m+1).DeltaFlux;
                f2 = D2lm.ell(l).emminus(-m+1).DeltaFlux;
            else
                f1 = D1lm.ell(l).emm(m+1).DeltaFlux;
                f2 = D2lm.ell(l).emm(m+1).DeltaFlux;
            end

            % Format: (l,m) left-aligned in width 7, values right-aligned in width 10 with 2 decimals
            fprintf('(%2d,%2d)   %10.2f%%   %10.2f%%\n', l, m, f1/D1*100, f2/D2*100)
        end
    end
    fprintf('------------------------------------\n')
    fprintf('Delta%s_%s    %7.4f    %10.4f\n',type,start,   D1,D2)

return