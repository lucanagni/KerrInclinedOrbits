function [dEdt,time,dEdt_lm] = DB_dEdt(s,varargin)
    ellmax = 4;

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'ellmax'
                    i      = i + 1;
                    ellmax = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    time = s.ell(2).emm(3).t;

    dEdt = zeros(length(time),1);

    for l=2:ellmax
        for m=-l:l
            %fprintf('l=%1f,m=%1f\n',l,m)
            n = abs(m);
            if m<0
                Dhlm = spline(s.ell(l).emminus(n+1).t,s.ell(l).emminus(n+1).Dhlm,time);
                %Dhlm = s.ell(l).emminus(n+1).Dhlm;
                dEdt_lm.ell(l).emminus(n+1).flux = abs(Dhlm).^2/16/pi;
            else
                Dhlm = spline(s.ell(l).emm(m+1).t,s.ell(l).emm(m+1).Dhlm,time);
                %Dhlm = s.ell(l).emm(m+1).Dhlm;
                dEdt_lm.ell(l).emm(m+1).flux = abs(Dhlm).^2/16/pi;
            end
            dEdt = dEdt + abs(Dhlm).^2;
        end
        dEdt_lm.ell(l).emminus(1).flux = dEdt_lm.ell(l).emm(1).flux;
    end

    dEdt = dEdt./16./pi;

return