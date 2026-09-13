function [s,t] = DB_read_modes(wavedir,varargin)
    ellmax = 2;
    tilt = 90;
    mmin   = 0;
    generic = 1;
    wave = 'h';

    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'tilt'
                    i      = i + 1;
                    tilt   = varargin{i};
                case 'generic'
                    i      = i + 1;
                    generic = varargin{i};
                case 'ellmax'
                    i      = i + 1;
                    ellmax = varargin{i};
                case 'wave'
                    i      = i+1;
                    wave = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    dirbase = sprintf('%s/teuk_HH10_geod_a0.0000_r07.000_th0%.3f_q1e+03_3601x161_m',wavedir,tilt);


    for l=2:ellmax
        for m=mmin:l
            dir = sprintf('%s%d/out0d/%s_Yl%dm%d_x10.0000.dat',dirbase,m,wave,l,m);
            wave = readmatrix(dir);
            s.ell(l).emm(m+1).hlm = wave(:,2)+1i.*wave(:,3);
            s.ell(l).emm(m+1).t = wave(:,1);
        end
    end

    t = wave(:,1);
    
    if generic==1
        for l=2:ellmax
            for n=1:l
                m=-n;
                dir = sprintf('%s%d/out0d/%s_Yl%dm%d_x10.0000.dat',dirbase,m,wave,l,m);
                wave = readmatrix(dir);
                s.ell(l).emminus(n+1).hlm = wave(:,2)+1i.*wave(:,3);
                s.ell(l).emminus(n+1).t = wave(:,1);
            end
        end
    end

return 