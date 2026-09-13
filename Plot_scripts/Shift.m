% Return the shift of the Teukode simulation
% Input: - t0, r0 , a (t0, r0 in BL coordinates)
%        - empty -> find t0, r0 loading the kerr_dynamic
%
function shift = Shift(varargin)

    S   = 10;
    rho = 10;
    M   = 1;

    if length(varargin)==3
        t0 = varargin{1};
        r0 = varargin{2};
        a  = varargin{3};
    elseif isempty(varargin)
        dyn_found = 0;
        files = dir();
        for i=1:length(files)
            fname = files(i).name;
            if startsWith(fname, 'kerr_dyn_a')
                X = load([fname,'/traj.dat']);
                t0 = X(1,1);
                r0 = X(1,2);
                dyn_info = split(fname, '_');
                a = str2double(erase(dyn_info(3), 'a'));
                dyn_found = 1;
                break
            end
        end
        if ~dyn_found
            error('kerr_dyn dir not found!')
        end
    else
        error('Wrong input')
    end
    
    shift_HH = TrajShift(r0, t0, a, M, S);
    shift = -rho -4*M*log( (S*rho+2*M*rho-2*M*S)/S )+2*M*log(2*M)-shift_HH;
     
return

function out = TrajShift(rBL, tBL, a, m ,s)

    if abs(m-a)<1.e-13 % for extremal cases
        rs = rBL -((2*m*(m + (-rBL + m)*log(rBL - m)))/(rBL - m));
    else
        sqrtma = sqrt(abs(m*m-a*a));
        rplus  = m + sqrtma;
        rmins  = m - sqrtma;
        oodr   = 2*m/(rplus-rmins);
        tmp    = ( rplus*log(rBL-rplus) - rmins*log(rBL-rmins) )*oodr;
        rs     = rBL + tmp;
    end
    tki = tBL - rBL + rs;
    R   = s*rBL/(s+rBL);
    out = - tki - 4*m*log(1.-R/s) + R*R/(s-R);
    
    % truncate the output to the 4th decimal, as in Teukode's screen-output
    % out = round(out*1e4)/1e4;
    
return
