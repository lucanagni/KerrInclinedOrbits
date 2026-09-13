function [h,t] = DB_sum_lmodes(wavedir,phi,iota,ellmax,m,varargin)

    mmin=0;
    tilt = 90;
    obj = 'h';
    type = 'none';
    orbit = 'eq';
    r0=7;


    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'type'
                    i      = i + 1;
                    type   = char(varargin{i});
                case 'tilt'
                    i      = i + 1;
                    tilt   = varargin{i};
                case 'orbit'
                    i     = i + 1;
                    orbit = varargin{i};
                case 'r0'
                    i     = i + 1;
                    r0 = varargin{i};
                case 'obj'
                    i     = i + 1;
                    obj = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    

    dirbase = sprintf('%s/teuk_HH10_geod_a0.0000_r0%d.000_th0%.3f_q1e+03_3601x161_m',wavedir,r0,tilt);

    if strcmp(orbit,'neq') %For off-equatorial orbits we read in both positive and negative modes

        for l=2:ellmax
            dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,m,l,m);
            wave = readmatrix(dir);
            s.ell(l).emm(m+1).hlm = wave(:,2)+1i.*wave(:,3);
            s.ell(l).emm(m+1).t = wave(:,1);
        end
        for l=2:ellmax
            n=m;
            dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,-m,l,-m);
            wave = readmatrix(dir);
            s.ell(l).emminus(n+1).hlm = wave(:,2)+1i.*wave(:,3);
            s.ell(l).emminus(n+1).t = wave(:,1);
        end

        for l=2:ellmax
            w.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm.*fix_mode(type,l,m);
            n=m;
            w.ell(l).emminus(n+1).hlm = s.ell(l).emminus(n+1).hlm.*fix_mode(type,l,-m);
        end
    
        t  = s.ell(2).emm(3).t;

        % sum over the multipoles
        sumr = 0;
        sumi = 0;
        for l=2:ellmax
            hlm = w.ell(l).emm(m+1).hlm;
            philm  = unwrap(angle(hlm));  %originally had a - (hence the sign of sp below was different)
            Alm    = abs(hlm);
            
            cp     = cos(philm);
            sp     = sin(philm);
            Ylm    = spinsphericalharm(-2, l, m, phi, iota);
            Y_real = real(Ylm);
            Y_imag = imag(Ylm);
            
            
            sumr   = sumr + Alm.*(cp.*Y_real - sp.*Y_imag); 
            sumi   = sumi + Alm.*(-cp.*Y_imag - sp.*Y_real); 
                
            % add m<0 modes
            n=m;
            hlm = w.ell(l).emminus(n+1).hlm;
            philm  = unwrap(angle(hlm)); %originally had a - (hence the sign of sp below was different)
            Alm    = abs(hlm);
            
            cp     = cos(philm);
            sp     = sin(philm);
            Ylm    = spinsphericalharm(-2, l, -m, phi, iota);
            Y_real = real(Ylm);
            Y_imag = imag(Ylm);
            
            sumr   = sumr + Alm.*(cp.*Y_real - sp.*Y_imag);
            sumi   = sumi + Alm.*(-cp.*Y_imag - sp.*Y_real);
        end
        h = sumr - 1i*sumi;
    end
return

function Y = spinsphericalharm(s, l, m, phi, i)
    if l<0 || m<-l || m>l
        error("wrong (l,m) inside spinsphericalharm(s,l,m,phi,i)");
    end
    c       = (-1).^(s)*sqrt( (2.*l+1.)./(4.*pi) );
    %c=1;
    dWigner = c.*wigner_d_function(l,m,-s,i);
    rY      = cos(m*phi).*dWigner;
    iY      = sin(m*phi).*dWigner; %originally was + sign (see https://arxiv.org/pdf/0709.0093 eq (11.7) THAT IS THE CORRECT DEFINITION!
    Y       = rY + 1i*iY;
return
    
function out = wigner_d_function(l,m,s,i)  %as written, this computes d^l_{m,s} (order of indices is important)
    cth  = cos((1+i)*0.5);
    sth  = sin((1-i)*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,s-m);
    kf   = min(l+s,l-m);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+s-k) * factorial(l-k-m) * factorial(k-s+m) );
        dWig = dWig + div*( (-1).^(k-s+m) * cth.^(2*l-m+s-2*k) * sth.^(2*k-s+m) );
    end
    out = dWig*norm;
return

function factor = fix_mode(type,l,m)
    if strcmp(type,'none')
        factor = 1;
    elseif strcmp(type,'eob')
        factor = (-1).^m;
    elseif strcmp(type, 'rwz')
        factor = (-1).^(l);
        if mod(l+m,2) % if odd
            factor = -factor*1i;
        end
    else
        error("type='%s' is not valid. Use 'rwz' or 'eob'", type)
    end
return 
    