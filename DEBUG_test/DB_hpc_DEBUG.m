function [h,t] = DB_hpc(wavedir,iota,phi,varargin)
    % [h,t]=compute_hpc(wavedir,iota,phi,varargin)
    % wavedir = directory containing Teukode's output folders (same dynamics, different lms)
    %
    % Options for varargin:
    %   - 'ellmax' : max value of l used in the sum, default is 8
    %   - 'mmin'   : min value of m used in the sum, default is 1
    %   - 'type'   : can be 'eob' (default) or 'rwz' and it's needed to fix 
    %                the sign convention of non-EOB waveforms.
    %   - 'tilt'   : initial theta angle of the orbit (i.e. th0. Needed for the names of the directories)
    %   - 'generic': set to 1 if the orbit is non-equatorial. Reads m<0 modes from file
    % ===================================================================================================================================================
    % REFERENCES
    % 1_ Wigner's small d-matrices: https://arxiv.org/pdf/2005.05338
    % 2_ SWSH: https://arxiv.org/pdf/0709.0093 (note that (II.7) here agrees with (16) above if one uses the property of the D matrices found in Ref. 1)
    %
    % These formulas agree with the ones found in Goldberg 1967 (except for the factor (-1)^s). Note that in G. the matrices D are D* of Ref. 1 and that 
    % the exponents of the exponentials are different, but the definition of SWSH in terms of D is coherent
    
    ellmax = 3;
    mmin   = 0;
    type   = 'eob';
    tilt = 90.000;
    orbit = 'eq';

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
                case 'ellmax'
                    i      = i + 1;
                    ellmax = varargin{i};
                case 'mmin'
                    i      = i + 1;
                    mmin   = varargin{i};
                case 'tilt'
                    i      = i + 1;
                    tilt   = varargin{i};
                case 'orbit'
                    i      = i + 1;
                    orbit = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    
    dirbase = sprintf('%s/teuk_HH10_geod_a0.0000_r07.000_th0%.3f_q1e+03_3601x161_m',wavedir,tilt);

    if strcmp(orbit,'neq') %For off-equatorial orbits we read in both positive and negative modes

        for l=2:ellmax
            for m=mmin:l
                dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,m,l,m);
                wave = readmatrix(dir);
                s.ell(l).emm(m+1).hlm = wave(:,2)+1i.*wave(:,3);
                s.ell(l).emm(m+1).t = wave(:,1);
            end
        end
        for l=2:ellmax
            for n=1:l
                m=-n;
                dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,m,l,m);
                wave = readmatrix(dir);
                s.ell(l).emminus(n+1).hlm = wave(:,2)+1i.*wave(:,3);
                s.ell(l).emminus(n+1).t = wave(:,1);
            end
        end

        for l=2:ellmax
            for m=mmin:l
                w.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm.*fix_mode(type,l,m);;
            end
            for n=1:l
                m=-n;
                w.ell(l).emminus(n+1).hlm = s.ell(l).emminus(n+1).hlm.*fix_mode(type,l,m);
            end
        end
    
        t  = s.ell(2).emm(3).t;

        % sum over the multipoles
        sumr = 0;
        sumi = 0;
        for l=2:ellmax
            for m=mmin:l
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
            end
            for n=1:l
                m=-n;
                hlm = w.ell(l).emminus(n+1).hlm;
                philm  = unwrap(angle(hlm)); %originally had a - (hence the sign of sp below was different)
                Alm    = abs(hlm);
                
                cp     = cos(philm);
                sp     = sin(philm);
                Ylm    = spinsphericalharm(-2, l, m, phi, iota);
                Y_real = real(Ylm);
                Y_imag = imag(Ylm);
                
                sumr   = sumr + Alm.*(cp.*Y_real - sp.*Y_imag);
                sumi   = sumi + Alm.*(-cp.*Y_imag - sp.*Y_real);
            end
        end

    elseif strcmp(orbit,'eq') %For equatorial orbits we only need positive modes
        for l=2:ellmax
            for m=mmin:l
                dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,m,l,m);
                wave = readmatrix(dir);
                s.ell(l).emm(m+1).hlm = wave(:,2)+1i.*wave(:,3);
                s.ell(l).emm(m+1).t = wave(:,1);
            end
        end
        t  = s.ell(2).emm(3).t;

        for l=2:ellmax
            for m=mmin:l
                w.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm.*fix_mode(type,l,m);;
            end
        end

        % sum over the multipoles
        sumr = 0;
        sumi = 0;
        for l=2:ellmax
            for m=mmin:l
                hlm    = w.ell(l).emm(m+1).hlm;
                philm  = -unwrap(angle(hlm));
                Alm    = abs(hlm);
                
                cp     = cos(philm);
                sp     = sin(philm);
                Ylm    = spinsphericalharm(-2, l, m, phi, iota);
                Y_real = real(Ylm);
                Y_imag = imag(Ylm);
                
                
                sumr   = sumr + Alm.*(cp.*Y_real + sp.*Y_imag);
                sumi   = sumi + Alm.*(cp.*Y_imag - sp.*Y_real);
                
                % add m<0 modes
                if m~=0
                    Ylm_mneg    = spinsphericalharm(-2, l, -m, phi, iota);
                    Y_real_mneg = real(Ylm_mneg);
                    Y_imag_mneg = imag(Ylm_mneg);

                    if mod(l,2)
                        sumr = sumr - Alm.*(cp.*Y_real_mneg - sp.*Y_imag_mneg);
                        sumi = sumi - Alm.*(cp.*Y_imag_mneg + sp.*Y_real_mneg);
                    else
                        sumr = sumr + Alm.*(cp.*Y_real_mneg - sp.*Y_imag_mneg);
                        sumi = sumi + Alm.*(cp.*Y_imag_mneg + sp.*Y_real_mneg);
                    end
                end
            end
        end
    end

    h = sumr - 1i*sumi;
return

% Definition of Ylm
function Y = spinsphericalharm(s, l, m, phi, i)
    if l<0 || m<-l || m>l
        error("wrong (l,m) inside spinsphericalharm(s,l,m,phi,i)");
    end
    c       = (-1).^(-s)*sqrt( (2.*l+1.)./(4.*pi) );
    dWigner = c.*wigner_d_function(l,m,-s,i);
    rY      = cos(m*phi).*dWigner;
    iY      = sin(m*phi).*dWigner;
    Y       = rY + 1i*iY;
return
    
function out = wigner_d_function(l,m,s,i)
    cth  = cos(i*0.5);
    sth  = sin(i*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,m-s);
    kf   = min(l+m,l-s);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+m-k) * factorial(l-s-k) * factorial(s-m+k) );
        dWig = dWig + div*( (-1).^k * cth.^(2*l+m-s-2*k) * sth.^(2*k+s-m) );
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
    
    