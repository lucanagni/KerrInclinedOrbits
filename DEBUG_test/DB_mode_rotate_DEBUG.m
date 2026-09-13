function [h,t] = DB_mode_rotate_DEBUG(wavedir,alpha,beta,l,M,varargin)
    % ==============================================================================================================================
    % NOTE
    % given an orientation of axes and some modes hlm computes in this reference frame, apply a rotation of euler angles alpha,beta
    % this script computes how the modes would look like in this new RF
    % https://arxiv.org/pdf/2005.05338 pag 8
    % ==============================================================================================================================
    % obj can be either 'h' or 'phi'
    % type can be either 'none', 'eob' or 'rwz'
    % orbit can be either 'eq' or 'neq'
    mmin=0;
    tilt = 90;
    obj = 'h';
    type = 'none';
    orbit = 'eq';
    r0 = 7; 


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

    if strcmp(orbit,'neq')
        for m=mmin:l
            dir = sprintf('%s%d/out0d/%s_Yl%dm%d_x10.0000.dat',dirbase,m,obj,l,m);
            wave = readmatrix(dir);
            s.ell(l).emm(m+1).hlm = (wave(:,2)+1i.*wave(:,3));
            s.ell(l).emm(m+1).t = wave(:,1);
        end
        for n=1:l
            m=-n;
            dir = sprintf('%s%d/out0d/h_Yl%dm%d_x10.0000.dat',dirbase,m,l,m);
            wave = readmatrix(dir);
            s.ell(l).emminus(n+1).hlm = (wave(:,2)+1i.*wave(:,3));
            s.ell(l).emminus(n+1).t = wave(:,1);
        end
        t  = s.ell(l).emm(3).t;


        %fix sign convention
        for m=mmin:l
            w.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm*fix_mode(type,l,m);
        end
        for n=1:l
            w.ell(l).emminus(n+1).hlm = s.ell(l).emminus(n+1).hlm*fix_mode(type,l,m);
        end
        
        hre = 0;
        him = 0;
        for m=mmin:l
            hlm    = w.ell(l).emm(m+1).hlm;
            philm  = unwrap(angle(hlm));
            Alm    = abs(hlm);
            
            cp       = cos(philm);
            sp       = sin(philm);

            dWigner  = wigner_d_function(l,m,M,beta);
            d_r      = cos(m*alpha).*dWigner;
            d_i      = -sin(m*alpha).*dWigner;
        
            hre    = hre + Alm.*(d_r.*cp - d_i.*sp);
            him    = him + Alm.*(d_r.*sp + d_i.*cp);
        end
        for n=1:l
            m=-n;
            hlm    = w.ell(l).emminus(n+1).hlm;
            philm  = unwrap(angle(hlm));
            Alm    = abs(hlm);
            
            cp     = cos(philm);
            sp     = sin(philm);
            dWigner  = wigner_d_function(l,m,M,beta);
            d_r      = cos(m*alpha).*dWigner;
            d_i      = -sin(m*alpha).*dWigner;
            
            hre    = hre + Alm.*(d_r.*cp - d_i.*sp);
            him    = him + Alm.*(d_r.*sp + d_i.*cp);
        end
    elseif strcmp(orbit,'eq')
        for m=mmin:l
            dir = sprintf('%s%d/out0d/%s_Yl%dm%d_x10.0000.dat',dirbase,m,obj,l,m);
            wave = readmatrix(dir);
            s.ell(l).emm(m+1).hlm = (wave(:,2)+1i.*wave(:,3));
            s.ell(l).emm(m+1).t = wave(:,1);
        end
        t  = s.ell(l).emm(3).t;


        %fix sign convention
        for m=mmin:l
            w.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm*fix_mode(type,l,m);
        end
        

        hre = 0;
        him = 0;
        for m=mmin:l
            hlm    = w.ell(l).emm(m+1).hlm;
            philm  = unwrap(angle(hlm));
            Alm    = abs(hlm);
            
            cp       = cos(philm);
            sp       = sin(philm);

            dWigner  = wigner_d_function(l,m,M,beta);
            d_r      = cos(m*alpha).*dWigner;
            d_i      = -sin(m*alpha).*dWigner;
        
            hre    = hre + Alm.*(d_r.*cp - d_i.*sp);
            him    = him + Alm.*(d_r.*sp + d_i.*cp);
            if m~=0
                dWigner  = wigner_d_function(l,-m,M,beta);
                d_r      = cos(-m*alpha).*dWigner;
                d_i      = -sin(-m*alpha).*dWigner;
                if mod(l,2) %if l is odd
                    hre    = hre - Alm.*(d_r.*cp + d_i.*sp);
                    him    = him - Alm.*(-d_r.*sp + d_i.*cp);
                else
                    hre    = hre + Alm.*(d_r.*cp + d_i.*sp);
                    him    = him + Alm.*(-d_r.*sp + d_i.*cp);
                end
            end
        end
    end
    h = hre + 1i.*him;
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