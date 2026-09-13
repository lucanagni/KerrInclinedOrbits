function [h,t] = DB_hpc(s,iota,phi,varargin)
    % [h,t]=compute_hpc(wavedir,iota,phi,varargin)
    % wavedir = directory containing Teukode's output folders (same dynamics, different lms)
    %
    % Options for varargin:
    %   - 'ellmax' : max value of l used in the sum, default is 8
    %   - 'tilt'   : initial theta angle of the orbit (i.e. th0. Needed for the names of the directories)
    %   - 'generic': set to 1 if the orbit is non-equatorial. Reads m<0 modes from file
    % ===================================================================================================================================================
    
    ellmax = 4;
    norm = 'teukode';

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
                case 'norm'
                    i      = i + 1;
                    norm = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end
    
    if strcmp(norm,'RWZ')
        s = normalize(s,ellmax);
    end
    t  = s.ell(2).emm(3).t;

    % sum over the multipoles
    h = 0;

    for l=2:ellmax
        for m=0:l
            hlm = s.ell(l).emm(m+1).hlm;
            Ylm = spinsphericalharm(-2, l, m, phi, iota);

            if length(h)>1
                if length(hlm)>length(h)
                    hlm = hlm(1:length(h));
                    t = t(1:length(h));
                elseif length(h)>length(hlm)
                    h = h(1:length(hlm));
                    t = t(1:length(hlm));
                end
            end


            h   = h + hlm.*Ylm;
        end
        for n=1:l
            m   = -n;
            hlm = s.ell(l).emminus(n+1).hlm;
            Ylm = spinsphericalharm(-2, l, m, phi, iota);
            if length(h)>1
                if length(hlm)>length(h)
                    hlm = hlm(1:length(h));
                    t = t(1:length(h));
                elseif length(h)>length(hlm)
                    h = h(1:length(hlm));
                    t = t(1:length(hlm));
                end
            end

            h   = h + hlm.*Ylm;
        end
    end
return

% Definition of Ylm
function Y = spinsphericalharm(s, l, m, phi, i)
    if l<0 || m<-l || m>l
        error("wrong (l,m) inside spinsphericalharm(s,l,m,phi,i)");
    end
    c       = (-1).^(s)*sqrt( (2.*l+1.)./(4.*pi) );
    dWigner = wigner_d_function(l,-s,m,i);
    DWigner = dWigner.*exp(1i*m*phi);
    
    Y       = c.*DWigner;
return
    
function out = wigner_d_function(l,m,s,beta)  %as written, this computes d^l_{m,s} (order of indices is important)
    cth  = cos(beta*0.5);
    sth  = sin(beta*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,s-m);
    kf   = min(l+s,l-m);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+s-k) * factorial(l-k-m) * factorial(k-s+m) );
        dWig = dWig + div.*( (-1).^(k) .* cth.^(2*l-m+s-2*k) .* sth.^(2*k-s+m) );
    end
    out = dWig*norm;
return
    
function s = normalize(s,ellmax)
    for l=2:ellmax
        norm = sqrt(factorial(l+2)./factorial(l-2))/2;
        for m=-l:l
            if m<0
                s.ell(l).emminus(-m+1).hlm = s.ell(l).emminus(-m+1).hlm./norm;
            else
                s.ell(l).emm(m+1).hlm = s.ell(l).emm(m+1).hlm./norm;
            end
        end
    end

return