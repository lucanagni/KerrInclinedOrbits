function r_LSO = DB_LSSO(a,iota,varargin)
    % For zero eccentricity works for all inclinations and spins
    % For non-zero eccentricity needs to be checked
    
    if ~isempty(varargin)
        e = varargin{1};
    else
        e = 0;
    end

    x = cos(deg2rad(iota));

    if a==0
        r_LSO=6;
        return
    end
    
    
    %while S(a,x,sep)*S(a,x,9)>0
    %    sep = sep - 1e-3;
    %    if sep<5
    %        error('Unable to locate sign change')
    %        %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
    %        r = linspace(4,7,10000);
    %        [~,idx] = min(S(a,x,r));
    %        r_LSO = r(idx);
    %        return
    %    end
    %end

    %fprintf('sep = %f\n',sep)
    if a*x<0 
        sep = 6;
        while S(e,a,x,sep)*S(e,a,x,9)>0
        sep = sep - 1e-3;
            if sep<5
                 error('Unable to locate sign change')
        %        %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
        %        r = linspace(4,7,10000);
        %        [~,idx] = min(S(a,x,r));
        %        r_LSO = r(idx);
        %        return
            end
        end
        r0 = [sep,9];
    elseif a*x>0
        sep = 2;
        while S(e,a,x,sep)*S(e,a,x,1)>0
            sep = sep + 1e-3;
            if sep>8
                error('Unable to locate sign change')
                %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
                r = linspace(1,9,10000);
                [~,idx] = min(S(e,a,x,r));
                r_LSO = r(idx);
                return
            end
        end
        r0 = [1,sep];
    else 
        %r_LSO = 6;
        %return
        r0 = [5,7];
    end



    %r_LSO = fzero(@(r) S(a,x,r),r0);

    %r = linspace(4,7,10000);

    %figure
    %plot(r,log(S(a,x,r)))
    %hold on
    %plot(r,DB_D1(S(a,x,r),r,4))
    %yscale log
    %yline(0)

    r_LSO = fzero(@(r) S(e,a,x,r),r0);


return

function S = S(e,a,x,r)

    %e = 0;

    S12 = 1;
    S11 = -4.*(3 + e);
    S10 = 4.*(3 + e).^2 + 2.*a.^2.*(3 + 2.*e + 3.*e.^2 - 2.*(3 + e.*(2 + e)) .* x.^2);
    S9  = 4 .* a.^2 .* (-7 + e .* (-7 + e .* (-13 - 5.*e + 4 .* (3 + e) .* x.^2)));
    S8  = -16 .* a.^2 .* (-1 + e) .* (1 + e).^2 .* (3 + e) .* (-1 + x.^2) + ...
        a.^4 .* (15 + 20.*e + 26.*e.^2 + 20.*e.^3 + 15.*e.^4 - ...
        4 .* (9 + e .* (12 + e .* (18 + e .* (12 + 5.*e)))) .* x.^2 + ...
        2 .* (15 + e .* (2 + e) .* (10 + 3.*e .* (2 + e))) .* x.^4);
    S7  = -8 .* a.^4 .* (1 + e).^2 .* (-1 + x) .* (1 + x) .* (-3 + e - e.^2 - 5 .* e.^3 + ...
        (15 + e .* (-5 + 3.*e .* (1 + e))) .* x.^2);
    S6  = -4 .* a.^4 .* (1 + e).^2 .* (-1 + x) .* (1 + x) .* ...
        (-2 .* (11 - 14.*e.^2 + 3.*e.^4) .* (-1 + x.^2) + ...
        a.^2 .* (5 + 6.*e.^2 + 5.*e.^4 - (5 + e.^2 .* (6 + e .* (8 + 5.*e))) .* x.^2 + ...
        (-1 + e) .* (3 + e) .* (3 + e .* (2 + e)) .* x.^4));
    S5  = 8 .* a.^6 .* (-1 + e) .* (1 + e).^3 .* (-1 + x.^2).^2 .* (3 + e + e.^2 - 5 .* e.^3 + ...
        2 .* (6 + e .* (2 + e + e.^2)) .* x.^2);
    S4  = a.^6 .* (1 + e).^4 .* (-1 + x.^2).^2 .* (-16 .* (-3 + e) .* (1 + e) .* (-1 + x.^2)) + ...
        a.^2 .* (15 + e .* (-20 + e .* (26 + 5.*e .* (-4 + 3.*e)))) + ...
        6 .* x.^2 - 2 .* e .* (2 + e) .* (2 + e .* (-6 + 5.*e)) .* x.^2 + ...
        (-1 + e).^2 .* (3 + e).^2 .* x.^4;
    S3  = -4 .* a.^8 .* (-1 + e) .* (1 + e).^5 .* (-1 + x.^2).^3 .* ...
        (7 - 7.*e + 13.*e.^2 - 5.*e.^3 + (-1 + e) .* (7 + e.^2) .* x.^2);
    S2  = 2 .* a.^8 .* (-1 + e).^2 .* (1 + e).^6 .* (-1 + x.^2).^3 .* ...
        (2 .* (-3 + e).^2 .* (-1 + x.^2) + a.^2 .* (-3 + 2.*e - 3.*e.^2 + ...
        (-1 + e) .* (3 + e) .* x.^2));
    S1  = -4 .* a.^10 .* (-3 + e) .* (-1 + e).^3 .* (1 + e).^7 .* (-1 + x.^2).^4;
    S0  = a.^12 .* (-1 + e).^4 .* (1 + e).^8 .* (-1 + x.^2).^4;
    
    S_vec = [S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12];
    S = 0;
    for i=1:13
        S = S + S_vec(i).*r.^(i-1);
    end

return

function S = S2(a,s,r)
    x = sqrt(1-s.^2);
    e = 0;

    S12 = 1;
    S11 = -4 .* (3 + e);
    S10 = 4 .* (3 + e).^2 + 2 .* a.^2 .* (-3 - 2.*e + e.^2 + 2 .* (3 + e .* (2 + e)).* x.^2);
    S9 = -4 .* a.^2 .* (7 + e .* (7 + e + e.^2 + 4 .* e .* (3 + e).* x.^2)) ;
    S8 = 16 .* a.^2 .* (-1 + e).*(1 + e).^2 .* (3 + e).*x.^2 + ...
        a.^4 .* ((-3 + e).^2 .* (1 + e).^2 + 2 .* (1 + e) .* (-15 + e .* (-5 + e .* (-1 + 5 .* e))).* x.^2 + ...
        4 .* (3 + e .* (2 + e)).^2 .* x.^4);
    S7 = 8 .* a.^4 .* (1 + e).^2 .* x.^2 .* (15 - 5 .* e + e.^2 - 3 .* e.^3 - 2 .* (9 + e .* (-3 + e + e.^2)) .* x.^2);
    S6 = 4 .* a.^4 .* (1 + e).^2 .* x.^2 .* (2 .* (-1 + e.^2) .* (7 + e.^2 + 2 .* (-9 + e.^2).*x.^2) + ...
        a.^2 .* (-4 + 9 .* x.^2 + e.^2 .* (-2 + 2 .* (-2 + e) .* e + (8 + e .* (4 + 3 .* e)) .* x.^2)));
    S5 = -8 .* a.^6 .* (-1 + e).*(1 + e).^3 .* x.^2 .* (3 + e - e.^2 + e.^3 + 2 .* (- 3 - e + 2 .* e.^3).* x.^2);
    S4 = a.^6 .* (1 + e).^4 .* x.^2 .* (16 .* (-3 + e) .* (-1 + e).^2 .* (1 + e).*x.^2 + ...
        a.^2 .* (2 .* (-3 + e) .* (-1 + e).^2 .* (1 + e) + ...
        (21 + e .* (-28 + e .* (22 + e .* (-12 + 13 .* e)))) .* x.^2)); 
    S3 = -4 .* a.^8 .* (-1 + e) .* (1 + e).^5 .* (-7 + e .* (7 + e .* (-13 + 5 .* e))) .* x.^4;
    S2 = 2 .* a.^8 .* (-1 + e).^2 .* (1 + e).^6 .* (2 .* (-3 + e).^2 + a.^2 .* (3 + e .* (-2 + 3 .* e))) .* x.^4;
    S1 = -4.*a.^10 .* (-3 + e) .* (-1 + e).^3 .* (1 + e).^7 .* x.^4;
    S0 = a.^12 .* (-1 + e).^4 .* (1 + e).^8 .* x.^4;

    S_vec = [S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12];
    S = 0;
    for i=1:13
        S = S + S_vec(i).*r.^(i-1);
    end

return