function r_LR = DB_LR(a,iota)

    x = cos(deg2rad(iota));

    SLR = 3;
    if a==0
        r_LR=SLR;
        return
    end
    
    while LR_pos(a,x,SLR)*LR_pos(a,x,4)>0
        SLR = SLR - 1e-3;
        %if SLR<2.5                                         NOT NEEDED: THE SEXTIC CAN ONLY HAVE 0,2,4 ROOTS
        %    error('Unable to locate sign change')
        %    %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
        %    r = linspace(2,4,10000);
        %    [~,idx] = min(S(a,x,r));
        %    r_LR = r(idx);
        %    return
        %end
    end
%
    if a*x<0 
        %while LR_pos(a,x,SLR)*LR_pos(a,x,4)>0
        %SLR = SLR - 1e-3;
        %    if SLR<1.5
        %        %error('Unable to locate sign change')
        %        %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
        %        r = linspace(2,4,10000);
        %        [~,idx] = min(LR_pos(a,x,r));
        %        r_LR = r(idx);
        %        return
        %    end
        %end
        r0 = [SLR,4];
    elseif a*x>0 
        %while LR_pos(a,x,1)*LR_pos(a,x,SLR)>0
        %SLR = SLR + 1e-3;
        %    if SLR>4
        %        %error('Unable to locate sign change')
        %        %r_LSO = fzero(@(r) DB_D1(S(a,x,r),r,4),6);
        %        r = linspace(2,4,10000);
        %        [~,idx] = min(LR_pos(a,x,r));
        %        r_LR = r(idx);
        %        return
        %    end
        %end
        r0 = [1.01,SLR];
    else 
        %r_LSO = 6;
        %return
        r0 = [2,4];
    end

    %r0
    %LR_pos(a,x,1.01)
    %LR_pos(a,x,r0(2))
    
    r_LR = fzero(@(r) LR_pos(a,x,r),r0);

    equat = 2.*(1+cos((2./3).*acos(-a)));
    %fprintf('computed: %.4f, equatorial: %.4f\n',r_LR,equat)

return

function LR_zero = LR_pos_old(a,x,r)

    l = (1-x.^2);
    LR_zero = r.^6 - 6.*r.^5 + (9 + 2.*a.^2.*l).*r.^4 - 4.*a.^2.*r.^3 + (a.^2 - 6).*a.^2.*l.*r.^2 + 2.*a.^4.*l.*r + a.^4.*l;

return

function LR_zero = LR_pos_sqrt(a,x,r)

    %l = (1-x.^2);
    LR_zero = r.^5 - 3.*r.^4 + 2.*a.^2.*r.^3.*(1 - x.^2) - 2.*a.^2.*r.^2 + ...
     a.^4.*r.*(1 - x.^2) + a.^4.*(1 - x.^2) + ...
     2.*a.*r.*x.*sqrt(3.*r.^4 + (1 - 3.*(1 - x.^2)).*a.^2.*r.^2 - a.^4.*(1 - x.^2));

return

function LR_zero = LR_pos(a,x,r)
    % x = cos(theta_inc)
    %cos^2(theta_inc) = sin^2(tehta_min)
    x2 = x.^2;
    s2 = 1-x2;

    Delta = r.^2 - 2.*r + a.^2;
    lz = 1./(a.*(r-1)).*((r.^2-a.^2) - r.*Delta);
    Q = r.^3./(a.^2.*(r-1).^2).*(4.*Delta - r.*(r-1).^2);

    LR_zero = Q - s2.*(-a.^2+lz.^2./x2);

return