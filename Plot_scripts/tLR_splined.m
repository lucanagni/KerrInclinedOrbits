function [t_LR,t2_LR,dt] = tLR_splined(dyn,varargin)
    a = dyn.chi1(3);
    iota = rad2deg(pi/2-dyn.th0);
    %if strcmp(which,'old')
    %    rs = 2;
    %    r_LR = rs.*(1+cos((2./3).*acos(-a)));
    %elseif strcmp(which,'new')
    r_LR = DB_LR(a,iota);
    %end
    R = dyn.r;
    t = dyn.t;

    LR_pos = find(R<r_LR,1);

    idx0 = find(t>=t(LR_pos)-10, 1, 'first');

    R_short = R(idx0:end);
    t_short = t(idx0:end);
    newT = (t(LR_pos)-1):1e-6:(t(LR_pos)+1);

    sR = spline(t_short, R_short, newT);
    [~, idx_LR] = find(sR<r_LR,1);
    t_LR = newT(idx_LR);

    if ~isempty(varargin)
        t2_LR = tLR_splined(varargin{1});
        dt = t2_LR - t_LR;
    else
        t2_LR = 0;
        dt = 0;
    end
        
return