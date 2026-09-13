function t_LR = KOS_tLR_splined(dyn)
    addpath '/home/luca/repos/kerrorbitsolver/KerrDynamics'
    a = dyn.a;

    rs = 2;
    r_LR = rs.*(1+cos((2./3).*acos(-a)));

    R = dyn.r;
    t = dyn.time;

    LR_pos = find(R<r_LR,1);

    idx0 = find(t>=t(LR_pos)-10, 1, 'first');

    R_short = R(idx0:end);
    t_short = t(idx0:end);
    newT = (t(LR_pos)-2):1e-6:(t(LR_pos)+2);

    sR = spline(t_short, R_short, newT);
    [~, idx_LR] = find(sR<r_LR,1);
    t_LR = newT(idx_LR);
        
return