function t_LSSO = KOS_tLSSO_splined(dyn)
    addpath '/home/luca/repos/kerrorbitsolver/KerrDynamics'

    a = dyn.a;
    r_LSSO = kerr_separatrix(0,a);
    R = dyn.r;
    t = dyn.time;

    LSSO_pos = find(R<r_LSSO,1);

    idx0 = find(t>=t(LSSO_pos)-10, 1, 'first');

    R_short = R(idx0:end);
    t_short = t(idx0:end);
    newT = (t(LSSO_pos)-2):1e-6:(t(LSSO_pos)+2);

    sR = spline(t_short, R_short, newT);
    idx_LSSO = find(sR<r_LSSO,1);
    t_LSSO = newT(idx_LSSO);    
        
return