function [t_LSSO,t2_LSSO,dt] = tLSSO_splined(dyn,varargin)
    debug = 0;

    a = dyn.chi1(3);
    iota = rad2deg(pi/2-dyn.th0);
    r_LSSO = DB_LSSO(a,iota);
    %r_LSSO = 5.9865;
    R = dyn.r;
    t = dyn.t;

    LSSO_pos = find(R<r_LSSO,1);
    if debug==1
        % =================================================================================================================
        % DEBUG
        % =================================================================================================================
        load('/home/luca/waveforms/S/fluxes/th0_90/wf.mat');
        %r = struct;
        %r = DB_mode_rotate_timedep(s,s.dyn,2,2,r);
        %clear s;
        %s = r;
        H = s.ell(2).emm(3).hlm;
        T = s.ell(2).emm(3).t;
        f22 = freq(H,T);
        idx0_D = find(T>400,1);
        f22_cut = f22(idx0_D:end);
        t_cut = T(idx0_D:end);

        idx1_D = find(f22_cut>0.167,1);
        t_freq = t_cut(idx1_D);

        f_short = f22_cut(idx1_D-100:idx1_D+100);
        t_short = t_cut(idx1_D-100:idx1_D+100);
        newT = t_short(1):1e-5:t_short(end);
        f_spl = spline(t_short,f_short,newT);
        idx2_D = find(f_spl>0.167,1);
        t_freq_spl = newT(idx2_D);
        t_LSSO = t_freq_spl;

        idx_dyn = find(dyn.t>t_LSSO,1);
        dyn.r(idx_dyn);
        % =================================================================================================================
        return
    end

    idx0 = find(t>=t(LSSO_pos)-10, 1, 'first');

    R_short = R(idx0:end);
    t_short = t(idx0:end);
    newT = (t(LSSO_pos)-2):1e-6:(t(LSSO_pos)+2);

    sR = spline(t_short, R_short, newT);
    idx_LSSO = find(sR<r_LSSO,1);
    t_LSSO = newT(idx_LSSO);

    if ~isempty(varargin)
        t2_LSSO = tLSSO_splined(varargin{1});
        dt = t2_LSSO - t_LSSO;
    else
        t2_LSSO = 0;
        dt = 0;
    end
        
return