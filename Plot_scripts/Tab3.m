function [Dt,Alm,Omglm] = Tab3(wf,l,m)
    dashes  = '--------------------------------------------------';
    s = struct;

    DB = wf.dyn;
    iota = 90 - rad2deg(DB.th0);

    s = DB_mode_rotate_timedep(wf,DB,l,m,s);

    h = s.ell(l).emm(m+1).hlm;
    t = s.ell(l).emm(m+1).t;

    Title = sprintf('h%d%d - iota = %d°',l,m,iota);

    [peak,t_peak,omg_peak] = Apeak(h,t);
    t_LR = tLR_splined(DB);

    norm = sqrt(factorial(l+2)./factorial(l-2))/2;

    Dt = t_peak - t_LR;
    Alm = peak./norm; %NEW%
    Omglm = omg_peak; %NEW%

    x = 4 - floor(log10(abs(Alm))) - 1;
    y = 4 - floor(log10(abs(Omglm))) - 1;

    fprintf('%s\n',dashes)
    fprintf('%s\n',Title)
    fprintf('%s\n',dashes)
    fprintf(['  Dt     Alm       Omglm \n %.2f   %.',num2str(x),'f   %.',num2str(y),'f \n'],Dt,Alm,Omglm)
    fprintf('%s\n',dashes)

return