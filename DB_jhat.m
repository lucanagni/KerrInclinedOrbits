function jhat = DB_jhat(dyn,flags)
    r = dyn.r;
    t = dyn.t;

    a = dyn.chi1(3);
    
    x = dyn.x;
    y = dyn.y;
    z = dyn.z;

    vx = DB_D1(x,t,4);
    vy = DB_D1(y,t,4);
    vz = DB_D1(z,t,4);

    px = dyn.px;
    py = dyn.py;
    pz = dyn.pz;

    jhat = zeros(length(r),1);

    for i=1:length(r)

        R = [x(i) y(i) z(i)];
        V = [vx(i) vy(i) vz(i)];
        P = [px(i) py(i) pz(i)];

        L = cross(R,P);
        Lnorm = norm(L);
        LN = cross(R,V);
        LNnorm = norm(LN);

        jhat(i) = Lnorm./LNnorm;
    end
    if ~flags.r
        jhat = jhat.*(1+a.*r.^(-3/2)).^(-4/3);
    end
return