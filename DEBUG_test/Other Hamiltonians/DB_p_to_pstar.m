function [pxs,pys,pzs] = DB_p_to_pstar(dyn)
    x = dyn.x;
    y = dyn.y;
    z = dyn.z;
    px = dyn.px;
    py = dyn.py;
    pz = dyn.pz;
    chi1 = dyn.chi1;

    n = length(x);
    for i=1:n
        X = [x(i);y(i);z(i)];
        p = [px(i);py(i);pz(i)];
        T = DB_Tmatrix(X,chi1);
        pstar = T*p;
        Pxs(i) = pstar(1);
        Pys(i) = pstar(2);
        Pzs(i) = pstar(3);
    end
    pxs = Pxs';
    pys = Pys';
    pzs = Pzs';
end