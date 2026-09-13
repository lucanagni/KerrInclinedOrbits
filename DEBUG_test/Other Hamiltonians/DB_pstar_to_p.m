function [px,py,pz] = DB_pstar_to_p(x,y,z,pxs,pys,pzs,chi1)
    n = length(x);
    Px = zeros(n,1);
    Py = zeros(n,1);
    Pz = zeros(n,1);
    for i=1:n
        X = [x(i);y(i);z(i)];
        pstar = [pxs(i);pys(i);pzs(i)];
        [~,Tinv] = DB_Tmatrix(X,chi1);
        P = Tinv*pstar;
        Px(i) = P(1);
        Py(i) = P(2);
        Pz(i) = P(3);
    end
    px = Px;
    py = Py;
    pz = Pz;
end