function dpth_analytical = DB_dpthdt(obj)
    % TO BE CHECKED
    
    x = obj.x;
    y = obj.y;
    z = obj.z;

    px = obj.px;
    py = obj.py;
    pz = obj.pz;

    [r,phi,th,pr,pphi,pth] = DB_coords_cart2spherical(x,y,z,px,py,pz);

    dHx = obj.dHdx';
    dHp = obj.dHdp';
    st = sin(th);
    ct = cos(th);
    sp = sin(phi);
    cp = cos(phi);
    s2t = sin(th).^2;
    dxdth = r.*ct.*cp;
    dydth = r.*ct.*sp;
    dzdth = -r.*st;
    dpxdth = pr.*ct.*cp + pphi.*(sp.*ct)./(r.*s2t);
    dpydth = pr.*ct.*sp - pphi.*(cp.*ct)./(r.*s2t);
    dpzdth = -pr.*st + pth.*(ct)./(r.^s2t);
    dpth_analytical = -(dHx(:,1).*dxdth + dHx(:,2).*dydth + dHx(:,3).*dzdth + dHp(:,1).*dpxdth + dHp(:,2).*dpydth + dHp(:,3).*dpzdth);

return 