function test_xPN(dyn)
    th = dyn.th;
    phi = dyn.phi;
    t = dyn.t;

    rOmg = ((sin(th)).^(1/3)-dyn.chi1(3).*(dyn.r).^(-3/2)).^(2/3).*dyn.r;

    dth = DB_D1(th,t,4);
    d2th = DB_D1(dth,t,4);

    dphi = DB_D1(phi,t,4);
    d2phi = DB_D1(dphi,t,4);

    out = 2.*cos(2.*th).*(2.*1i.*dth.*dphi - d2th) + sin(2.*th).*(4.*dth.^2 + dphi.^2 + 1i.*d2phi);

    %figure
    %plot(t,out)

    %f1 = figure;
    %f1.Position(1) = f1.Position(1) - f1.Position(3)-100;
    %plot(t,abs(out))

    %f2 = figure;
    %f1.Position(1) = f1.Position(1) - f1.Position(3)-100;
    %plot(t,real(out))

    %f3 = figure;
    %f3.Position(1) = f3.Position(1) + f3.Position(3)+100;
    %plot(t,imag(out))

    f4 = figure;
    f4.Position(1) = f4.Position(1) + 2.*(f4.Position(3)+100);
    plot(t,rOmg.^3.*dphi.^2)

    pphi = dyn.pphi;
    pth = dyn.pth;
    r = dyn.r;
    A = 1-2./r;
    H = dyn.Heff;
    out = (A./(H.*r.^2)).^2.*((pphi./sin(th)).^2+pth.^2);


    figure
    plot(t,dyn.Omg.^2)

    figure
    plot(t,1./(dyn.r.^3))