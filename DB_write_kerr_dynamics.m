function DB_write_kerr_dynamics(obj)

%==========================================================================
% Write dynamics for Teukode. 
% Structure similar to KerrOrbitSolver
%==========================================================================

q  = obj.q;
nu = q/(1+q)^2;

if strcmp(obj.hamiltonian,'balmelli')
    error('Trying to write Kerr dynamics with Balmelli Hamiltonian')
elseif strcmp(obj.hamiltonian, 'kerr')
    a = DB_testmass_checks(obj,'Trying to write dynamics, but ');

    T = obj.t;

    x = obj.x;
    y = obj.y;
    z = obj.z;

    px = obj.px;
    py = obj.py;
    pz = obj.pz;

    r = obj.r;
    phi = obj.phi;
    th = obj.th;

    pr = obj.pr;
    pphi = obj.pphi;
    pth = obj.pth;

    H = obj.Heff; % in the test-mass limit the Hamiltonian is the effective one 

    prs = pr*0;

    %compute time derivatives analytically
    dHdp = obj.dHdp;
    dHdx = obj.dHdx;
    [dX,dY,dZ,dPx,dPy,dPz] = DB_coords_spherical_derivatives(r,phi,th,pr,pphi,pth);
    
    X = [x y z];
    P = [px py pz];
    L = length(obj.t);
    F = zeros(L,3);
    for i = 1:L
        F(i,:) = DB_flux2(X(i,:),P(i,:),dHdp(i,:),q,obj.chi1,0);
    end
    Fx = F(:,1);
    Fy = F(:,2);
    Fz = F(:,3);
    [~,~,~,Fr,Fphi,Fth] = DB_coords_cart2spherical(x,y,z,Fx,Fy,Fz);

    dphi = dHdp(:,1).*dPx.dpphi + dHdp(:,2).*dPy.dpphi + dHdp(:,3).*dPz.dpphi;

    dpr = -(dHdx(:,1).*dX.dr + dHdx(:,2).*dY.dr + dHdx(:,3).*dZ.dr + dHdp(:,1).*dPx.dr + dHdp(:,2).*dPy.dr + dHdp(:,3).*dPz.dr) + Fr;
    dpth = -(dHdx(:,1).*dX.dth + dHdx(:,2).*dY.dth + dHdx(:,3).*dZ.dth + dHdp(:,1).*dPx.dth + dHdp(:,2).*dPy.dth + dHdp(:,3).*dPz.dth) + Fth;
    dpphi = -(dHdx(:,1).*dX.dphi + dHdx(:,2).*dY.dphi + dHdx(:,3).*dZ.dphi + dHdp(:,1).*dPx.dphi + dHdp(:,2).*dPy.dphi + dHdp(:,3).*dPz.dphi) + Fphi;
    
    %DEBUG -> Check agreement of analytical and numerical value 
    %{
    dphin = DB_D1(phi,T,4);
    dprn = DB_D1(pr,T,4);
    dpthn = DB_D1(pth,T,4);
    dpphin = DB_D1(pphi,T,4);

    anal = [dphi dpr dpth dpphi];
    num = [dphin dprn dpthn dpphin];

    figure
    hold on
    labels = ["$\varphi$" "$p_r$" "$p_\theta$" "$p_\varphi$"];
    for i=1:4
        plot(T,anal(:,i)-num(:,i),'DisplayName',labels(i))
    end
    legend('Interpreter','latex','FontSize',14)
    xlim([0,obj.t(end)-20])

    N=3;
    figure
    hold on
    plot(T,anal(:,N),'DisplayName',sprintf('%s, analytical',labels(N)))
    plot(T,num(:,N),'DisplayName',sprintf('%s, numerical',labels(N)))
    legend('Location','northwest','Interpreter','latex','FontSize',14,'BackgroundAlpha',0.7,'NumColumns',2)
    error('debug')
    %}

    %dphi = DB_D1(phi,T,4);
    d2phi = DB_D1(dphi,T,4);

    dth = DB_D1(th,T,4);
    d2th = DB_D1(dth,T,4);

    dH     = DB_D1(H,T,4);
    d2H    = DB_D1(dH,T,4);

    %dpr    = DB_D1(pr,T,4);
    d2pr   = DB_D1(dpr,T,4);

    %dpth   = DB_D1(pth,T,4);
    d2pth  = DB_D1(dpth,T,4);

    %dpphi  = DB_D1(pphi,T,4);
    d2pphi = DB_D1(dpphi,T,4);

    dprs   = DB_D1(prs,T,4);
    d2prs  = DB_D1(dprs,T,4);

    OneColumn = [ T; r; th; phi; dphi; d2phi; H; dH; d2H; pr; dpr; d2pr; pth; dpth; d2pth; pphi; dpphi; d2pphi; prs; dprs; d2prs;];
    nvars = 21;

    %DEBUG
    %{
    figure
    plot(T,dpr)
    title('$d p_r/dt$','Interpreter','latex')
    figure
    plot(T,d2pr)
    title('$d ^2 p_r/dt^2$','Interpreter','latex')
    figure
    plot(T,dpth)
    title('$d p_\theta/dt$','Interpreter','latex')
    figure
    plot(T,d2pth)
    title('$d ^2 p_\theta/dt^2$','Interpreter','latex')
    figure
    plot(T,dpphi)
    title('$d p_\phi/dt$','Interpreter','latex')
    figure
    plot(T,d2pphi)
    title('$d ^2 p_\phi/dt^2$','Interpreter','latex')

    figure
    plot(T,dphi)
    title('$d \phi/dt$','Interpreter','latex')
    figure
    plot(T,d2phi)
    title('$d ^2 \phi/dt^2$','Interpreter','latex')

    figure
    plot(T,dH)
    title('$d H/dt$','Interpreter','latex')
    figure
    plot(T,d2H)
    title('$d ^2 H/dt^2$','Interpreter','latex')
    %}

    if(length(OneColumn) ~= nvars*length(T) )
        error('!!!!! Missing vars !!!!! exit..');
    end 

    % Definitions for filename, useful for file names  
    r0 = obj.r0;
    th0 = rad2deg(obj.th0);

    ID = sprintf('a%.4f_th0%.3f_r0%.3f_q%.0e',a,th0,r0,q);
    if obj.geodesics == 1
        ID = ['geod_',ID];
    end

    dirname = sprintf('%s/DB_dyn_%s',obj.matlab_outdir_folder,ID);

    if exist(dirname,'dir')
        if obj.verbose
            disp(['Existing dir: ',dirname]);
        end
    else
        system(['mkdir ',dirname]);
        if obj.verbose
            disp(['Creating dir: ',dirname]);
        end
    end

    % file pointers
    fid(1)  = fopen([dirname,'/traj.dat']   ,'w');
    fid(2)  = fopen([dirname,'/pr.dat']     ,'w');
    fid(3)  = fopen([dirname,'/pr_star.dat'],'w');
    fid(4)  = fopen([dirname,'/pph.dat']    ,'w');
    fid(5)  = fopen([dirname,'/pth.dat']    ,'w');
    fid(6)  = fopen([dirname,'/H.dat']      ,'w');

    if obj.verbose
        disp('Writing files...')
    end
    
    for n=1:length(T)
        fprintf(fid(1), '%15.12f %15.12f %15.12f %15.12f  %15.12f %15.12f %15.12f %15.12f\n',[T(n) r(n) th(n) phi(n) dth(n) d2th(n) dphi(n) d2phi(n)]);
        fprintf(fid(2), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pr(n)   dpr(n)  d2pr(n)] );
        fprintf(fid(3), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) prs(n) dprs(n) d2prs(n)] );
        fprintf(fid(4), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pphi(n) dpphi(n) d2pphi(n)] );
        fprintf(fid(5), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pth(n) dpth(n) d2pth(n)] );
        fprintf(fid(6), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) H(n)     dH(n)   d2H(n)] ); 
    end

    nfid = length(fid);
    for k=1:nfid
        fclose(fid(k));
    end

    i = nfid+1;
    file = sprintf('%s/DB_OneCol_%s.dat',dirname,ID);
    fid(i)  = fopen(file,'w');

    fprintf(fid(i),'#nt %i\n',length(T));
    fprintf(fid(i),'#nv 21\n');
    fprintf(fid(i),'#abh %.6f\n',a);    
    for n=1:nvars*length(T)
        fprintf(fid(i),'%.15e\n',OneColumn(n));
    end
    fclose(fid(i));

elseif strcmp(obj.hamiltonian,'kerr_eq')
    a = DB_testmass_checks(obj,'Trying to write dynamics, but ');

    T = obj.t;

    x = obj.x;
    y = obj.y;
    z = obj.z;

    px = obj.px;
    py = obj.py;
    pz = obj.pz;

    [r,phi,th,pr,pphi,pth] = DB_coords_cart2spherical(x,y,z,px,py,pz);

    rc2  = r.^2 + a.^2.*(1+2./r);
    rc = sqrt(rc2);
    uc = 1./rc;

    H     = obj.Heff; % in the test-mass limit the Hamiltonian is the effective one 
    Horb  = obj.Horb;
    Gs    = obj.Gs./(r.*rc2);
    A     = obj.A;
    Bp    = obj.Bp;
    Bnp   = obj.Bnp;
    B     = 1./(Bp+Bnp);

    prs = sqrt(A./B).*pr;

    Omg = (Gs.*a+A.*uc.^2.*pphi./(Horb));

    dH     = DB_D1(H,T,4);
    d2H    = DB_D1(dH,T,4);

    dpr    = DB_D1(pr,T,4);
    d2pr   = DB_D1(dpr,T,4);

    dpth   = DB_D1(pth,T,4);
    d2pth  = DB_D1(dpth,T,4);   

    dpphi  = DB_D1(pphi,T,4);
    d2pphi = DB_D1(dpphi,T,4);

    dprs   = DB_D1(prs,T,4);
    d2prs  = DB_D1(dprs,T,4);

    dOmg   = DB_D1(Omg,T,4);

    OneColumn = [ T; r; th; phi; Omg; dOmg; H; dH; d2H; pr; dpr; d2pr; pth; dpth; d2pth; pphi; dpphi; d2pphi; prs; dprs; d2prs;];
    nvars = 21;

    if(length(OneColumn) ~= nvars*length(T) )
        error('!!!!! Missing vars !!!!! exit..');
    end 

    % Definitions for filename, useful for file names  
    rr_flag = 'cart';
    ham_flag = obj.hamiltonian;

    a_s  = sprintf('%.4f',a);
    r_s  = sprintf('%.3f',obj.r0);
    mu_s = sprintf('%.3e',nu);
    dirname = ['DB_dyn_',ham_flag,'_a',a_s,'_r0',r_s,'_nu',mu_s,'_',rr_flag];

    if exist(dirname,'dir')
        if obj.verbose
            disp(['Existing dir: ',dirname]);
        end
    else
        system(['mkdir ',dirname]);
        if obj.verbose
            disp(['Creating dir: ',dirname]);
        end
    end

    % file pointers
    fid(1)  = fopen([dirname,'/traj.dat']   ,'w');
    fid(2)  = fopen([dirname,'/pr.dat']     ,'w');
    fid(3)  = fopen([dirname,'/pr_star.dat'],'w');
    fid(4)  = fopen([dirname,'/pph.dat']    ,'w');
    fid(5)  = fopen([dirname,'/H.dat']      ,'w');

    if obj.verbose
        disp('Writing files...')
    end

    for n=1:length(T)
        fprintf(fid(1), '%15.12f %15.12f %15.12f  %15.12f %15.12f\n',[T(n) r(n)    phi(n)  Omg(n) dOmg(n)]);
        fprintf(fid(2), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pr(n)   dpr(n)  d2pr(n)] );
        fprintf(fid(3), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) prs(n) dprs(n) d2prs(n)] );
        fprintf(fid(4), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pphi(n) dpphi(n) d2pphi(n)] );
        fprintf(fid(5), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) H(n)     dH(n)   d2H(n)] ); 
    end

    nfid = length(fid);
    for k=1:nfid
        fclose(fid(k));
    end

    i = nfid+1;
    fid(i)  = fopen([dirname,'/DB_',rr_flag,'_a',a_s,'_r0',r_s,'_mu',num2str(nu, '%.3e'),'.dat'],'w');
    fprintf(fid(i),'#nt %i\n',length(T));
    fprintf(fid(i),'#nv 21\n');
    fprintf(fid(i),'#abh %.6f\n',a);
    for n=1:nvars*length(T)
        fprintf(fid(i),'%.15e\n',OneColumn(n));
    end
    fclose(fid(i));

elseif strcmp(obj.hamiltonian,'kerr_star')
    a = DB_testmass_checks(obj,'Trying to write dynamics, but ');
    disp('retreiving dynamics')
    T = obj.t;

    x = obj.x;
    y = obj.y;
    z = obj.z;

    px = obj.px;
    py = obj.py;
    pz = obj.pz;
    pxs = obj.pxs;
    pys = obj.pys;
    pzs = obj.pzs;

    [r,phi_w,th_w,pr,pphi,pth] = DB_coords_cart2spherical(x,y,z,px,py,pz);
    [~,~,~,prs] = DB_coords_cart2spherical(x,y,z,pxs,pys,pzs);
    
    phi = unwrap(phi_w);
    th = unwrap(th_w);
    H     = obj.Heff; % in the test-mass limit the Hamiltonian is the effective one 
    %{
    Horb = obj.Horb;
    Gs = obj.Gs;
    A = obj.A;
    rc2  = r.^2 + a.^2.*(1+2./r);
    rc = sqrt(rc2);
    uc = 1./rc;
    %}
    dHeff = obj.dHeff';
    Omg = -dHeff(:,1).*sin(phi)./(r.*sin(th))+dHeff(:,2).*cos(phi)./(r.*sin(th));

    dH     = DB_D1(H,T,4);
    d2H    = DB_D1(dH,T,4);

    dpr    = DB_D1(pr,T,4);
    d2pr   = DB_D1(dpr,T,4);

    dpth   = DB_D1(pth,T,4);
    d2pth  = DB_D1(dpth,T,4);

    dpphi  = DB_D1(pphi,T,4);
    d2pphi = DB_D1(dpphi,T,4);

    dprs   = DB_D1(prs,T,4);
    d2prs  = DB_D1(dprs,T,4);

    dOmg   = DB_D1(Omg,T,4);

    OneColumn = [ T; r; th; phi; Omg; dOmg; H; dH; d2H; pr; dpr; d2pr; pth; dpth; d2pth; pphi; dpphi; d2pphi; prs; dprs; d2prs;];
    nvars = 21;

    if(length(OneColumn) ~= nvars*length(T) )
        error('!!!!! Missing vars !!!!! exit..');
    end 

    % Definitions for filename, useful for file names  
    rr_flag = 'cart';
    ham_flag = obj.hamiltonian;

    a_s  = sprintf('%.4f',a);
    r_s  = sprintf('%.3f',obj.r0);
    mu_s = sprintf('%.3e',nu);
    dirname = ['DB_dyn_',ham_flag,'_a',a_s,'_r0',r_s,'_nu',mu_s,'_',rr_flag];

    if exist(dirname,'dir')
        if obj.verbose
            disp(['Existing dir: ',dirname]);
        end
    else
        system(['mkdir ',dirname]);
        if obj.verbose
            disp(['Creating dir: ',dirname]);
        end
    end

    % file pointers
    fid(1)  = fopen([dirname,'/traj.dat']   ,'w');
    fid(2)  = fopen([dirname,'/pr.dat']     ,'w');
    fid(3)  = fopen([dirname,'/pr_star.dat'],'w');
    fid(4)  = fopen([dirname,'/pph.dat']    ,'w');
    fid(5)  = fopen([dirname,'/H.dat']      ,'w');

    if obj.verbose
        disp('Writing files...')
    end

    for n=1:length(T)
        fprintf(fid(1), '%15.12f %15.12f %15.12f  %15.12f %15.12f\n',[T(n) r(n)    phi(n)  Omg(n) dOmg(n)]);
        fprintf(fid(2), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pr(n)   dpr(n)  d2pr(n)] );
        fprintf(fid(3), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) prs(n) dprs(n) d2prs(n)] );
        fprintf(fid(4), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) pphi(n) dpphi(n) d2pphi(n)] );
        fprintf(fid(5), '%15.12f %15.12f %15.12f  %15.12f \n',[T(n) H(n)     dH(n)   d2H(n)] ); 
    end

    nfid = length(fid);
    for k=1:nfid
        fclose(fid(k));
    end

    i = nfid+1;
    fid(i)  = fopen([dirname,'/DB_',rr_flag,'_a',a_s,'_r0',r_s,'_mu',num2str(nu, '%.3e'),'.dat'],'w');
    fprintf(fid(i),'#nt %i\n',length(T));
    fprintf(fid(i),'#nv 21\n');
    fprintf(fid(i),'#abh %.6f\n',a);
    for n=1:nvars*length(T)
        fprintf(fid(i),'%.15e\n',OneColumn(n));
    end
    fclose(fid(i));

end

if contains(obj.hamiltonian, 'kerr')
    DB_write_teuk_parfile(obj, dirname,ID)
end

system(['mv input_fields.dat ',dirname]);
dyn = obj;
save([dirname,'/dynamics.mat'],'dyn')

end