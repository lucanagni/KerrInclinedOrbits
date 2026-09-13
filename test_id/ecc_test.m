function ecc_test
    inputDB.geodesics = 0;
    inputDB.e0        = 0.5;
    inputDB.sr0       = 6.6; %used to be 7.1
    inputDB.ICs       = 'eccentric';
    inputDB.Tmax      = 10e+3;
    inputDB.anomaly0  = 0.;
    inputDB.pth0_ecc_ICs = 0.1; %used to be 0.1
    inputDB.teuk_output = 1;
    %inputDB.chi1(3) = 0.5; %used to be 0.2 (default)
    
    db1 = DB_class(inputDB);
    
    inputDB.anomaly0 = 0.5;
    inputDB.teuk_output = 0;
    %db2 = DB_class(inputDB);
    db2 = db1;
    
    r1 = sqrt(db1.x.^2 + db1.y.^2 + db1.z.^2);
    r2 = sqrt(db2.x.^2 + db2.y.^2 + db2.z.^2);
    
    [~,ipks1] = findpeaks(r1);
    [~,ipks2] = findpeaks(r2);
    
    tau_vec = {db1.t-db1.t(ipks1(1)), db2.t-db2.t(ipks2(1))};
    db_vec = {db1, db2};
    
    figure 
    subplot(1,2,1)
    plot(db1.t, r1)
    subplot(1,2,2)
    plot3(db1.x, db1.y, db1.z)
    drawnow 

    figure
    for i=1:2
        db = db_vec{i};
        tau = tau_vec{i};
        r = sqrt(db.x.^2 + db.y.^2 + db.z.^2);
        subplot(2,2,1)
        plot(tau, r)
        subplot(2,2,2)
        plot(db.x, db.y)
        subplot(2,2,3)
        plot3(db.x, db.y, db.z)
        subplot(2,2,4)
        plot(tau, db.z)
        for j=1:4
            subplot(2,2,j)
            hold on
        end
        fprintf('spin = %.2f\n',db.chi1(3))
    end
    drawnow
return
