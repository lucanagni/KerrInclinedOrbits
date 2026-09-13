function TEMP_peakchek(a,iota)
    inputDB.verbose = 0;
    inputDB.chi1(3) = a.*sign(90-iota);
    inputDB.th0 = deg2rad(abs(90-iota));
    inputDB.r0 = DB_LSSO(a,iota) + 0.5;
    DB = DB_class(inputDB);

    DB_compare_freq(DB);
    figure
    findpeaks(DB.Omg,DB.t)
    xline(tLR_splined(DB),'Color',[.7 .7 .7],'LineStyle','--')
return