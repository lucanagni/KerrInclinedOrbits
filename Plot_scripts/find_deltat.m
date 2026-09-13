function delta = find_deltat(h1,t1,h2,t2)
    % find dt to match peaks of waveforms h1 and h2
    % add dt to t1 to match

    [~,tp1] = Apeak(h1,t1);
    [~,tp2] = Apeak(h2,t2);
    delta = tp2-tp1;

    if abs(delta)<1e-4
        delta = 0;
    end
return

