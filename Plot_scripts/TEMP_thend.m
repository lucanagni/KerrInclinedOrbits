function TEMP_thend
    masterdir = '~/waveforms/K/plunge';
    spins = ["02" "05" "09" "-02" "-05" "-09"];
    th0s = [90 60 45 30];
    for i=1:length(spins)
        for j=1:length(th0s)
            fulldir = sprintf('%s/th0_%d/a%s/wf.mat',masterdir,th0s(j),spins(i));
            if isfile(fulldir)
                load(fulldir)
                th_end = rad2deg(s.dyn.th(end));
                delta = s.dyn.th(end) - s.dyn.th(end-100);
                if delta>0
                    th_dot = '>0';
                else
                    th_dot = '<0';
                end
                iota = 90 - th0s(j);
                fprintf('a = %s, iota = %d: th_end = %.0f°, th_dot = %s \n',spins(i),iota,th_end,th_dot)
            end
        end
    end

    return