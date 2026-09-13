function [hpeak,tpeak,omg_peak] = Apeak(h,t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find waveform's amplitude peak position (tpeak) and value (hpeak)
    % omg_peak is the frequency value at peak: odd that it's computed here, it's only used in Tab3.m and whis was the easiest thing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idx_start = find(t>200,1); %only scan t>200 to remove initial junk 

    tmp_t = t(idx_start:end);
    tmp_amp_teuk = abs(h(idx_start:end));
    [~, Amax_idx] = max(tmp_amp_teuk);
    tpeakA_approx = tmp_t(Amax_idx);

    % splined wave
    idx0 = find(t>=tpeakA_approx-10, 1, 'first');
    A_short = abs(h(idx0:end));
    t_short = t(idx0:end);
    newT = (tpeakA_approx-2):1e-5:(tpeakA_approx+2);

    sA = spline(t_short, A_short, newT);
    [~, idx_peakA] = max(abs(sA));
    tpeak = newT(idx_peakA);
    hpeak = sA(idx_peakA);

    h_short = h(idx0:end);
    sh = spline(t_short, h_short, newT);
    omg = freq(sh,newT);
    omg_peak = omg(idx_peakA);
return

