function [Omgpeak,tpeak] = OmgPeak(dyn,which)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find waveform's amplitude peak position (tpeak) and value (hpeak)
    % omg_peak is the frequency value at peak: odd that it's computed here, it's only used in Tab3.m and whis was the easiest thing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = dyn.t;
    if strcmp(which,'orb')
        Omg = dyn.Omg_orb;
    else
        Omg = dyn.Omg;
    end
    Omg_phi = dyn.Omg_phi;

    %idx1 = find(t>tLR_splined(dyn)-2800,1);
    %idx2 = find(t>tLR_splined(dyn)+10,1);
    idx1 = 1;
    idx2 = length(t);

    Omg = Omg(idx1:idx2);
    t = t(idx1:idx2);

    [~,idx_max1] = max(Omg);
    newT = (t(idx_max1)-2):1e-6:(t(idx_max1)+2);
    sOmg = spline(t, Omg, newT);
    [~, idx_peakOmg] = max(sOmg);
    tpeak = newT(idx_peakOmg);
    Omgpeak = sOmg(idx_peakOmg);

    %{
    idx_end = length(t);

    if Omg_phi(end)<0
        while Omg(idx_end)==max(Omg)
            idx_end= idx_end-1;
        end
        %idx = find(flip(Omg_phi)>0,1);
        %idx_end = length(t)-idx
    end
    %idx_end;

    tmp_t = t(1:idx_end);
    tmp_Omg = Omg(1:idx_end);
    [~, max_idx] = max(tmp_Omg);
    tpeak_approx = tmp_t(max_idx);

    % splined wave
    idx0 = find(t>=tpeak_approx-10, 1, 'first');
    Omg_short = tmp_Omg(idx0:end);
    t_short = tmp_t(idx0:end);
    newT = (tpeak_approx-2):1e-5:(tpeak_approx+2);

    sOmg = spline(t_short, Omg_short, newT);
    [~, idx_peakOmg] = max(sOmg);
    tpeak = newT(idx_peakOmg);
    Omgpeak = sOmg(idx_peakOmg);
    %}
return

