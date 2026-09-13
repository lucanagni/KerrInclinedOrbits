function lbl = KerrLabel(dyn)
    % ==========================================================================================================
    % Label Dynamics for Plots
    % ==========================================================================================================
    a = dyn.chi1(3);
    A = abs(a)*10;

    if a~=0
        I = rad2deg(pi/2 - sign(a).*dyn.th0);
    else
        I = rad2deg(pi/2 - dyn.th0);
    end

    if dyn.e0 ~=0
        I = 20;
    end
    lbl = sprintf('{\\tt a0%1.0fi%1.f}',A,I);
return
