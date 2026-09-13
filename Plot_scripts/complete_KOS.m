function new_s = complete_KOS(s)
    % ==========================================================================================================
    % Input: KOS structure. Completes it with missing stuff for some script (I don't remember which)
    % ========================================================================================================== 
    t = s.T;

    for l=2:4
        norm = sqrt((l+2)*(l+1)*l*(l-1));
        for m=0:l
            s.ell(l).emm(m+1).hlm = norm.*s.ell(l).emm(m+1).psi;
        end
    end

    for l=2:4
        for m=-l:-1
            hlm = s.ell(l).emm(-m+1).hlm;
            s.ell(l).emminus(-m+1).hlm = (-1).^l.*conj(hlm);
        end

        s.ell(l).emminus(1).hlm = s.ell(l).emm(1).hlm;
    end

    for l=2:4
        for m=-l:l
            if m<0 
                hlm = s.ell(l).emminus(-m+1).hlm;
                Dhlm = DB_D1(hlm,t,4);
                s.ell(l).emminus(-m+1).Dhlm = Dhlm;
            else
                hlm = s.ell(l).emm(m+1).hlm;
                Dhlm = DB_D1(hlm,t,4);
                s.ell(l).emm(m+1).Dhlm = Dhlm;
            end
        end
        
        s.ell(l).emminus(1).Dhlm = s.ell(l).emm(1).Dhlm;
    end

    new_s = s;

return
