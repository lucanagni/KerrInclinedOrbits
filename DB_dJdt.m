function [dJdt,time,dJdtx_lm,dJdty_lm,dJdtz_lm] = DB_dJdt(s,varargin)
    % ==================================================================================================================================================================
    % Compute angular momentum flux along x,y,z
    % Takes in input a structure s where
    % s.ell(l).emm(m).hlm are the waveform modes
    % s.ell(l).emm(m).Dhlm are the time derivatives of the modes (as given by Teukode)
    % The factor 2 comes from the conventions defining Psi4 in Teukode
    % ==================================================================================================================================================================
    ellmax = 4;
    splined = 1;
    parity = 0; % applies parity transformation to modes before computing fluxes. Useful for paper plots

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'ellmax'
                    i       = i + 1;
                    ellmax  = varargin{i};
                case 'splined'
                    i       = i + 1;
                    splined = varargin{i};
                case 'parity'
                    i       = i + 1;
                    parity = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    %time = s.T;
    time = s.ell(2).emm(3).t;

    dJdt_x = zeros(length(time),1);
    dJdt_y = zeros(length(time),1);
    dJdt_z = zeros(length(time),1);

    for l = 2:ellmax
        for m = -l:l
            %fprintf('l = %d, m = %d\n',l,m)
            n = abs(m);

            if splined %when hlm arrays have different lengths, spline to lm=22

                if m<0
                    Dhlm = spline(s.ell(l).emminus(n+1).t,s.ell(l).emminus(n+1).Dhlm,time);
                    hlm = spline(s.ell(l).emminus(n+1).t,s.ell(l).emminus(n+1).hlm,time);
                else
                    Dhlm = spline(s.ell(l).emm(m+1).t,s.ell(l).emm(m+1).Dhlm,time);
                    hlm = spline(s.ell(l).emm(m+1).t,s.ell(l).emm(m+1).hlm,time);
                end

                if m<0
                    Dhlm_p1 = spline(s.ell(l).emminus(n).t,s.ell(l).emminus(n).Dhlm,time); %emminus 3 because usually the problematic modes are m=4;
                    if abs(m)==l
                        Dhlm_m1 = 0;
                    else
                        Dhlm_m1 = spline(s.ell(l).emminus(n+2).t,s.ell(l).emminus(n+2).Dhlm,time);
                    end
                elseif m>0
                    Dhlm_m1 = spline(s.ell(l).emm(m).t,s.ell(l).emm(m).Dhlm,time);
                    if m==l
                        Dhlm_p1 = 0;
                    else
                        Dhlm_p1 = spline(s.ell(l).emm(m+2).t,s.ell(l).emm(m+2).Dhlm,time);
                    end
                else
                    Dhlm_p1 = spline(s.ell(l).emm(2).t,s.ell(l).emm(2).Dhlm,time);
                    Dhlm_m1 = spline(s.ell(l).emminus(2).t,s.ell(l).emminus(2).Dhlm,time);
                end

            else

                if m<0
                    Dhlm = s.ell(l).emminus(n+1).Dhlm;
                    hlm = s.ell(l).emminus(n+1).hlm;
                else
                    Dhlm = s.ell(l).emm(m+1).Dhlm;
                    hlm = s.ell(l).emm(m+1).hlm;
                end

                if m<0
                    Dhlm_p1 = s.ell(l).emminus(n).Dhlm;
                    if abs(m)==l
                        Dhlm_m1 = 0;
                    else
                        Dhlm_m1 = s.ell(l).emminus(n+2).Dhlm;
                    end
                elseif m>0
                    Dhlm_m1 = s.ell(l).emm(m).Dhlm;
                    if m==l
                        Dhlm_p1 = 0;
                    else
                        Dhlm_p1 = s.ell(l).emm(m+2).Dhlm;
                    end
                else
                    Dhlm_p1 = s.ell(l).emm(2).Dhlm;
                    Dhlm_m1 = s.ell(l).emminus(2).Dhlm;
                end
            end

            if parity
                hlm = (-1).^m.*conj(hlm);
                Dhlm = (-1).^m.*conj(Dhlm);
                Dhlm_m1 = (-1).^m.*conj(Dhlm_m1);
                Dhlm_p1 = (-1).^m.*conj(Dhlm_p1);
            end

            Jx = imag(hlm.*(f(l,m).*conj(Dhlm_p1) + f(l,-m).*conj(Dhlm_m1)));
            Jy = - real(hlm.*(f(l,m).*conj(Dhlm_p1) - f(l,-m).*conj(Dhlm_m1)));
            Jz = m*imag(hlm.*conj(Dhlm)); 

            dJdt_x = dJdt_x + Jx;
            dJdt_y = dJdt_y + Jy; 
            dJdt_z = dJdt_z + Jz;

            if m<0
                dJdtx_lm.ell(l).emminus(n+1).flux = Jx./32./pi;
                dJdty_lm.ell(l).emminus(n+1).flux = Jy./32./pi;
                dJdtz_lm.ell(l).emminus(n+1).flux = Jz./16./pi;
            else
                dJdtx_lm.ell(l).emm(n+1).flux = Jx./32./pi;
                dJdty_lm.ell(l).emm(n+1).flux = Jy./32./pi;
                dJdtz_lm.ell(l).emm(n+1).flux = Jz./16./pi;
            end
        end 
        dJdtx_lm.ell(l).emminus(1).flux =  dJdtx_lm.ell(l).emm(1).flux;
        dJdty_lm.ell(l).emminus(1).flux =  dJdtx_lm.ell(l).emm(1).flux;
        dJdtz_lm.ell(l).emminus(1).flux =  dJdtx_lm.ell(l).emm(1).flux;
    end

    dJdt_x = dJdt_x/32/pi;
    dJdt_y = dJdt_y/32/pi;
    dJdt_z = dJdt_z/16/pi;

    dJdt = [dJdt_x dJdt_y dJdt_z];
return

function f = f(l,m)
    f = sqrt(l*(l+1) - m*(m+1));
return