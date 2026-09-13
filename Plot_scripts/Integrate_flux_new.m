function [DeltaFlux,DeltaFlux_lm] = Integrate_flux_new(s,id,int_start,varargin)
    % ==================================================================================================================================================================
    % s: struct containing wf modes (including Dhlm)
    % id: can be either 'E' (energy flux) or 'J' (angular momentum flux)
    % int_start: start of integration. Can be either 'LR' or 'LSSO'
    % varargin: ellmax and multipoles. multipoles = 1 computes integrated flux of all multipoles
    % For 'J' only the norm is integrated (not the single components)
    % ==================================================================================================================================================================
    ellmax = 4;
    multipoles = 0;
    parity = 0; % applies parity transformation to modes before computing fluxes. Useful for paper p    lots

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'ellmax'
                    i      = i + 1;
                    ellmax = varargin{i};
                case 'multipoles'
                    i      = i + 1;
                    multipoles = varargin{i};
                case 'parity'
                    i      = i + 1;
                    parity = varargin{i};
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    if strcmp(id,'J')
        %multipoles = 0;
        [flux_v,t,flux_xlm,flux_ylm,flux_zlm] = DB_dJdt(s,'ellmax',ellmax,'parity',parity); %compute ang mom flux
        flux = sqrt(flux_v(:,1).^2+flux_v(:,2).^2+flux_v(:,3).^2); %ang mom flux norm
        for l=2:ellmax
            for m=-l:l
                if m<0 
                    flux_lm.ell(l).emminus(-m+1).flux = sqrt(flux_xlm.ell(l).emminus(-m+1).flux.^2 + flux_ylm.ell(l).emminus(-m+1).flux.^2 + flux_zlm.ell(l).emminus(-m+1).flux.^2);
                else
                    flux_lm.ell(l).emm(m+1).flux = sqrt(flux_xlm.ell(l).emm(m+1).flux.^2 + flux_ylm.ell(l).emm(m+1).flux.^2 + flux_zlm.ell(l).emm(m+1).flux.^2);
                end
            end
            flux_lm.ell(l).emminus(1).flux = flux_lm.ell(l).emm(1).flux;
        end
    elseif strcmp(id,'E')
        [flux,t,flux_lm] = DB_dEdt(s,'ellmax',ellmax);
    else
        error('id %s not recognised: use either J or E',id)
    end

    if strcmp(int_start,'LSSO')
        t_start = tLSSO_splined(s.dyn);
    elseif strcmp(int_start,'LR')
        t_start = tLR_splined(s.dyn);
    elseif strcmp(int_start,'BN')
        time = s.ell(2).emm(3).t;
        o22 = freq(s.ell(2).emm(3).hlm,time);
        idx0 = find(time>200,1);
        time_cut = time(idx0:end);
        idx1 = find(o22(idx0:end)>0.167,1);
        t_start = time_cut(idx1);
    else
        error('t_start %s not recognized: integration can start either at LR or LSSO',int_start)
    end

    idx0 = find(t>=t_start-10, 1, 'first');
    t_dense = t(idx0):0.01:t(end);

    flux_splined = spline(t,flux,t_dense);
    idx1 = find(t_dense>=t_start,1,'first');

    DeltaFlux_vec = cumtrapz(t_dense(idx1:end),flux_splined(idx1:end));

    DeltaFlux = DeltaFlux_vec(end);

    if multipoles~=0
        %if strcmp(id,'E')
            for l=2:ellmax
                for m=-l:l
                    if m<0
                        flux_lm_0 = flux_lm.ell(l).emminus(-m+1).flux;
                    else
                        flux_lm_0 = flux_lm.ell(l).emm(m+1).flux;
                    end
                    flux_lm_splined = spline(t,flux_lm_0,t_dense);
                    DeltaFlux_lm_vec = cumtrapz(t_dense(idx1:end),flux_lm_splined(idx1:end));
                    if m<0
                        DeltaFlux_lm.ell(l).emminus(-m+1).DeltaFlux = DeltaFlux_lm_vec(end);
                    else
                        DeltaFlux_lm.ell(l).emm(m+1).DeltaFlux = DeltaFlux_lm_vec(end);
                    end
                end
            end
            %{
        elseif strcmp(id,'J')
            for l=2:ellmax
                for m=-l:l
                    if m<0
                        flux_lm_x0 = flux_lm_x.ell(l).emminus(-m+1).flux;
                        flux_lm_y0 = flux_lm_y.ell(l).emminus(-m+1).flux;
                        flux_lm_z0 = flux_lm_z.ell(l).emminus(-m+1).flux;
                    else
                        flux_lm_x0 = flux_lm_x.ell(l).emm(-m+1).flux;
                        flux_lm_y0 = flux_lm_y.ell(l).emm(-m+1).flux;
                        flux_lm_z0 = flux_lm_z.ell(l).emm(-m+1).flux;
                    end
                    flux_lm_x_splined = spline(t,flux_lm_x0,t_dense);
                    flux_lm_y_splined = spline(t,flux_lm_y0,t_dense);
                    flux_lm_z_splined = spline(t,flux_lm_z0,t_dense);

                    DeltaFlux_lm_xvec = cumtrapz(t_dense(idx1:end),flux_lm_x_splined(idx1:end));
                    DeltaFlux_lm_yvec = cumtrapz(t_dense(idx1:end),flux_lm_y_splined(idx1:end));
                    DeltaFlux_lm_zvec = cumtrapz(t_dense(idx1:end),flux_lm_z_splined(idx1:end));
                    if m<0
                        DeltaFlux_lm.ell(l).emminus(-m+1).DeltaFlux = [DeltaFlux_lm_xvec(end) DeltaFlux_lm_yvec(end) DeltaFlux_lm_zvec(end)];
                    else
                        DeltaFlux_lm.ell(l).emm(m+1).DeltaFlux = [DeltaFlux_lm_xvec(end) DeltaFlux_lm_yvec(end) DeltaFlux_lm_zvec(end)];
                    end
                end
            end
        end
            %}
    else
        DeltaFlux_lm = 0;
    end

return
