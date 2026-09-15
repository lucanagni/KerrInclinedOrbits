classdef DB_class < handle

    properties
        r0     = 6;
        iota   = 0;        % Inclination angle, in radians
        th0    = pi/2;     % Historically gave IC with th0, but it's better to avoid it now. Keep a>0 and use iota.
        phi0   = 0;
        Tmin   = 0;
        Tmax   = 1e+6;
        dt     = 0.25;


        e0           = 0.;
        sr0          = 7.;
        pth0_ecc_ICs = 0.;
        anomaly0     = 0.;

        reltol = 1e-13;
        abstol = 1e-13;
        q      = 1e3;       % Use q>=1
        nu     = -1;        % If not given in input, updated from q.
                            % If given in input, update q accordingly.
        chi1   = [0;0;0.2]; % Kerr BH spin in the test-mass limit (if q>=1). MUST BE ALIGNED WITH z AXIS. Fixed parameter, not evolved.

        DBvKOS = 0; % DEBUG: used to compare some stuff with KOS

        verbose = 1;      % Switch-off to kill printings on console
        dashes  = '--------------------------------------------------';

        plots        = {}; % cell-array for plots (see DB_array_plotter)
        geodesics    = 0;  % geodesic motion
        ICs          = 'post-spherical';
        flux_KOS     = 0; %DEBUG: KOS flux for comparison
        flux_nucorrections = 0; %EOB flux with corrections in nu
        so_coupling  = 1; % bool. Turn off spin-orbit coupling in the Hamiltonian. Useful for debugging/checks

        %matlab stuff
        matlab_outdir_folder = '.'; %path for Matlab dynamics

        % Teukode stuff
        teuk_output  = 0;  % bool to write output for Teukode
        ringdown_extra = 500; % sets Teukode Tmax to allow ringdown
        mpi_xsize = 4; %number of processors to use in Teukode. Needed in parfile
        mpi_ysize = 2;
        grid_nx = 3601; %Number of points in the radial Teukode grid
        grid_ny = 321; %Number of points in the angular Teukode grid
        evolve_cfl = 2.0; %CFL factor for Teukode evolution
        spinweight = -2; %spin weight for Teukode perturbations
        out1d = 'none'; %(theta,t) output for Teukode. For fluxes use 'rphi iphi'
        outdir_folder = 'dyn'; %folder in /data/prometeo/luca.nagni/ where output is stored
        lsum = 'yes'; %To run horizon fluxes one needs to set this to "no". Also changes other settings in parfile
        emmax = 4; % Max m in parfiles
        teuk_outdir_folder = '/data/prometeo/luca.nagni'; %path where Teukode output will be saved
        trajectory_path = '../trajectories'; %path where the dynamics is located


        t
        x
        y
        z
        px
        py
        pz

        r
        th
        phi
        pr
        pth
        pphi
        L
        Lnorm
        F

        tLR
        tLSSO

        Heff
        Horb
        dHdp
        dHdx
        dHorbdp
        dHorbdx
        dHsodp
        dHsodx

        dHschwdx
        dHschwdp

        C
        Omg
        VOmg
        Omg_orb
        VOmg_orb
        Omg_so
        VOmg_so
        Omg_schw
        VOmg_schw
        VOmg_spin_only

        Omg_phi
        Omg_th
        Omg_orb_th
        Omg_orb_phi
        Omg_so_phi
        Omg_so_th

    end

    methods
        function obj = DB_class(varargin)
            % Constructor of the class
            input_fields = {};
            if ~isempty(varargin)
                input_struct = varargin{1};
                input_fields = fields(input_struct);
                for i=1:numel(input_fields)
                    f = input_fields{i};
                    obj.(f) = input_struct.(f);
                end
                if obj.teuk_output==1
                    write_input_file(obj.dt,input_struct);
                end
            end

            if obj.nu<=0
                obj.nu = obj.q/(1+obj.q)^2;
            else
                obj.q = (1+sqrt(1-4*obj.nu)-2*obj.nu)/(2*obj.nu);
                if obj.verbose
                    fprintf('+++ Updating q according to nu=%.3e: q=%.5f\n', obj.nu, obj.q)
                end
            end

            obj.check_input(input_fields); % also resolves {th0, iota} consistently

            [r,p] = DB_initial_conditions(obj);
            %I = pi/3; %Used to produce inclined eccentric dynamics. Should be removed at some point and implemented in eccentric initial conditions
            %R = [[cos(I), 0, -sin(I)];[0,1,0];[sin(I),0,cos(I)]];
            %r = R*r;
            if obj.verbose
                fprintf('%s\nInitial data\n%s\n', obj.dashes, obj.dashes)
                prec = 8;
                [~, X1, X2] = DB_nuX1X2(obj.q); % FIXME: we might get rid of DB_nuX1X2
                names  = {'q','nu','X1','X2','chi1','r','p'};
                values = {obj.q, obj.nu, X1, X2, obj.chi1, r, p};
                for i = 1:numel(names)
                    obj.print_with_prec(names{i}, values{i}, prec);
                end
            end

            %if obj.chi1(3)<0
            %    r(1) = -r(1);
            %    p(1) = -p(1);
            %    obj.chi1(3) = -obj.chi1(3);
            %end
            if obj.verbose
                fprintf('%s\nSolving ODEs\n%s\n', obj.dashes, obj.dashes)
            end
            ode_start = tic;
            obj.solve_ODEs(r,p);
            if obj.verbose
                toc(ode_start)
            end

            obj.vectorial_Hamiltonian();
            obj.Other_Quantities();

            if ~obj.geodesics && ~strcmp(obj.ICs,'elliptic')
                obj.find_tLR();
                obj.find_tLSSO();
            end

            if ~isempty(obj.plots)
                DB_array_plotter(obj, obj.plots)
            end

            if obj.teuk_output
                DB_write_kerr_dynamics(obj)
            end
        end

        function obj = check_input(obj, given_fields)
            if nargin<2
                given_fields = {};
            end
            gave_th0  = any(strcmp(given_fields, 'th0'));
            gave_iota = any(strcmp(given_fields, 'iota'));

            if gave_th0 && gave_iota
                error('Specify only one of th0 or iota, not both -- one determines the other.')
            end

            if obj.q<1
                error('Provided q<1! Use q>1 for test mass limit (q=m1/m2)')
            end

            nu_chk = obj.q/(1+obj.q)^2;
            if abs(nu_chk-obj.nu)>1e-15
                error('Inconsistency between nu and q!')
            end

            % check that spin is a (3,1) array
            if all(size(obj.chi1)==[1,3])
                obj.chi1 = transpose(obj.chi1);
            elseif ~all(size(obj.chi1)==[3,1])
                error('Invalid size for chi1!')
            end

            % check that ID are consistent
            DB_testmass_checks(obj, 'Using Kerr, but ');

            % ------------------------------------------------------------
            % Resolve {th0, iota}. iota is the user-facing convention
            % (a=chi1(3)>0 always; 0<iota<pi/2 prograde, pi/2<iota<pi
            % retrograde -- see e.g. Apostolatos et al / Hughes et al).
            % th0 is the internal/legacy polar angle used to build the
            % Cartesian initial position; from here on it is always a
            % DERIVED quantity, never independently authoritative.
            % ------------------------------------------------------------
            if obj.chi1(3) < 0
                % Legacy convention: retrograde encoded by a negative spin
                % instead of iota>pi/2. th0 (given, or its default) is
                % authoritative here; retrograde motion comes from
                % chi1(3)'s sign entering the Hamiltonian directly, NOT
                % from iota -- so iota is derived only for display/
                % bookkeeping and does not feed back into th0/phi0.
                if gave_iota
                    error(['iota was given together with a negative spin. The iota ', ...
                           'convention assumes chi1(3)>0 always; use a positive spin ', ...
                           'with iota>pi/2 for retrograde orbits instead.'])
                end
                if obj.th0<0 || obj.th0>pi/2
                    error('th0 must be in [0, pi/2].')
                end
                if obj.verbose
                    disp('Provided negative spin: the code will work, but a>0 with iota>pi/2 is the recommended way to set up retrograde orbits.')
                end
                obj.iota = pi/2 - obj.th0;
                obj.phi0 = 0;

            else
                % obj.chi1(3) >= 0: standard convention.
                if gave_th0
                    if obj.th0<0 || obj.th0>pi/2
                        error('th0 must be in [0, pi/2].')
                    end
                    obj.iota = pi/2 - obj.th0;  % th0 given with a>=0 is always prograde
                    if obj.verbose
                        disp('Provided a>0 with th0 as initial conditions. Iota will be computed consequently. It is recommended to specify a and iota instead of a and th0.')
                    end
                end
                % (if iota was given, or neither was given, obj.iota is
                % already correct: either user-supplied, or the default 0.)

                if obj.iota > pi/2
                    obj.th0  = obj.iota - pi/2;
                    obj.phi0 = pi;
                elseif 0 < obj.iota < pi/2
                    obj.th0  = pi/2 - obj.iota;
                    obj.phi0 = 0;
                else
                    obj.th0  = pi/2;
                    obj.phi0 = 0;
                end
            end

            % check that orbit starts outside of LSSO
            if obj.r0 < DB_LSSO(obj.chi1(3),obj.iota)
                fprintf('WARNING: initial distance r0 = %.3f is smaller than r_LSSO = %.3f. The orbit will likely be unstable\n',obj.r0,DB_LSSO(obj.chi1(3),obj.iota))
            end

             % check radiation reaction flags
            if obj.geodesics && strcmp(obj.ICs, 'post-spherical')
                disp('Ignoring PA flag for geodesic dynamics, using spherical ICs')
                obj.ICs = 'spherical';
            end

            if ~obj.geodesics && strcmp(obj.ICs, 'spherical')
                disp('WARNING: running non-conservative dynamics without PA correction!!!')
            end

            if obj.flux_nucorrections==1
                fprintf('WARNING: Using EOB flux with nu corrections!!!\n')
            end

            % check that trajectories don't end with / and if they do remove last character
            if obj.trajectory_path(end)=='/'
                obj.trajectory_path = obj.trajectory_path(1:end-1);
            end
            if obj.teuk_outdir_folder(end)=='/'
                obj.teuk_outdir_folder = obj.teuk_outdir_folder(1:end-1);
            end
            if obj.matlab_outdir_folder(end)=='/'
                obj.matlab_outdir_folder = obj.matlab_outdir_folder(1:end-1);
            end
        end

        function obj = solve_ODEs(obj, r, p)
            y0(1)  = r(1);
            y0(2)  = r(2);
            y0(3)  = r(3);
            y0(4)  = p(1);
            y0(5)  = p(2);
            y0(6)  = p(3);

            rend = 1+sqrt(1-obj.chi1(3)^2)+1e-4; % FIXME

            rhs     = @(t,y)  DB_rhs(obj, t, y);
            options = odeset('events',@(T,Y) DB_ode_stop(T,Y,rend), ...
                            'RelTol', obj.reltol, 'AbsTol', obj.abstol);
            [T,Y]   = ode113(rhs, obj.Tmin:obj.dt:obj.Tmax, y0, options);

            obj.x     = Y(:,1);
            obj.y     = Y(:,2);
            obj.z     = Y(:,3);
            obj.px    = Y(:,4);
            obj.py    = Y(:,5);
            obj.pz    = Y(:,6);
            obj.t = T;

            [obj.r,obj.phi,obj.th,obj.pr,obj.pphi,obj.pth] = DB_coords_cart2spherical(obj.x,obj.y,obj.z,obj.px,obj.py,obj.pz);
        end

        function obj = print_with_prec(obj, name, var, prec)
            str  = sprintf('%-8s : ', name);
            n    = length(var);
            if n>1
                str = [str, '['];
            else
                str = [str, ' '];
            end
            for i=1:length(var)
                str = [str, sprintf('%*.*f', prec+5, prec, var(i))]; %#ok<AGROW>
                if i<n
                    str = [str, ', ']; %#ok<AGROW>
                end
            end
            if n>1
                str = [str, ']'];
            end
            disp(str)
        end

        function obj = vectorial_Hamiltonian(obj)
            % F is only filled in below when radiation reaction is on
            % (obj.geodesics==0); for geodesic motion it stays zero.
            obj.F = zeros(length(obj.x),3);

            if obj.verbose
                fprintf('%s\nMemorizing Dynamics...\n%s\n', obj.dashes, obj.dashes)
            end
            write_start = tic;

            % Evaluate the whole trajectory at once (3xN) instead of
            % looping point-by-point: DB_metric_Kerr/DB_Hamiltonian_Kerr/
            % DB_flux2 are elementwise in their inputs, so a single
            % vectorized call replaces what used to be up to millions of
            % individual calls (each carrying its own function-call
            % overhead) -- this is what used to make memorization so slow.
            R = [obj.x.'; obj.y.'; obj.z.'];
            p = [obj.px.'; obj.py.'; obj.pz.'];

            L = cross(R,p,1);
            obj.L = L.';
            obj.Lnorm = vecnorm(L,2,1).';

            [Heff,Horb,dHeff,dHorb,dHso] = DB_Hamiltonian_Kerr(obj,R,p,obj.chi1);
            obj.Heff = Heff.';
            obj.Horb = Horb.';

            obj.dHdp = dHeff.dp.';
            obj.dHdx = dHeff.dx.';
            obj.dHorbdp = dHorb.dp.';
            obj.dHorbdx = dHorb.dx.';
            obj.dHsodp = dHso.dp.';
            obj.dHsodx = dHso.dx.';

            [~,~,~,dHschw] = DB_Hamiltonian_Kerr(obj,R,p,0*obj.chi1);
            obj.dHschwdx = dHschw.dx.';
            obj.dHschwdp = dHschw.dp.';

            if obj.geodesics==0
                obj.F = DB_flux2(R,p,dHeff.dp,obj.q,obj.chi1).';
            end

            if obj.verbose
                toc(write_start)
            end
        end

        function obj = Other_Quantities(obj)
            theta = obj.th;
            cth = cos(obj.th);
            sth = sin(obj.th);
            csth = csc(obj.th);

            cph = cos(obj.phi);
            sph = sin(obj.phi);

            a = norm(obj.chi1);
            u = 1./obj.r;

            %Carter's Constant
            obj.C = obj.pth.^2 + cth.^2.*(a.^2.*(1-obj.Heff.^2)+csth.^2.*obj.pphi.^2);

            %Angular Frequency
            dHeff = obj.dHdp;
            dHorb = obj.dHorbdp;
            dHso = obj.dHsodp;
            dHschw = obj.dHschwdp;
            %DEBUG -> different way to compute Omg (check if definitions are correct)
            %vr    = dHeff(:,1).*sin(theta).*cos(phi)+ dHeff(:,2).*sin(theta).*sin(phi) + dHeff(:,3).*cos(theta);
            %omega = u.*sqrt(dHeff(:,1).^2 + dHeff(:,2).^2 + dHeff(:,3).^2 - vr.^2);
            %obj.Omg = omega;

            dphdt = -dHeff(:,1).*sin(obj.phi).*csth.*u + dHeff(:,2).*cos(obj.phi).*csth.*u;
            dthdt = u.*cth.*(dHeff(:,1).*cos(obj.phi)+dHeff(:,2).*sin(obj.phi)) - dHeff(:,3).*u.*sin(theta);
            obj.VOmg = [-(sph.*dthdt+cth.*cph.*sth.*dphdt), cph.*dthdt - cth.*sth.*sph.*dphdt, sth.^2.*dphdt];
            obj.Omg = sqrt(dthdt.^2+sin(obj.th).^2.*dphdt.^2);

            dphdt_orb = -dHorb(:,1).*sin(obj.phi).*csth.*u + dHorb(:,2).*cos(obj.phi).*csth.*u;
            dthdt_orb = u.*cth.*(dHorb(:,1).*cos(obj.phi)+dHorb(:,2).*sin(obj.phi)) - dHorb(:,3).*u.*sin(theta);
            obj.VOmg_orb = [-(sph.*dthdt_orb+cth.*cph.*sth.*dphdt_orb), cph.*dthdt_orb - cth.*sth.*sph.*dphdt_orb, sth.^2.*dphdt_orb];
            obj.Omg_orb = sqrt(dthdt_orb.^2+sin(obj.th).^2.*dphdt_orb.^2);

            dphdt_so = -dHso(:,1).*sin(obj.phi).*csth.*u + dHso(:,2).*cos(obj.phi).*csth.*u;
            dthdt_so = u.*cth.*(dHso(:,1).*cos(obj.phi)+dHso(:,2).*sin(obj.phi)) - dHso(:,3).*u.*sin(theta);
            obj.VOmg_so = [-(sph.*dthdt_so+cth.*cph.*sth.*dphdt_so), cph.*dthdt_so - cth.*sth.*sph.*dphdt_so, sth.^2.*dphdt_so];
            obj.Omg_so = sqrt(dthdt_so.^2+sin(obj.th).^2.*dphdt_so.^2);

            dphdt_schw = -dHschw(:,1).*sin(obj.phi).*csth.*u + dHschw(:,2).*cos(obj.phi).*csth.*u;
            dthdt_schw = u.*cth.*(dHschw(:,1).*cos(obj.phi)+dHschw(:,2).*sin(obj.phi)) - dHschw(:,3).*u.*sin(theta);
            obj.VOmg_schw = [-(sph.*dthdt_schw+cth.*cph.*sth.*dphdt_schw), cph.*dthdt_schw - cth.*sth.*sph.*dphdt_schw, sth.^2.*dphdt_schw];
            obj.Omg_schw = sqrt(dthdt_schw.^2+sin(obj.th).^2.*dphdt_schw.^2);

            obj.VOmg_spin_only = obj.VOmg - obj.VOmg_schw;

            %vx = dHeff(:,1);
            %vy = dHeff(:,2);
            %vz = dHeff(:,3);
            %v = [vx,vy,vz];
            %R = [obj.x,obj.y,obj.z];

            %r_times_v = zeros(length(obj.t),1);
            %for i=1:length(obj.t)
            %    r_times_v(i) = norm(cross(R(i,:),v(i,:)));
            %    obj.Omg1(i,:) = u(i).^2.*norm(cross(R(i,:),v(i,:)));
            %end

            obj.Omg_phi = dphdt;
            obj.Omg_th = dthdt;
            obj.Omg_orb_phi = dphdt_orb;
            obj.Omg_orb_th = dthdt_orb;
            obj.Omg_so_phi = dphdt_so;
            obj.Omg_so_th = dthdt_so;

        end

        function obj = find_tLR(obj)
            a = obj.chi1(3);
            Iota = rad2deg(pi/2-obj.th0);
            r_LR = DB_LR(a,Iota);
            %end
            R = obj.r;
            T = obj.t;

            LR_pos = find(R<r_LR,1);

            idx0 = find(T>=T(LR_pos)-10, 1, 'first');

            R_short = R(idx0:end);
            t_short = T(idx0:end);
            newT = (T(LR_pos)-1):1e-6:(T(LR_pos)+1);

            sR = spline(t_short, R_short, newT);
            [~, idx_LR] = find(sR<r_LR,1);
            t_LR = newT(idx_LR);

            obj.tLR = t_LR;
        end

        function obj = find_tLSSO(obj)
            a = obj.chi1(3);
            iota = rad2deg(pi/2-obj.th0);
            r_LSSO = DB_LSSO(a,iota);
            %r_LSSO = 5.9865;
            R = obj.r;
            T = obj.t;

            LSSO_pos = find(R<r_LSSO,1);

            idx0 = find(T>=T(LSSO_pos)-10, 1, 'first');

            R_short = R(idx0:end);
            t_short = T(idx0:end);
            newT = (T(LSSO_pos)-2):1e-6:(T(LSSO_pos)+2);

            sR = spline(t_short, R_short, newT);
            idx_LSSO = find(sR<r_LSSO,1);
            t_LSSO = newT(idx_LSSO);

            obj.tLSSO = t_LSSO;
        end

    end
end

function write_input_file(dt,input_struct)
    input_fields = fields(input_struct);
    fid = fopen('input_fields.dat', 'w');
    dt_flag = 0;
    for i = 1:numel(input_fields)
        f = input_fields{i};  % Get the value from varargin
        if strcmp(f,'dt')
            dt_flag = 1;
        end
        value = input_struct.(f);
        fprintf(fid, 'inputDB.%s = ', input_fields{i});

        if isnumeric(value)
            if numel(value)==1
                if strcmp(input_fields{i},'th0')
                    den = rat(pi/value);
                    fprintf(fid, 'pi/(%s)\n', den);
                else
                    fprintf(fid, '%g;\n', value);
                end
                if strcmp(input_fields{i},'iota')
                    if value>pi/2
                        den = rat(pi/(value - pi/2));
                        fprintf(fid, 'pi/2 + pi/(%s)\n', den);
                    else
                        den = rat(pi/value);
                        fprintf(fid, 'pi/(%s)\n', den);
                    end
                else
                    fprintf(fid, '%g;\n', value);
                end
            else
                fprintf(fid, '[');  % Start cell array representation
                for j = 1:numel(value)
                    v = value(j);  % Extract cell content
                    if j~=numel(value)
                        fprintf(fid, '%g ',v);
                    else
                        fprintf(fid, '%g',v);
                    end
                end
                fprintf(fid,'];\n');
            end
        elseif ischar(value)
            fprintf(fid, '''%s'';\n', value);  % String
        elseif islogical(value)
            fprintf(fid, '%d;\n', value);  % Boolean (prints 0 or 1)
        elseif iscell(value)
            fprintf(fid, '{');  % Start cell array representation
            for j = 1:numel(value)
                v = value{j};  % Extract cell content
                if isnumeric(v)
                    fprintf(fid, '%g ', v);  % General numeric
                elseif ischar(v)
                    fprintf(fid, '''%s'' ', v);  % String
                elseif islogical(v)
                    fprintf(fid, '%d ', v);  % Boolean
                else
                    fprintf(fid, '[Unsupported Type] ');
                end
            end
            fprintf(fid, '};\n');  % End cell array representation
        else
            fprintf(fid, '[Unsupported Type]\n');
        end
    end
    if dt_flag==0
        fprintf(fid,'inputDB.dt = %f \n', dt);
    end
    fclose(fid);
end
