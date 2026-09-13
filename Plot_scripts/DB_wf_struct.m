function s = DB_wf_struct(dir,dyn,varargin)
    % ==========================================================================================================
    % Reads in data in a given directory (computed using dir and dyn) and generates a structure containing both
    % the dynamics and the modes specified in varargin
    % ==========================================================================================================

    ellmin = 2;
    ellmax = 4;
    mmax = -1;
    mmin = 0;
    obj.h = 1;
    obj.dh = 1;
    if dyn.geodesics == 1
        geod_flag = 'geod_';
    else
        geod_flag = '';
    end

    s.dyn = dyn;

    % Compute time shift to relate DB dynamics and teuk waveforms
    shift = Shift(dyn.t(1),dyn.r(1),dyn.chi1(3));

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'ellmax'
                    i     = i + 1;
                    ellmax = varargin{i};
                case 'ellmin'
                    i     = i + 1;
                    ellmin = varargin{i};
                case 'mmax'
                    i     = i + 1;
                    mmax = varargin{i};
                case 'mmin'
                    i     = i + 1;
                    mmin = varargin{i};
                case 'obj'
                    i     = i + 1;
                    obj = varargin{i};
                case 'newfiles'
                    newfiles_flag = 1;
                    geod_flag = 'sw-2_';
                case 'geod'
                    geod_flag = sprintf('%sgeod_',geod_flag);
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    % Get directory ID as computed in DB_write_teuk_parfile.m
    if ~newfiles_flag
        ID = sprintf('a%.4f_r0%.3f_th0%.3f_q%.0e',dyn.chi1(3),dyn.r0,rad2deg(dyn.th0),dyn.q);
    else
        ID = sprintf('a%.4f_th0%.3f_r0%.3f_q%.0e',dyn.chi1(3),rad2deg(dyn.th0),dyn.r0,dyn.q);
    end
    input_dir = sprintf('%s/teuk_HH10_%s%s_%dx%d_m',dir,geod_flag,ID,dyn.grid_nx,dyn.grid_ny);
    fprintf('%s\n',input_dir)
    %if exist(full_dir,'dir')==7
    %    input_dir = full_dir;
    %else
    %    nogrid_dir = sprintf('%s/teuk_HH10_%s%s_*_m',dir,geod_flag,ID,dyn.grid_nx,dyn.grid_ny);
    %    fprintf('directory not found: %s\n ')


    if obj.h==1
        % Useful if one only needs a single mode
        if ellmin==ellmax && mmin==mmax
            l = ellmin;
            m = mmin;
            if m>=0
                data = read_file(input_dir,'h',l,m);
                s.ell(l).emm(m+1).hlm = 2.*(data(:,2) + 1i.*data(:,3));
                s.ell(l).emm(m+1).t = data(:,1) + shift;
            else
                % If dynamics is equatorial compute -m from +m, otherwise look for teuk output
                n = abs(m);
                if abs(dyn.th0-pi/2)<1e-10
                    s.ell(l).emminus(n+1).hlm = (-1).^l.*conj(s.ell(l).emm(n+1).hlm);
                    s.ell(l).emminus(n+1).t = s.ell(l).emm(n+1).t;
                else
                    data = read_file(input_dir,'h',l,m);
                    s.ell(l).emminus(n+1).hlm = 2.*(data(:,2) + 1i.*data(:,3));
                    s.ell(l).emminus(n+1).t = data(:,1) + shift;
                end
            end

        else
            % Read all modes (at all m) from ellmin to ellmax
            for l=ellmin:ellmax
                for m=flip(-l:l)
                    if m>=0
                        data = read_file(input_dir,'h',l,m);
                        s.ell(l).emm(m+1).hlm = 2.*(data(:,2) + 1i.*data(:,3));
                        s.ell(l).emm(m+1).t = data(:,1) + shift;
                    else
                        n = abs(m);
                        if abs(dyn.th0-pi/2)<1e-10
                            %disp('computing negative m modes from positive ones')
                            s.ell(l).emminus(n+1).hlm = (-1).^l.*conj(s.ell(l).emm(n+1).hlm);
                            s.ell(l).emminus(n+1).t = s.ell(l).emm(n+1).t;
                        else
                            %disp('reading output for negative m modes')
                            data = read_file(input_dir,'h',l,m);
                            s.ell(l).emminus(n+1).hlm = 2.*(data(:,2) + 1i.*data(:,3));
                            s.ell(l).emminus(n+1).t = data(:,1) + shift;
                        end
                    end
                end
                s.ell(l).emminus(1).hlm = s.ell(l).emm(1).hlm;
                s.ell(l).emminus(1).t = s.ell(l).emm(1).t;
            end
        end
    end
    if obj.dh==1 % Same as above but reads dh instead of h

        if ellmin==ellmax && mmin==mmax
            l = ellmin;
            m = mmin;
            if m>=0
                data = read_file(input_dir,'dh',l,m);
                Hlm = 2.*(data(:,2) + 1i.*data(:,3));
                dT = data(:,1) + shift;
                time = s.ell(l).emm(m+1).t;

                if time_check(dT,time)
                    s.ell(l).emm(m+1).Dhlm = Hlm;
                else
                    fprintf('WARNING: splining dhlm for the l=%1f m=%1f mode\n',l,m)
                    s.ell(l).emm(m+1).Dhlm = spline(dT,Hlm,time);
                end

            else
                n = abs(m);
                if abs(dyn.th0-pi/2)<1e-10
                    s.ell(l).emminus(n+1).Dhlm = (-1).^l.*conj(s.ell(l).emm(n+1).Dhlm);
                else
                    data = read_file(input_dir,'dh',l,m);
                    Hlm = 2.*(data(:,2) + 1i.*data(:,3));
                    dT = data(:,1) + shift;
                    time = s.ell(l).emminus(n+1).t;

                    if time_check(dT,time)
                        s.ell(l).emminus(n+1).Dhlm = Hlm;
                    else
                        fprintf('WARNING: splining dhlm for the l=%1f m=%1f mode\n',l,m)
                        s.ell(l).emminus(n+1).Dhlm = spline(dT,Hlm,time);
                    end
                end
            end

        else
            for l=ellmin:ellmax
                for m=flip(-l:l)
                    if m>=0
                        data = read_file(input_dir,'dh',l,m);
                        Hlm = 2.*(data(:,2) + 1i.*data(:,3));
                        dT = data(:,1) + shift;
                        time = s.ell(l).emm(m+1).t;

                        if time_check(dT,time)
                            s.ell(l).emm(m+1).Dhlm = Hlm;
                        else
                            fprintf('WARNING: splining dhlm for the l=%1f m=%1f mode\n',l,m)
                            s.ell(l).emm(m+1).Dhlm = spline(dT,Hlm,time);
                        end
                    else
                        n = abs(m);
                        if abs(dyn.th0-pi/2)<1e-10
                            %disp('computing negative m modes from positive ones')
                            s.ell(l).emminus(n+1).Dhlm = (-1).^l.*conj(s.ell(l).emm(n+1).Dhlm);
                        else
                            %disp('reading output for negative m modes')
                            data = read_file(input_dir,'dh',l,m);
                            Hlm = 2.*(data(:,2) + 1i.*data(:,3));
                            dT = data(:,1) + shift;
                            time = s.ell(l).emminus(n+1).t;

                            if time_check(dT,time)
                                s.ell(l).emminus(n+1).Dhlm = Hlm;
                            else
                                fprintf('WARNING: splining dhlm for the l=%1f m=%1f mode\n',l,m)
                                s.ell(l).emminus(n+1).Dhlm = spline(dT,Hlm,time);
                            end
                        end
                    end
                end
                s.ell(l).emminus(1).Dhlm = s.ell(l).emm(1).Dhlm;
            end
        end

    end

return

function data = read_file(basedir,obj,l,m)
    n = abs(m);
    if m<0
        fulldir = sprintf('%s-%d/out0d/%s_Yl%dm-%d_x10.0000.dat',basedir,n,obj,l,n);
    else
        fulldir = sprintf('%s+%d/out0d/%s_Yl%dm%d_x10.0000.dat',basedir,n,obj,l,n);
    end

    if exist(fulldir,'file')==2
        data = readmatrix(fulldir);
    else
        %fprintf('Directory \n%s\n not found. Looking for directories with different grid size\n',fulldir)
        L = length(basedir);
        %gridsize = basedir(strfind(basedir,'m')-9:strfind(basedir,'m')-1);
        gridsize = basedir(L-9:L-1);
        basedir2 = strrep(basedir,gridsize,'*');
        if m<0
            gridless_fulldir = sprintf('%s-%d/out0d/%s_Yl%dm-%d_x10.0000.dat',basedir2,n,obj,l,n);
        else
            gridless_fulldir = sprintf('%s+%d/out0d/%s_Yl%dm%d_x10.0000.dat',basedir2,n,obj,l,n);
        end
        dirstruct = dir(gridless_fulldir);

        if exist(dirstruct.folder,'dir')==7
            fulldir2 = sprintf('%s/%s',dirstruct.folder,dirstruct.name);
            data = readmatrix(fulldir2);
            fprintf('Cannot find directory\n %s \nUsing\n %s \ninstead.\n',fulldir,fulldir2)
        else
            fprintf('%s\n',fulldir)
            error('Unable to locate the %d%d mode with any grid size',l,m)
        end
    end
return

function out = time_check(dt,t)
    if length(dt)~=length(t)
        out = 0;
    else
        if dt==t
            out = 1;
        else
            out = 0;
        end
    end
return
