function DB_write_teuk_parfile(obj, kerr_dir, ID)

if kerr_dir(end)=='/'
    kerr_dir = kerr_dir(1:end-1); %this is the directory of the Matlab dynamics
end

if obj.geodesics
    tlim = obj.t(end);
else
    % for non-geo motion, assume plunge and add what specified for ringdown (default is 500)
    tlim = obj.t(end)+obj.ringdown_extra;
end

% standard options
a = DB_testmass_checks(obj); % get kerr spin
grid_nx = obj.grid_nx;
grid_ny = obj.grid_ny;
sw = obj.spinweight;
r_plus    = get_outer_horizon_HH10(a);
cfl       = obj.evolve_cfl;
mpi_xsize = obj.mpi_xsize;
mpi_ysize = obj.mpi_ysize;

outdir_template = sprintf('%s/teuk_HH10_sw%+d_%s_%dx%d_m@m@',obj.teuk_outdir_folder,sw,ID,grid_nx,grid_ny);
outfile_template = sprintf('%s/teuk_HH10_sw%+d_%s_%dx%d_proc%dx%d_m@m@.par', kerr_dir,sw,ID,grid_nx,grid_ny,mpi_xsize,mpi_ysize);
out1d = obj.out1d;
lsum = obj.lsum;

one_col_file    = sprintf('DB_OneCol_%s.dat',ID);
traj_path = sprintf('%s/%s/%s', obj.trajectory_path,kerr_dir,one_col_file);
%traj_path = ['../trajectories/', kerr_dir, one_col_file];

% define common replacements
common_replacements = containers.Map;
common_replacements('kerr_abh')       = sprintf('%.6f', a);
common_replacements('kerr_spinfield') = sprintf('%d', sw);
common_replacements('grid_nx')        = sprintf('%d', grid_nx);
common_replacements('grid_ny')        = sprintf('%d', grid_ny);
common_replacements('grid_xmin')      = sprintf('%.32e', r_plus);
common_replacements('grid_xmin_.4f')  = sprintf('%.4f',  r_plus); 
common_replacements('mpi_xsize')      = sprintf('%d', mpi_xsize);
common_replacements('mpi_ysize')      = sprintf('%d', mpi_ysize);
common_replacements('grid_nx+2')      = sprintf('%d', grid_nx+2);
common_replacements('grid_ny+1')      = sprintf('%d', grid_ny+1);
common_replacements('evolve_cfl')     = sprintf('%.6f', cfl);
common_replacements('evolve_tmax')    = sprintf('%.2f', tlim);
common_replacements('particle_trajectory_file') = traj_path;
common_replacements('out1d')       = sprintf('%s', out1d);
common_replacements('lsum')       = sprintf('%s', lsum);


% read dummy and apply replacements common to all m
lines  = read_dummy();
nlines = numel(lines);
for i=1:nlines
    lines{i} = replace_placeholder(lines{i},common_replacements);
end

% Create parfiles for each Fourier m-mode
mmodes = -obj.emmax:1:obj.emmax;
mode_replacements = containers.Map;
fmt_m = '%+d'; % +d to force sign also for positive values

for m=mmodes
    % replace 'm' in the parfile
    mode_replacements('kerr_mmode') = sprintf('%d', m);
    mode_replacements('outdir')     = replace(outdir_template, '@m@', num2str(m, fmt_m) );
    lines_m = lines;
    for i=1:nlines
        lines_m{i} = replace_placeholder(lines_m{i},mode_replacements);
    end
    % write m-parfile
    outfile = replace(outfile_template, '@m@', num2str(m, fmt_m) );
    %outfile = sprintf('~/Desktop/a09dyn/%s',outfile); %output on Desktop, sometimes useful 
    fid = fopen(outfile, 'w');
    for j=1:nlines
        fprintf(fid, '%s\n', lines_m{j});
    end
    fprintf('Teukode parfile created: %s\n', outfile)
end
return

%##########################################################################
% External functions
%##########################################################################
function line = replace_placeholder(line,common_replacements)
% Function to apply replacements
if ~contains(line,'@')
    return
else
    keys = common_replacements.keys();
    for j=1:numel(keys)
        key = keys{j};
        toreplace = ['@',key,'@'];  
        if contains(line,toreplace)
            line = replace(line, toreplace, common_replacements(key));
            break;
        end
    end
end
return

function lines = read_dummy()
% Read dummy parfile for Teukode
filename = 'teuk.dummy';
fid = fopen(filename, 'r');
if fid == -1
    error('File could not be opened');
end
lines = {};
tline = fgets(fid);
while ischar(tline)
    line = strtrim(tline);  
    if isempty(line)
        line = ' ';
    end
    lines{end+1} = line; %#ok<AGROW> 
    tline = fgets(fid);
end
fclose(fid);
return

function r_plus = get_outer_horizon_HH10(a)
% Horizon in HH10 coordinates
M = 1;
S = 10;
r_plus = (S*a^2 + M*S^2 + sqrt(-(a^2)*(S^4) + (M^2)*(S^4)))/(a^2 + 2*M*S + S^2);
return
