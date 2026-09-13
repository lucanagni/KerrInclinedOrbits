function DB_UpdateSavedStruct(varargin)
    if isempty(varargin)
        initial_angles = ["90" "60" "45" "30"];
        spins = ["02" "-02" "05" "-05" "09" "-09"];
        N = length(initial_angles);
        M = length(spins);
    else
        th0 = varargin{1};
        a = varargin{2};
        N = 1;
        M = 1;
    end

    for i=1:N
        if N~=1
            th0 = initial_angles(i);
        end
        for j = 1:M
            if M~=1
                a = spins(j);
            end

        inputDB = DB_LoadFromDat(th0,a);
        inputDB.verbose = 0;
        inputDB.teuk_output = 0;
        inputDB.plots = {};
        inputDB = rmfield(inputDB,'PA');

        DB = DB_class(inputDB);
        basedir = sprintf('/home/luca/waveforms/K/plunge/th0_%s/a%s',th0,a);

        s = DB_wf_struct(basedir,DB,'ellmax',4);
        structdir = sprintf('%s/wf.mat',basedir);
        save(structdir,'s');

        fprintf('Created %s\n',structdir)
        end
    end

return

function inputDB = DB_LoadFromDat(th0,a)

    fulldir = sprintf('/home/luca/waveforms/K/plunge/th0_%s/a%s/input_fields.dat',th0,a);

    fid = fopen(fulldir,'r');

    tline = fgetl(fid);
    while ischar(tline)
        eval(tline)
        tline = fgetl(fid);
    end

return
