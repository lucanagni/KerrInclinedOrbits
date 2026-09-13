function [DeltaE_LR,DeltaE_LSSO,DeltaJ_LR,DeltaJ_LSSO] = DB_flux_table(varargin)
    % ==========================================================================================================
    % Compute fluxes for table V of PaperI
    % ==========================================================================================================
    
    ellmax = 4;
    parity_flag = 0;

    % read varargin and eventually update default values
    i = 1;
    while i<=length(varargin)
        flag = varargin{i};
        lv   = lower(flag);
        if isstring(lv) || ischar(lv)
            switch lv
                case 'invert_parity'
                    parity_flag = 1;
                otherwise
                    error("'%s' is not a valid flag!", flag)
            end
        end
        i = i + 1;
    end

    masterdir = '~/waveforms/K/plunge';
    spins = ["02" "05" "09" "-02" "-05" "-09"];
    %spins = ["02" "-02"];
    th0s = [90 60 45 30];

    DeltaJ_LR = zeros(length(spins),length(th0s));
    DeltaJ_LSSO = DeltaJ_LR;
    DeltaE_LR = DeltaJ_LR;    
    DeltaE_LSSO = DeltaJ_LR;

    for i=1:length(spins)
        for j=1:length(th0s)
            fulldir = sprintf('%s/th0_%d/a%s/wf.mat',masterdir,th0s(j),spins(i));
            fprintf('a = %s, th0 = %d\n',spins(i),th0s(j))
            %if contains(spins(i),"09")
            %    ellmax=3;
            %else
            %    ellmax=4;
            %end

            if contains(spins(i),"02")
                splined = 1;
            else
                splined = 0;
            end
                                
            if isfile(fulldir)
                load(fulldir)
                DeltaE_LR(i,j) = Integrate_flux_new(s,'E','LR','ellmax',ellmax);
                DeltaE_LSSO(i,j) = Integrate_flux_new(s,'E','LSSO','ellmax',ellmax);
                DeltaJ_LR(i,j) = Integrate_flux_new(s,'J','LR','ellmax',ellmax,'parity',parity_flag);
                DeltaJ_LSSO(i,j) = Integrate_flux_new(s,'J','LSSO','ellmax',ellmax,'parity',parity_flag);
                %[DeltaJ_LR(i,j),DeltaJ_LSSO(i,j)] = Integrate_flux(s);
            else
                DeltaJ_LR(i,j) = 0;
                DeltaJ_LSSO(i,j) = 0;
            end
        end
    end

    n=0;
    for i=1:length(spins)/2
        for j=1:length(th0s)
            n = n + 1;
            fulldir = sprintf('%s/th0_%d/a%s/wf.mat',masterdir,th0s(j),spins(i));
            load(fulldir)
            Label = KerrLabel(s.dyn);
            fprintf('%d & %s & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ \n',n,Label,DeltaE_LR(i,j),DeltaE_LSSO(i,j),DeltaJ_LR(i,j),DeltaJ_LSSO(i,j))
        end
        for j=1:length(th0s)
            n = n + 1;
            fth0s = flip(th0s);
            fulldir = sprintf('%s/th0_%d/a-%s/wf.mat',masterdir,fth0s(j),spins(i));
            load(fulldir)
            Label = KerrLabel(s.dyn);
            fprintf('%d & %s & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ \n',n,Label,DeltaE_LR(i+(length(spins)/2),length(th0s)+1-j),DeltaE_LSSO(i+(length(spins)/2),length(th0s)+1-j),DeltaJ_LR(i+(length(spins)/2),length(th0s)+1-j),DeltaJ_LSSO(i+(length(spins)/2),length(th0s)+1-j))
        end
    end


return