function basedir = DB_GetBaseDir(flag, varargin)

    % Return basedir as a char-arrays. If return_as_string is 1, return it as a
    % Matlab-string. The code flag identify which basedir return. 
    % If in the future there will be more user, they will have to add a 
    % 'user_root' and to fill 'AllUserPaths'. If you insert a path that does 
    % not need the root befor, add them to 'add_no_root' in AllUserPaths().
    %
    % if verbose is 1, then print the final dir 
    
    if isempty(varargin)
        return_as_string = 0;
        verbose          = 0;
    elseif length(varargin)==1
        return_as_string = varargin{1};
        verbose          = 0;
    elseif length(varargin)==2
        return_as_string = varargin{1};
        verbose          = varargin{2};
    else
        error("Wrong input")
    end
    
    % check which user-root exist
    user_roots{1} = '/home/simone/';
    user_roots{2} = '/Users/nagar/';
    user_roots{3} = '/Users/andrea/';
    user_roots{4} = '/home/simonealbanesi/'; % stanzino 
    user_roots{5} = '/Users/simonealbanesi/';
    user_roots{6} = 'Fake/Users/Mattia'; % stanzino da aggiornare
    user_roots{7} = '\Users\Mattia'; % mypc
    user_roots{8} = '/home/luca/';
    
    
    code_flags{1}  = 'SimulationsTeukode';
    code_flags{2}  = 'separatrix';
    code_flags{3}  = 'EccRingdown';
    code_flags{4}  = 'SimulationsRWZ';
    code_flags{5}  = 'TestKOSMaps';
    code_flags{6}  = 'TestKOSPlots';
    code_flags{7}  = 'RingdownFit';
    code_flags{8}  = 'MartaSimulations';
    code_flags{9}  = 'RingdownRWZ';
    code_flags{10} = 'RelDiffPlots';
    code_flags{11} = 'datafiles';
    code_flags{12} = 'RainbowRingdown';
    code_flags{13} = 'TestDeltalm';
    code_flags{14} = 'rhoFD';
    code_flags{15} = 'MatlabTools_rwz';
    code_flags{16} = 'RingdownAnalysis';
    code_flags{17} = 'HyperbolicTests';
    code_flags{18} = 'SpinRemnant';
    code_flags{19} = 'NonCircular2PNTest';
    code_flags{20} = 'TestmassHypData';
    code_flags{21} = 'CompareInsplungeCache';
    code_flags{22} = 'EccRingdownMats';
    code_flags{23} = 'EccRingdownMats_neweccFphi';
    code_flags{24} = 'TestingRRcache';
    code_flags{25} = 'fluxesFD';
    code_flags{26} = 'hyperbolic';
    code_flags{27} = 'spin2PN_cache';
    code_flags{28} = 'spinCircInsplungeMats';
    code_flags{29} = 'dataQNMs';
    code_flags{30} = 'eccmerger_cache';
    code_flags{31} = 'RWZmats';
    
    nusers = numel(user_roots);
    nflags = numel(code_flags);
    allpaths = cell(nusers, nflags);
    for j=1:nusers
        for i=1:nflags
            allpaths{i,j} = '/';
        end
    end
    
    % user 1 - Simone
    allpaths{1,1}  = '/home/simone/data/simulations_teukode/';
    allpaths{2,1}  = '/home/simone/data/simulations_teukode/separatrix/';
    allpaths{3,1}  = '/home/simone/data/simulations_teukode/EccRingdown/';
    allpaths{4,1}  = '/home/simone/data/simulations_rwz/';
    allpaths{5,1}  = '/home/simone/data/test_miscellanea/TestKOSMaps/';
    allpaths{6,1}  = '/home/simone/data/test_miscellanea/TestKOSPlots/';
    allpaths{7,1}  = '/home/simone/data/ringdown_fit/';
    allpaths{8,1}  = '/media/simone/Extreme SSD/data/simulations_teuk/teuk_circ_proc6x2/';
    allpaths{9,1}  = '/home/simone/data/simulations_rwz/danilo/mats/ringdown/';
    allpaths{10,1} = '/home/simone/data/test_miscellanea/DiffPlots/';
    allpaths{11,1} = '/home/simone/repos/eob_eccentric_bis/TeukodeTests/MatlabScripts/datafiles/';
    allpaths{12,1} = '/home/simone/data/test_miscellanea/RainbowRingdown/';
    allpaths{13,1} = '/home/simone/data/test_miscellanea/TestDeltalm/';
    allpaths{14,1} = '/home/simone/repos/kerrorbitsolver/Mathematica/rholm_iResum_infty/Data_rho_Inftys/';
    allpaths{15,1} = '/home/simone/repos/kerrorbitsolver/MatlabTools_rwz/';
    allpaths{16,1} = '/home/simone/repos/kerrorbitsolver/RingdownAnalysis/';
    allpaths{17,1} = '/home/simone/data/test_miscellanea/Hyperbolic/';
    allpaths{18,1} = '/home/simone/data/spinRemnant/';
    allpaths{19,1} = '/home/simone/data/test_miscellanea/2PN_noncirc/';
    allpaths{20,1} = '/home/simone/data/simulations_teukode/hyperbolic/mats/';
    allpaths{21,1} = '/home/simone/data/caches/compare_insplunge_cache/';
    allpaths{22,1} = '/home/simone/data/simulations_teukode/EccRingdown/mats/';
    allpaths{23,1} = '/home/simone/data/simulations_teukode/EccRingdown/mats_neweccFphi/';
    allpaths{24,1} = '/home/simone/data/caches/testingRR_cache/';
    allpaths{25,1} = '/home/simone/data/dataFD/DataHughes/';
    allpaths{26,1} = '/home/simone/data/simulations_teukode/hyperbolic/';
    allpaths{27,1} = '/home/simone/data/caches/spin2PN/';
    allpaths{28,1} = '/home/simone/data/simulations_teukode/spinCircularInsplunge/mats/';
    allpaths{29,1} = '/home/simone/repos/kerrorbitsolver/KerrDynamics/';
    allpaths{30,1} = '/home/simone/data/caches/eccmerger/';
    allpaths{31,1} = '/home/simone/data/simulations_rwz/mats/';
    
    % user 2 - Alessandro
    allpaths{11,2} = '/Users/nagar/BITBUCKET_Repository/eob_eccentric_bis/TeukodeTests/MatlabScripts/datafiles/';
    allpaths{25,2} = '/Users/nagar/Data/DataHughes/';
    
    % user 3 - Andrea 
    allpaths{1,3}  = '/Users/andrea/GWdata/separatrix';
    allpaths{2,3}  = '/Users/andrea/GWdata/separatrix';
    allpaths{19,3} = '/Users/andrea/GWdata/test_miscellanea/2PN_noncirc/';
    allpaths{22,3} = '/Users/andrea/GWdata/simulations_teukode/EccRingdown/mats/';
    allpaths{23,3} = '/Users/andrea/GWdata/simulations_teukode/EccRingdown/mats_neweccFphi/';
    allpaths{28,3} = '/Users/andrea/GWdata/simulations_teukode/spinCircularInsplunge/mats/';
    allpaths{29,3} = '/Users/andrea/Git_Repositories/kerrorbitsolver/KerrDynamics/';
    allpaths{30,3} = '/Users/andrea/GWdata/caches/eccmerger';
    
    % user 4 - Simone again
    for i=1:nflags
        allpaths{i,4} = strrep(allpaths{i,1}, '/home/simone/', '/home/simonealbanesi/');
    end
    
    % user 5 - Simone again
    for i=1:nflags
        allpaths{i,5} = strrep(allpaths{i,1}, '/home/simone/', '/Users/simonealbanesi/');
    end
    
    % user 6 - Mattia stanzino
    allpaths{25,6} = '/home/Mattia/Panzieri/Data/DataHughes/';
    
    % user 7 - Mattia
    allpaths{25,7} = '\Users\Mattia\OneDrive\Documents\DataHughes\DataHughes\';
    
    % user 8 - Luca
    allpaths{25,8} = '/home/luca/Desktop/DataHughes/DataHughes/';
    
    user_idx = 0;
    flag_idx = 0;
    
    % find user
    for i=1:numel(user_roots)
        if isfolder(user_roots{i})
            user_idx = i;
            break;
        end
    end
    
    % find code flag
    for i=1:numel(code_flags)
        if strcmpi(flag, code_flags{i})
            flag_idx = i;
            break
        end
    end
    
    if ~user_idx
        error('User not found')
    end
    
    if ~flag_idx
        error('Flag not found')
    end
    
    basedir = allpaths{flag_idx, user_idx};
    
    if basedir(end)~='/' && basedir(end)~='\'
        basedir = [basedir, '/'];
    end
    
    if return_as_string
        basedir = string(basedir);
    end
    
    if verbose
        fprintf('flag:    %s\nbasedir: %s\n\n', flag, basedir)
    end
    
    %if ~isfolder(basedir)
    %    error("%s is not a valid path!\n", basedir)
    %end
    
    return
    
    
    