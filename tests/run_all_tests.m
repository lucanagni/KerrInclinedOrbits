% RUN_ALL_TESTS  Run the whole regression suite and print a summary.
%
%   Usage: cd to (or add to path) the folder containing this file and
%   the core DB_*.m library, then run:
%       run_all_tests
%
%   Exits (errors) with a nonzero-failure message if anything failed, so
%   this can also be dropped into a CI step later if you ever want one.

here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
addpath(repoRoot);

results = runtests(here);

disp(table(results))

nFailed = sum([results.Failed]);
if nFailed > 0
    error('%d of %d test(s) failed.', nFailed, numel(results));
else
    fprintf('All %d tests passed.\n', numel(results));
end
