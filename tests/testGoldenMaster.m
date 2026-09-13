classdef testGoldenMaster < matlab.unittest.TestCase
%TESTGOLDENMASTER Regression test against a saved reference snapshot.
%
%   Reruns the fixed set of cases captured by generate_golden_master.m
%   and checks the output still matches golden_master.mat. This is the
%   safety net you want in place BEFORE doing the DB_metric_Kerr /
%   DB_Hamiltonian_Kerr CSE rewrite: run generate_golden_master.m first
%   on the current (pre-rewrite) code, then keep running this test
%   throughout the rewrite. As long as it stays green, the rewrite is
%   behavior-preserving to numerical precision.
%
%   Run with: runtests('testGoldenMaster')

    properties
        golden
    end

    methods (TestClassSetup)
        function loadGolden(testCase)
            here = fileparts(mfilename('fullpath'));
            goldenFile = fullfile(here, 'golden_master.mat');
            testCase.assumeTrue(isfile(goldenFile), ...
                'golden_master.mat not found -- run generate_golden_master.m first.');
            s = load(goldenFile, 'golden');
            testCase.golden = s.golden;
        end
    end

    methods (Test)
        function matchesGoldenCases(testCase)
            keys = fieldnames(testCase.golden);
            for i = 1:numel(keys)
                key = keys{i};
                ref = testCase.golden.(key);
                dyn = DB_class(ref.input);

                testCase.verifyEqual(dyn.r,    ref.r,    'RelTol', 1e-8, [key ': r mismatch']);
                testCase.verifyEqual(dyn.th,   ref.th,   'RelTol', 1e-8, [key ': th mismatch']);
                testCase.verifyEqual(dyn.phi,  ref.phi,  'RelTol', 1e-8, [key ': phi mismatch']);
                testCase.verifyEqual(dyn.Heff, ref.Heff, 'RelTol', 1e-8, [key ': Heff mismatch']);
                testCase.verifyEqual(dyn.C,    ref.C,    'RelTol', 1e-8, [key ': C mismatch']);
                testCase.verifyEqual(dyn.pphi, ref.pphi, 'RelTol', 1e-8, [key ': pphi mismatch']);
            end
        end
    end
end
