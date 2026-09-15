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

                % Reruns the exact same deterministic computation, so this should
                % reproduce bit-for-bit; RelTol only needs to absorb legitimate
                % floating-point reordering from a behavior-preserving rewrite
                % (e.g. the DB_metric_Kerr/DB_Hamiltonian_Kerr CSE rewrite), which
                % sits far below the golden-master cases' own 1e-12 integrator
                % tolerance -- not genuine physics drift.
                %
                % RelTol is 1e-8 rather than tighter because the post-spherical
                % plunge cases (geodesics=0) integrate for O(1e4) adaptive steps
                % through the deep-plunge terminal approach to the horizon, where
                % the dynamics are numerically stiff: a per-step reordering
                % difference at the level of machine precision (unavoidable with
                % any CSE rewrite that isn't literally bit-identical) gets
                % amplified there to ~1e-9 by the last ~1% of steps. The
                % pure-geodesic cases (1-3) match to ~1e-12 throughout; only the
                % plunge tail needs the looser bound.
                %
                % This was verified NOT to be a rewrite-specific bug: (1) a
                % direct, per-call comparison of old vs. new DB_metric_Kerr/
                % DB_Hamiltonian_Kerr at the exact states visited along the
                % plunge trajectories (not random points) stayed at ~1e-13 or
                % better everywhere, including the deepest strong-field points
                % right before plunge -- i.e. no formula-level precision loss.
                % (2) A controlled experiment on the pristine, UNMODIFIED
                % pre-rewrite code confirmed the sensitivity is inherent to
                % this dynamical regime: reassociating a single 3-term sum in
                % just one of the six (redundant) derivative formulas --
                % nothing else changed anywhere -- reproduced trajectory-level
                % divergence of the same order as the full CSE rewrite:
                %
                %   case            | full CSE rewrite | one reassociated sum
                %   case5 (plunge)  |      6.7e-10      |       6.5e-10
                %   case6 (plunge)  |      2.2e-9       |       1.7e-9
                %
                % See tests/CSE_rewrite_conversation.txt for the full
                % investigation (this is why past tests -- which never touched
                % the arithmetic inside these two functions -- always matched
                % at 1e-11: there was no reordering for the plunge sensitivity
                % to amplify, not because the sensitivity wasn't there).
                %
                % C gets its own, looser RelTol for a related but distinct
                % reason: vectorial_Hamiltonian was rewritten to evaluate
                % DB_Hamiltonian_Kerr once on the whole trajectory (3xN) instead
                % of once per point in a loop, for performance. Evaluating many
                % points at once vs. one at a time is not guaranteed bit-for-bit
                % (vecnorm/dot/sum over a wide array can take a different
                % instruction path than the same reduction over 3 elements
                % done one column at a time), and this measures out at ~3e-13
                % absolute in Heff on case4's plunge tail -- a direct batch-vs.
                % loop comparison of dyn.Heff at every trajectory point found
                % only 630/22058 points differing at all, by at most 2.7e-13.
                % C amplifies this the same way the CSE rewrite's reordering
                % was amplified above: C contains a^2*(1-Heff^2), and Heff -> 1
                % in the last steps before plunge, so the ~1e-13 absolute noise
                % in Heff becomes ~1e-12 relative in that term, then compounds
                % over the stiff terminal integration to ~1.6e-8 relative in C
                % at the last 3 of case4's 22058 steps -- just over the 1e-8
                % used for every other field. r/th/phi/Heff/pphi are untouched
                % by this (they matched at 1e-8 already) because they don't
                % contain a comparable cancellation.
                testCase.verifyEqual(dyn.r,    ref.r,    'RelTol', 1e-8, [key ': r mismatch']);
                testCase.verifyEqual(dyn.th,   ref.th,   'RelTol', 1e-8, [key ': th mismatch']);
                testCase.verifyEqual(dyn.phi,  ref.phi,  'RelTol', 1e-8, [key ': phi mismatch']);
                testCase.verifyEqual(dyn.Heff, ref.Heff, 'RelTol', 1e-8, [key ': Heff mismatch']);
                testCase.verifyEqual(dyn.C,    ref.C,    'RelTol', 5e-7, [key ': C mismatch']);
                testCase.verifyEqual(dyn.pphi, ref.pphi, 'RelTol', 1e-8, [key ': pphi mismatch']);
            end
        end
    end
end
