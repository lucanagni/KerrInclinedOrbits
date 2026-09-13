classdef testKnownLimits < matlab.unittest.TestCase
%TESTKNOWNLIMITS DB_LR / DB_LSSO against closed-form literature results.
%
%   Schwarzschild (a=0): photon sphere at r=3M, ISCO at r=6M.
%   Kerr, equatorial (iota=0 prograde / iota=180 retrograde): closed-form
%   ISCO and photon-orbit radii from Bardeen, Press & Teukolsky (1972).
%
%   The photon-orbit reference formula used here for the prograde case
%   is the same one you already have sitting as a sanity-check comment
%   at DB_LR.m:61 ("equat = 2.*(1+cos((2./3).*acos(-a)));") -- this test
%   just turns that comment into an assertion that runs automatically,
%   and adds the retrograde branch and the ISCO formula alongside it.
%
%   Run with: runtests('testKnownLimits')

    properties (TestParameter)
        a = {0.1, 0.3, 0.5, 0.7, 0.9, 0.998};
    end

    methods (Test)

        function schwarzschildLimits(testCase)
            testCase.verifyEqual(DB_LR(0, 0), 3, 'AbsTol', 1e-10, ...
                'Schwarzschild photon sphere should be at r=3M.');
            testCase.verifyEqual(DB_LSSO(0, 0), 6, 'AbsTol', 1e-10, ...
                'Schwarzschild ISCO should be at r=6M.');
        end

        function equatorialISCO(testCase, a)
            r_pro   = DB_LSSO(a, 0);
            r_retro = DB_LSSO(a, 180);

            [r_pro_ref, r_retro_ref] = bardeenISCO(a);

            testCase.verifyEqual(r_pro, r_pro_ref, 'RelTol', 1e-6, ...
                sprintf('Prograde equatorial ISCO mismatch at a=%.3f', a));
            testCase.verifyEqual(r_retro, r_retro_ref, 'RelTol', 1e-6, ...
                sprintf('Retrograde equatorial ISCO mismatch at a=%.3f', a));
        end

        function equatorialPhotonOrbit(testCase, a)
            r_pro   = DB_LR(a, 0);
            r_retro = DB_LR(a, 180);

            r_pro_ref   = 2*(1 + cos((2/3)*acos(-a)));
            r_retro_ref = 2*(1 + cos((2/3)*acos( a)));

            testCase.verifyEqual(r_pro, r_pro_ref, 'RelTol', 1e-8, ...
                sprintf('Prograde equatorial photon orbit mismatch at a=%.3f', a));
            testCase.verifyEqual(r_retro, r_retro_ref, 'RelTol', 1e-8, ...
                sprintf('Retrograde equatorial photon orbit mismatch at a=%.3f', a));
        end
    end
end

function [r_pro, r_retro] = bardeenISCO(a)
    % Bardeen, Press & Teukolsky (1972), equatorial ISCO, M=1.
    Z1 = 1 + (1 - a^2)^(1/3) * ((1+a)^(1/3) + (1-a)^(1/3));
    Z2 = sqrt(3*a^2 + Z1^2);
    r_pro   = 3 + Z2 - sqrt((3 - Z1) * (3 + Z1 + 2*Z2));  % prograde  (minus sign)
    r_retro = 3 + Z2 + sqrt((3 - Z1) * (3 + Z1 + 2*Z2));  % retrograde (plus sign)
end
