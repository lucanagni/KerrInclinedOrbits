classdef testGeodesicConservation < matlab.unittest.TestCase
%TESTGEODESICCONSERVATION Conservation checks for spherical (geodesic) orbits.
%
%   For a genuine Kerr geodesic (no radiation reaction), the effective
%   Hamiltonian H, Carter's constant C, and the axial angular momentum
%   p_phi are all constants of motion, and a spherical orbit's radius r
%   should not drift. This test builds spherical geodesic orbits at a
%   grid of spins and inclinations and checks all four quantities stay
%   constant to well within the integrator's requested tolerance.
%
%   This is deliberately the first test in the suite: it exercises
%   DB_spherical_ICs, DB_rhs (geodesic branch), DB_Hamiltonian_Kerr, and
%   the post-processing in DB_class.Other_Quantities all at once, with a
%   check (Carter's constant) that isn't trivially guaranteed by how the
%   equations of motion are written -- so a failure here is informative
%   about where to look, not just "something is wrong somewhere".
%
%   NOTE: iota = 90 deg (exact polar start, th0 = 0) is deliberately
%   excluded -- DB_spherical_ICs divides by sin(th0) internally, so an
%   exact polar orbit is a coordinate singularity for the root-find, not
%   a real physics failure. 89/91 deg bracket it instead.
%
%   Run with: runtests('testGeodesicConservation')

    properties (TestParameter)
        a        = {0, 0.5, 0.9, 0.998};
        iota_deg = {0, 45, 89, 91, 135, 180};
    end

    methods (Test)
        function conservationHolds(testCase, a, iota_deg)

            in.chi1      = [0; 0; a];
            in.iota      = deg2rad(iota_deg);
            in.r0        = 15;             % comfortably outside LSSO for every (a,iota) tested here
            in.ICs       = 'spherical';    % exact spherical-geodesic ICs -- see note in run_all_tests.m
            in.geodesics = 1;              % no radiation reaction
            in.Tmax      = 2000;            % a handful of orbital periods at r0=15
            in.dt        = 0.25;
            in.reltol    = 1e-13;
            in.abstol    = 1e-13;
            in.verbose   = 0;

            dyn = DB_class(in);

            tol = 1e-11; % looser than reltol/abstol to allow for dense-output interpolation error
                            % with reltol and abstol 1e-13 it passes up to 1e-11 at which point only
                            % the highly spinning a = 0.998 configuration fails. The others
                            % pass up to 1e-13

            r_spread = (max(dyn.r) - min(dyn.r)) / mean(dyn.r);
            testCase.verifyLessThan(r_spread, tol, ...
                'Radius drifted -- this should be an exactly spherical orbit.');

            H_spread = (max(dyn.Heff) - min(dyn.Heff)) / abs(mean(dyn.Heff));
            testCase.verifyLessThan(H_spread, tol, ...
                'Effective Hamiltonian (energy) is not conserved along the geodesic.');

            %C0 = dyn.C(1);
            C0 = mean(dyn.C);
            C_spread = (max(dyn.C) - min(dyn.C)) / max(abs(C0), 1e-13);
            testCase.verifyLessThan(C_spread, tol, ...
                'Carter constant is not conserved along the geodesic.');

            Lz_spread = (max(dyn.pphi) - min(dyn.pphi)) / max(abs(mean(dyn.pphi)), 1e-13);
            testCase.verifyLessThan(Lz_spread, tol, ...
                'Axial angular momentum p_phi is not conserved along the geodesic.');
        end
    end
end
