classdef test_radar_mutual_gain_cdf_rev1 < matlab.unittest.TestCase
    %TEST_RADAR_MUTUAL_GAIN_CDF_REV1 Unit tests for the radar coupling functions.
    %   Covers radar_mutual_gain_cdf_rev1, antenna_prepare_rev1,
    %   antenna_gain_rev1 and m1851_rect_field_rev1.
    %   Run with: runtests('test_radar_mutual_gain_cdf_rev1')

    properties (TestParameter)
        % First sidelobe levels that select each M.1851-2 Table 9 taper.
        rect_sll = struct('n0', -15, 'n1', -25, 'n2', -35, 'n3', -42, 'n4', -47)
        % First sidelobe levels that select each M.1851-2 Table 14 taper.
        parab_sll = struct('n0', -16, 'n1', -22, 'n2', -29, 'n3', -34, 'n4', -40)
        % One peak gain in each Statgain category of NTIA TM-13-489 Table 4-1.
        sg_gain = struct('low', 15, 'medium', 33.5, 'high', 50)
    end

    methods (Test)
        function rectFieldMatchesTable4ClosedForms(testCase)
            % Known value: M.1851-2 Table 4 closed forms, n = 1..3.
            mu = linspace(0.01, 20, 400)';
            f1 = (pi/2)*cos(mu)./((pi/2)^2 - mu.^2) / (2/pi);
            f2 = (pi^2./(2*mu)).*sin(mu)./(pi^2 - mu.^2) / 0.5;
            f3 = (3*pi*cos(mu)/8).*(1./((pi/2)^2 - mu.^2) - 1./((3*pi/2)^2 - mu.^2)) / (4/(3*pi));
            testCase.verifyEqual(m1851_rect_field_rev1(1, mu), f1, "AbsTol", 1e-10);
            testCase.verifyEqual(m1851_rect_field_rev1(2, mu), f2, "AbsTol", 1e-10);
            testCase.verifyEqual(m1851_rect_field_rev1(3, mu), f3, "AbsTol", 1e-10);
        end

        function rectHalfPowerAtHalfBeamwidth(testCase, rect_sll)
            % Physics: gain is ~3 dB down at half the 3 dB beamwidth
            % (M.1851-2 beamwidth factors are fits, so 0.1 dB tolerance).
            a = struct('type', 'rect', 'g_max_dBi', 30, 'az_bw_deg', 2, ...
                'az_sll_dB', rect_sll, 'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'peak');
            p = antenna_prepare_rev1(a);
            g = antenna_gain_rev1(p, 1, 0, one_state());
            testCase.verifyEqual(g - 30, -3, "AbsTol", 0.1);
        end

        function peakMaskContinuousAtBreak(testCase, rect_sll)
            % Invariant: theoretical pattern and peak mask meet at the
            % M.1851-2 Table 6 break point.
            a = struct('type', 'rect', 'g_max_dBi', 0, 'az_bw_deg', 2, ...
                'az_sll_dB', rect_sll, 'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'peak');
            p = antenna_prepare_rev1(a);
            th = p.az_cut.brk_deg + [-1e-6; 1e-6];
            g = antenna_gain_rev1(p, th, 0, two_states());
            testCase.verifyLessThan(abs(diff(g)), 0.5);
        end

        function csc2PeakAtTiltAndStartIsMinus3dB(testCase)
            % Known value: M.1851-2 Table 7, gain is 0 dB at the tilt and
            % 3 dB down at the CSC^2 start (tilt + bw/2).
            a = csc2_antenna();
            p = antenna_prepare_rev1(a);
            st = one_state();
            st.tilt = 2;
            testCase.verifyEqual(antenna_gain_rev1(p, 0, 2, st), 33.5, "AbsTol", 1e-9);
            testCase.verifyEqual(antenna_gain_rev1(p, 0, 2 + 4.8/2, st) - 33.5, -3, "AbsTol", 0.1);
            % Eq. (27) hand calculation at 10 deg, inside the CSC^2 region:
            % 20log10(sin(4.4)/sin(10)) + sinc gain at the start angle.
            g_start = 20*log10(abs(sin(pi*50.8*sind(2.4)/4.8)/(pi*50.8*sind(2.4)/4.8)));
            expected = 20*log10(sind(4.4)/sind(10)) + g_start;
            testCase.verifyEqual(antenna_gain_rev1(p, 0, 10, st) - 33.5, expected, "AbsTol", 1e-9);
        end

        function arrayPeakAndScanLoss(testCase)
            % Hand calculation: uniform array peak = Ge + 10log10(N);
            % steering to 60 deg costs the cos^2 element factor.
            p = antenna_prepare_rev1(array_antenna(16, "uniform"));
            g0 = antenna_gain_rev1(p, 0, 0, one_state());
            st = one_state();
            st.steer_az = 60;
            g60 = antenna_gain_rev1(p, 60, 0, st);
            testCase.verifyEqual(g0, 5 + 10*log10(256), "AbsTol", 1e-6);
            testCase.verifyEqual(g60 - g0, 20*log10(cosd(60)), "AbsTol", 0.01);
        end

        function taylorSidelobeLevel(testCase)
            % Independent source: taylorwin(16, 4, -25) gives -25 dB peak
            % sidelobes, so the azimuth cut's highest sidelobe is ~-25 dB.
            p = antenna_prepare_rev1(array_antenna(16, "taylor"));
            az = (0:0.05:40)';
            st = struct('mech_az', zeros(size(az)), 'tilt', zeros(size(az)), ...
                'steer_az', zeros(size(az)), 'steer_el', zeros(size(az)));
            g = antenna_gain_rev1(p, az, 0, st) - p.g_max_dBi - 20*log10(cosd(az));
            [pk, loc] = findpeaks(g);
            first_sl = pk(find(az(loc) > 3, 1));
            testCase.verifyEqual(first_sl, -25, "AbsTol", 0.7);
        end

        function couplingNeverExceedsSumOfPeaks(testCase)
            % Invariant: mutual gain <= Gmax1 + Gmax2.
            [s1, s2] = two_radars(75, 90);
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), 33.5 + 34 + 1e-9);
        end

        function reciprocalWhenTxEqualsRx(testCase)
            % Invariant: g12 == g21 when each radar uses one pattern.
            [s1, s2] = two_radars(75, 90);
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1));
            testCase.verifyEqual(out.g12_dBi, out.g21_dBi, "AbsTol", 1e-9);
        end

        function sameSeedSameResult(testCase)
            % Reproducibility.
            [s1, s2] = two_radars(75, 90.3);
            a = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(7));
            b = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(7));
            c = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(8));
            testCase.verifyEqual(a.g12_dBi, b.g12_dBi);
            testCase.verifyNotEqual(a.g12_dBi, c.g12_dBi);
        end

        function timeSeriesConvergesToRandom(testCase)
            % Method check (M.2069-0): for incommensurate rates, a long
            % time series matches the random-sampling CDF.
            [s1, s2] = two_radars(75, 90.37);
            o = quick_opts(3);
            o.n_samples = 2e6;
            r = radar_mutual_gain_cdf_rev1(s1, s2, o);
            o.mode = "time_series";
            o.duration_s = 7200;
            o.n_trials = 2;
            t = radar_mutual_gain_cdf_rev1(s1, s2, o);
            p = [99 99.9];
            testCase.verifyEqual(prctile(t.g12_dBi, p), prctile(r.g12_dBi, p), "AbsTol", 1.5);
        end

        function pathElevationIsUsed(testCase)
            % Location-agnostic geometry: the path elevation input changes
            % the gain each antenna presents along the link.
            [s1, s2] = two_radars(75, 90.3);
            o = quick_opts(2);
            flat = radar_mutual_gain_cdf_rev1(s1, s2, o);
            o.el12_deg = 15;
            o.el21_deg = -15;
            tilted = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(tilted.geom.el12_deg, 15);
            testCase.verifyEqual(flat.geom.az21_deg, 180);
            testCase.verifyNotEqual(max(tilted.g12_dBi), max(flat.g12_dBi));
        end

        function commensurateRatesWarn(testCase)
            % 75 and 90 deg/s are 5:6 and repeat every 24 s.
            [s1, s2] = two_radars(75, 90);
            testCase.verifyWarning(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)), ...
                "radar_coupling:commensurateRates");
        end

        function equalRatesWarnWhenUnsynchronized(testCase)
            % Unsynchronized equal rates run as MACD (unknown phase) with
            % a warning that points to the synchronized option.
            [s1, s2] = two_radars(75, 75);
            testCase.verifyWarning(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)), ...
                "radar_coupling:equalRates");
        end

        function syncMaxIsMainbeamPlusBacklobe(testCase)
            % Known value: synchronized at offset 0, each main beam meets
            % the other antenna's back lobe. Statgain back lobe for
            % 22 <= G < 48 is 11 - G/2 (NTIA TM-13-489 Table 4-1), so
            % 34/34 -> 34 + (11 - 17) = 28 dBi and
            % 34/46 -> max(34 + (11 - 23), 46 + (11 - 17)) = 40 dBi.
            [s1, s2] = sg_radars(34, 34, 75);
            o = sync_opts(1, 0);
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(out.sync.max_g12_dBi, 28, "AbsTol", 1e-9);
            testCase.verifyEqual(out.sync.max_g21_dBi, 28, "AbsTol", 1e-9);
            [s1, s2] = sg_radars(34, 46, 75);
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(out.sync.max_g12_dBi, 40, "AbsTol", 1e-9);
        end

        function syncOffset180IsMainbeamToMainbeam(testCase)
            % Known value: offset 180 aligns both main beams once per turn.
            [s1, s2] = sg_radars(34, 46, 75);
            out = radar_mutual_gain_cdf_rev1(s1, s2, sync_opts(1, 180));
            testCase.verifyEqual(out.sync.max_g12_dBi, 80, "AbsTol", 1e-9);
            testCase.verifyEqual(out.sync.az1_at_max12_deg, 0, "AbsTol", 1e-9);
        end

        function syncRandomSamplesMatchSweep(testCase)
            % Invariant: synchronized samples never exceed the sweep maximum
            % and reach it closely; unsynchronized (MACD) samples of the
            % same pair go far above it.
            [s1, s2] = sg_radars(34, 34, 75);
            sync = radar_mutual_gain_cdf_rev1(s1, s2, sync_opts(4, 0));
            testCase.verifyLessThanOrEqual(max(sync.g12_dBi), sync.sync.max_g12_dBi + 1e-9);
            testCase.verifyGreaterThan(max(sync.g12_dBi), sync.sync.max_g12_dBi - 0.05);
            o = quick_opts(4);
            unsync = testCase.verifyWarning(@() radar_mutual_gain_cdf_rev1(s1, s2, o), ...
                "radar_coupling:equalRates");
            testCase.verifyEmpty(unsync.sync);
            testCase.verifyGreaterThan(prctile(unsync.g12_dBi, 99.9), sync.sync.max_g12_dBi + 10);
        end

        function syncTimeSeriesMatchesSweep(testCase)
            % Time series of synchronized antennas stays at or below the
            % sweep maximum and reaches it within the time-step loss.
            [s1, s2] = sg_radars(34, 34, 75);
            o = sync_opts(5, 0);
            o.mode = "time_series";
            o.duration_s = 30;
            o.n_trials = 2;
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), out.sync.max_g12_dBi + 1e-9);
            testCase.verifyGreaterThan(max(out.g12_dBi), out.sync.max_g12_dBi - 0.1);
        end

        function f1245MatchesHandValues(testCase)
            % Known values: Rec. ITU-R F.1245 as coded in
            % AntennaModel_ITU_1245_1.cs, D/lambda = 10^((G - 7.7)/20).
            % D/lambda <= 100 branch (34 dBi, D/lambda = 20.65):
            dl = 10^((34 - 7.7)/20);
            p = antenna_prepare_rev1(struct('type', 'f1245', 'g_max_dBi', 34));
            g = antenna_gain_rev1(p, [1; 10; 180], 0, states(3));
            expected = [34 - 2.5e-3*dl^2; 39 - 5*log10(dl) - 25; -3 - 5*log10(dl)];
            testCase.verifyEqual(g, expected, "AbsTol", 1e-9);
            testCase.verifyEqual(g(3), -9.575, "AbsTol", 1e-9);   % paper MAG = 34 + this
            % D/lambda > 100 branch (50 dBi): G1 plateau, 29 - 25log, -13.
            dl = 10^((50 - 7.7)/20);
            g1 = 2 + 15*log10(dl);
            p = antenna_prepare_rev1(struct('type', 'f1245', 'g_max_dBi', 50));
            testCase.assumeGreaterThan(p.f.phi_r, p.f.phi_m);
            phi_mid = mean([p.f.phi_m, p.f.phi_r]);
            g = antenna_gain_rev1(p, [phi_mid; 30; 90], 0, states(3));
            testCase.verifyEqual(g, [g1; 29 - 25*log10(30); -13], "AbsTol", 1e-9);
        end

        function f1245ReproducesPaperSyncMAG(testCase)
            % Known values: 2.7 GHz link budget working paper Table 4.6-1
            % MAGs (ASR 34 dBi, WSR-88D 46 dBi, F.1245) are the synchronized
            % mainbeam-plus-backlobe values 24.425, 36.425 and 33.425 dBi.
            pairs = [34 34 24.425; 34 46 36.425; 46 46 33.425];
            for k = 1:3
                [s1, s2] = fg_radars(pairs(k, 1), pairs(k, 2), 75);
                out = radar_mutual_gain_cdf_rev1(s1, s2, sync_opts(1, 0));
                testCase.verifyEqual(out.sync.max_g12_dBi, pairs(k, 3), "AbsTol", 1e-3);
            end
        end

        function f1245MacdMatchesCsharp(testCase)
            % Independent source: the C# MACD method (independent uniform
            % pointing) gives 24.35 dBi at 99% for two 34 dBi F.1245
            % antennas (crosscheck_csharp_rev1, continuous reference).
            [s1, s2] = fg_radars(34, 34, 75);
            s2.scan.rate_deg_s = 90.37;
            o = quick_opts(6);
            o.n_samples = 1e6;
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(prctile(out.g12_dBi, 99), 24.35, "AbsTol", 0.1);
        end

        function f1245RangeErrors(testCase)
            % Gmax must exceed G1 = 2 + 15log10(D/lambda).
            a = struct('type', 'f1245', 'g_max_dBi', 20, 'd_over_lambda', 100);
            testCase.verifyError(@() antenna_prepare_rev1(a), "radar_coupling:f1245Range");
        end

        function syncInputErrors(testCase)
            [s1, s2] = sg_radars(34, 34, 75);
            s2.scan.rate_deg_s = 90;
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, sync_opts(1, 0)), ...
                "radar_coupling:syncNeedsEqualRates");
            [s1, s2] = sg_radars(34, 34, 75);
            s2.scan.rate_jitter_deg_s = 0.5;
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, sync_opts(1, 0)), ...
                "radar_coupling:syncJitter");
            o = sync_opts(1, 0);
            o.scan_sync = "locked";
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, o), ...
                "radar_coupling:badScanSync");
        end

        function missingSeedErrors(testCase)
            [s1, s2] = two_radars(75, 90.3);
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, struct()), ...
                "radar_coupling:missingSeed");
        end

        function badScanTypeErrors(testCase)
            [s1, s2] = two_radars(75, 90.3);
            s2.scan.type = "spiral";
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)), ...
                "radar_coupling:badScanType");
        end

        function outputContract(testCase)
            [s1, s2] = two_radars(75, 90.3);
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1));
            testCase.verifySize(out.g12_dBi, [2e5 1]);
            testCase.verifyClass(out.cdf12, "table");
            testCase.verifyTrue(all(diff(out.cdf12.cdf) >= 0));
            testCase.verifyEqual(out.cdf12.cdf(end), 1, "AbsTol", 1e-12);
            testCase.verifyTrue(all(diff(out.pct12.gain_dBi) >= 0));
        end

        % ---- Statgain (default): known values from NTIA TM-13-489 ------

        function statgainIsDefaultType(testCase)
            % A model with no type field uses Statgain.
            p = antenna_prepare_rev1(struct('g_max_dBi', 33.5));
            testCase.verifyEqual(p.type, "statgain");
        end

        function statgainMediumGainKnownValues(testCase)
            % Table 4-1, 22 <= Gmax < 48 (Gmax = 33.5 dBi):
            %   0 deg: Gmax; 5 deg (theta_m < 5 < theta_r): 0.75*Gmax - 7;
            %   10 deg: 53 - Gmax/2 - 25log10(10); 60 deg (> 48): 11 - Gmax/2.
            G = 33.5;
            p = antenna_prepare_rev1(struct('g_max_dBi', G));
            th = [0; 5; 10; 60];
            g = antenna_gain_rev1(p, th, 0, states(numel(th)));
            expected = [G; 0.75*G - 7; 53 - G/2 - 25; 11 - G/2];
            testCase.verifyEqual(g, expected, "AbsTol", 1e-9);
            testCase.verifyEqual(p.sg.theta_m, 50*sqrt(0.25*G + 7)/10^(G/20), "RelTol", 1e-12);
            testCase.verifyEqual(p.sg.theta_r, 250/10^(G/20), "RelTol", 1e-12);
        end

        function statgainHighGainKnownValues(testCase)
            % Table 4-1, Gmax >= 48 (Gmax = 50 dBi): 10 deg gives
            % 29 - 25log10(10) = 4 dBi; 100 deg gives -13 dBi;
            % theta_r1 = 27.466*10^(-0.3*Gmax/10).
            p = antenna_prepare_rev1(struct('g_max_dBi', 50));
            g = antenna_gain_rev1(p, [10; 100], 0, states(2));
            testCase.verifyEqual(g, [4; -13], "AbsTol", 1e-9);
            testCase.verifyEqual(p.sg.theta_r, 27.466*10^(-1.5), "RelTol", 1e-12);
        end

        function statgainLowGainKnownValues(testCase)
            % Table 4-1, 10 <= Gmax < 22 (Gmax = 15 dBi): theta_r3 =
            % 250/10^(15/20) = 44.5 deg and theta_b3 = 131.8257*10^(-0.3)
            % = 66.07 deg. 30 deg is on the 0.75*Gmax - 7 plateau, 50 deg
            % on 53 - Gmax/2 - 25log10(theta), and 90 deg is 0 dBi.
            p = antenna_prepare_rev1(struct('g_max_dBi', 15));
            testCase.verifyEqual(p.sg.theta_b, 131.8257*10^(-0.3), "RelTol", 1e-12);
            g = antenna_gain_rev1(p, [30; 50; 90], 0, states(3));
            testCase.verifyEqual(g, [0.75*15 - 7; 53 - 7.5 - 25*log10(50); 0], "AbsTol", 1e-9);
        end

        function statgainContinuousAtMainBeamEdge(testCase, sg_gain)
            % Invariant: theta_m is where the main beam meets 0.75*Gmax - 7.
            p = antenna_prepare_rev1(struct('g_max_dBi', sg_gain));
            th = p.sg.theta_m + [-1e-7; 1e-7];
            g = antenna_gain_rev1(p, th, 0, states(2));
            testCase.verifyEqual(g(1), g(2), "AbsTol", 1e-4);
        end

        function statgainUsesOffAxisAngle(testCase)
            % Rotationally symmetric: 5 deg in elevation equals 5 deg in
            % azimuth; a 5 deg tilt points the beam at 5 deg elevation.
            p = antenna_prepare_rev1(struct('g_max_dBi', 33.5));
            g_az = antenna_gain_rev1(p, 5, 0, one_state());
            g_el = antenna_gain_rev1(p, 0, 5, one_state());
            st = one_state();
            st.tilt = 5;
            testCase.verifyEqual(g_el, g_az, "AbsTol", 1e-9);
            testCase.verifyEqual(antenna_gain_rev1(p, 0, 5, st), 33.5, "AbsTol", 1e-9);
        end

        function statgainMonotoneNonIncreasing(testCase, sg_gain)
            % Invariant: Statgain does not rise with off-axis angle, except
            % the steps of up to ~0.05 dB built into the published Table 4-1
            % constants at theta_r and theta_b (for example, at 48 deg the
            % 25log10 segment ends at 10.97 - Gmax/2, below 11 - Gmax/2).
            p = antenna_prepare_rev1(struct('g_max_dBi', sg_gain));
            th = (0:0.01:180)';
            g = antenna_gain_rev1(p, th, 0, states(numel(th)));
            testCase.verifyLessThanOrEqual(max(diff(g)), 0.06);
        end

        function statgainBelow10dBiErrors(testCase)
            testCase.verifyError(@() antenna_prepare_rev1(struct('g_max_dBi', 8)), ...
                "radar_coupling:statgainRange");
        end

        function statgainRadarsSimulate(testCase)
            % End to end with default antennas: coupling stays at or
            % below the sum of peaks, and reaches near it in the tail.
            s1 = struct('tx', struct('g_max_dBi', 33.5), ...
                'scan', struct('type', 'mech_rotating', 'rate_deg_s', 75));
            s2 = struct('tx', struct('g_max_dBi', 34), ...
                'scan', struct('type', 'mech_rotating', 'rate_deg_s', 90.37));
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(12));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), 67.5 + 1e-9);
            testCase.verifyGreaterThan(max(out.g12_dBi), 67.5 - 6);
            % Each antenna's lowest Table 4-1 value is the far-sidelobe
            % segment at 48 deg, 53 - Gmax/2 - 25log10(48), just below the
            % 11 - Gmax/2 back lobe. Sampling rarely puts both antennas
            % at that exact angle at once, so the sampled minimum lies
            % between the joint minimum and the sum of the back lobes.
            lo = @(G) min(11 - G/2, 53 - G/2 - 25*log10(48));
            g_min = min(out.g12_dBi);
            testCase.verifyGreaterThanOrEqual(g_min, lo(33.5) + lo(34) - 1e-9);
            testCase.verifyLessThanOrEqual(g_min, (11 - 33.5/2) + (11 - 34/2) + 1e-9);
        end

        % ---- Antenna patterns: known values from M.1851-2 ----------------

        function parabolicHalfPowerAtHalfBeamwidth(testCase, parab_sll)
            % Physics: ~3 dB down at half the 3 dB beamwidth. The Table 11
            % beamwidth factor for n = 4 gives -3.28 dB, so 0.3 dB tolerance.
            b = struct('type', 'parabolic', 'g_max_dBi', 30, 'az_bw_deg', 2, ...
                'az_sll_dB', parab_sll, 'envelope', 'peak');
            p = antenna_prepare_rev1(b);
            testCase.verifyEqual(antenna_gain_rev1(p, 0, 0, one_state()), 30, "AbsTol", 1e-9);
            testCase.verifyEqual(antenna_gain_rev1(p, 1, 0, one_state()) - 30, -3, "AbsTol", 0.3);
        end

        function rectPeakAndAverageMaskTable6(testCase)
            % Known value, M.1851-2 Table 6 (uniform): at 3 beamwidths the
            % peak mask is -8.584*ln(2.876*3); the average mask adds -3.72.
            a = struct('type', 'rect', 'g_max_dBi', 0, 'az_bw_deg', 2, 'az_sll_dB', -15, ...
                'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'peak');
            pk = antenna_gain_rev1(antenna_prepare_rev1(a), 6, 0, one_state());
            a.envelope = 'average';
            av = antenna_gain_rev1(antenna_prepare_rev1(a), 6, 0, one_state());
            testCase.verifyEqual(pk, -8.584*log(2.876*3), "AbsTol", 1e-9);
            testCase.verifyEqual(av, -8.584*log(2.876*3) - 3.72, "AbsTol", 1e-9);
        end

        function parabolicPeakAndAverageMaskTable13(testCase)
            % Known value, M.1851-2 Table 13 (n = 0): at 2 beamwidths the
            % peak mask is -28.9*log10(2) - 11.9; average is 4 dB lower.
            b = struct('type', 'parabolic', 'g_max_dBi', 0, 'az_bw_deg', 2, ...
                'az_sll_dB', -16, 'envelope', 'peak');
            pk = antenna_gain_rev1(antenna_prepare_rev1(b), 4, 0, one_state());
            b.envelope = 'average';
            av = antenna_gain_rev1(antenna_prepare_rev1(b), 4, 0, one_state());
            testCase.verifyEqual(pk, -28.9*log10(2) - 11.9, "AbsTol", 1e-9);
            testCase.verifyEqual(av, -28.9*log10(2) - 11.9 - 4, "AbsTol", 1e-9);
        end

        function theoreticalPatternHitsFloorAtNull(testCase)
            % Known value: the uniform taper's first null (mu = pi) is at
            % asind(bw/50.8); the theoretical envelope clamps to the -30 dB
            % Table 6 floor there.
            a = struct('type', 'rect', 'g_max_dBi', 0, 'az_bw_deg', 2, 'az_sll_dB', -15, ...
                'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'theoretical');
            g = antenna_gain_rev1(antenna_prepare_rev1(a), asind(2/50.8), 0, one_state());
            testCase.verifyEqual(g, -30, "AbsTol", 1e-9);
        end

        function averageNeverAbovePeak(testCase)
            % Invariant: the average mask never exceeds the peak mask.
            a = struct('type', 'rect', 'g_max_dBi', 0, 'az_bw_deg', 2, 'az_sll_dB', -25, ...
                'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'peak');
            az = (0:0.1:180)';
            st = struct('mech_az', zeros(size(az)), 'tilt', zeros(size(az)), ...
                'steer_az', zeros(size(az)), 'steer_el', zeros(size(az)));
            pk = antenna_gain_rev1(antenna_prepare_rev1(a), az, 0, st);
            a.envelope = 'average';
            av = antenna_gain_rev1(antenna_prepare_rev1(a), az, 0, st);
            testCase.verifyTrue(all(av <= pk + 1e-9));
        end

        function measuredPatternInterpolates(testCase)
            % Known value: linear interpolation of a digitized cut.
            m = struct('type', 'measured', 'g_max_dBi', 30, ...
                'az_deg', [-180 -10 0 10 180], 'az_gain_dB', [-40 -20 0 -20 -40], 'floor_dB', -50);
            g = antenna_gain_rev1(antenna_prepare_rev1(m), 5, 0, one_state());
            testCase.verifyEqual(g, 30 - 10, "AbsTol", 1e-9);
        end

        function arrayPatternSymmetric(testCase)
            % Invariant: a broadside-steered array is symmetric in azimuth.
            p = antenna_prepare_rev1(array_antenna(16, "taylor"));
            st = two_states();
            g = antenna_gain_rev1(p, [12; -12], 0, st);
            testCase.verifyEqual(g(1), g(2), "AbsTol", 1e-9);
        end

        function arrayBackHemisphereUsesFrontToBack(testCase)
            % Known value: behind the face, gain = peak + front_back_dB.
            a = array_antenna(16, "uniform");
            a.front_back_dB = -45;                  % differs from floor_dB = -60
            p = antenna_prepare_rev1(a);
            g = antenna_gain_rev1(p, 180, 0, one_state());
            testCase.verifyEqual(g, p.g_max_dBi - 45, "AbsTol", 1e-9);
        end

        function arrayGratingLobeWarns(testCase)
            % M.1851-2 Section 6.1: spacing above lambda/2 risks grating lobes.
            a = array_antenna(16, "uniform");
            a.dx_m = 0.7*physconst('LightSpeed')/3e9;
            testCase.verifyWarning(@() antenna_prepare_rev1(a), "radar_coupling:gratingLobes");
        end

        function badSidelobeLevelsError(testCase)
            a = struct('type', 'rect', 'g_max_dBi', 0, 'az_bw_deg', 2, 'az_sll_dB', -10, ...
                'el_bw_deg', 5, 'el_sll_dB', -15, 'envelope', 'peak');
            testCase.verifyError(@() antenna_prepare_rev1(a), "radar_coupling:badSLL");
            b = struct('type', 'parabolic', 'g_max_dBi', 0, 'az_bw_deg', 2, ...
                'az_sll_dB', -12, 'envelope', 'peak');
            testCase.verifyError(@() antenna_prepare_rev1(b), "radar_coupling:badSLL");
        end

        function csc2StartBelowHorizonErrors(testCase)
            p = antenna_prepare_rev1(csc2_antenna());
            st = one_state();
            st.tilt = -5;
            testCase.verifyError(@() antenna_gain_rev1(p, 0, 0, st), "radar_coupling:csc2Start");
        end

        function arrayDefaultsGive41dBiArray(testCase)
            % Defaults: 41 dBi target (M.1464-2 Radar F1), lambda/2 at
            % 2900 MHz, unit-cell element 10log10(pi), Taylor -25 dB, and
            % floor/back level at the TM-13-489 back level 11 - G/2.
            p = antenna_prepare_rev1(struct('type', 'array_planar'));
            testCase.verifyGreaterThanOrEqual(p.g_max_dBi, 41);
            testCase.verifyLessThan(p.g_max_dBi, 41.3);       % smallest square that reaches 41
            lam = physconst('LightSpeed')/2.9e9;
            testCase.verifyEqual(p.dx_m, lam/2, "RelTol", 1e-12);
            testCase.verifyEqual(p.elem_gain_dBi, 10*log10(pi), "AbsTol", 1e-12);
            testCase.verifyEqual(p.taper, "taylor");
            testCase.verifyEqual(p.floor_dB, (11 - p.g_max_dBi/2) - p.g_max_dBi, "AbsTol", 1e-12);
            testCase.verifyEqual(p.front_back_dB, p.floor_dB);
            % One element fewer per side would fall short of 41 dBi.
            q = antenna_prepare_rev1(struct('type', 'array_planar', 'nx', p.nx - 1));
            testCase.verifyLessThan(q.g_max_dBi, 41);
            % Broadside gain equals g_max; behind the face is the back level.
            g = antenna_gain_rev1(p, [0; 180], 0, states(2));
            testCase.verifyEqual(g(1), p.g_max_dBi, "AbsTol", 1e-6);
            testCase.verifyEqual(g(2), 11 - p.g_max_dBi/2, "AbsTol", 1e-9);
        end

        function arrayScanDefaultsRun(testCase)
            % Minimal inputs: default array on a 4-face mount and on a
            % rotating mount, both horizon only. Coupling stays below the
            % sum of peaks; the 4-face case reaches near peak against a
            % flat antenna (beam azimuth is random over 360 deg).
            [~, flat] = pencil_vs_flat();
            arr = struct('type', 'array_planar');
            gmax = antenna_prepare_rev1(arr).g_max_dBi;
            mf = struct('tx', arr, 'scan', struct('type', 'array_multiface'));
            out = radar_mutual_gain_cdf_rev1(flat, mf, quick_opts(21));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), gmax + 1e-9);
            testCase.verifyGreaterThan(max(out.g12_dBi), gmax - 4);
            rot = struct('tx', arr, 'scan', struct('type', 'mech_rotating', 'rate_deg_s', 60));
            out = radar_mutual_gain_cdf_rev1(flat, rot, quick_opts(22));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), gmax + 1e-9);
            rot.scan.beam_mode = 'random';
            out = radar_mutual_gain_cdf_rev1(flat, rot, quick_opts(23));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), gmax + 1e-9);
        end

        function arraySizeConflictErrors(testCase)
            a = struct('type', 'array_planar', 'nx', 20, 'g_max_dBi', 41);
            testCase.verifyError(@() antenna_prepare_rev1(a), "radar_coupling:arraySizeConflict");
        end

        function missingAntennaFieldErrors(testCase)
            a = rmfield(csc2_antenna(), 'az_bw_deg');
            testCase.verifyError(@() antenna_prepare_rev1(a), "radar_coupling:missingField");
        end

        % ---- Scan and simulation behavior --------------------------------

        function volumeScanTimeWeighting(testCase)
            % Known value: steps at 10 and 30 deg/s spend 36 s and 12 s per
            % turn, so tilt 0 holds 75% of the time. Against a flat
            % antenna, the main beam (within -3 dB) points at the other
            % radar for 0.75 * 2/360 of the samples.
            [s1, s2] = pencil_vs_flat();
            s1.scan = struct('type', 'mech_rotating', 'rate_deg_s', [10 30], 'tilt_deg', [0 20]);
            o = quick_opts(11);
            o.n_samples = 1e6;
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            p = mean(out.g12_dBi >= 40 - 3);
            expected = 0.75*2/360;
            testCase.verifyEqual(p, expected, "AbsTol", 4*sqrt(expected/o.n_samples) + 2e-4);
        end

        function eventStatisticsKnownValue(testCase)
            % Known value: a 2 deg beam rotating at 36 deg/s passes the
            % other radar once every 10 s for about 2/36 s.
            [s1, s2] = pencil_vs_flat();
            o = struct('seed', 4, 'mode', 'time_series', 'duration_s', 100, ...
                'n_trials', 1, 'event_thresh_dBi', 40 - 3);
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(out.events12.n_events, 10, "AbsTol", 1);
            testCase.verifyEqual(out.events12.mean_duration_s, 2/36, "AbsTol", 2*out.opts.dt_s);
        end

        function sectorScanKnownValues(testCase)
            % Known values: a 2 deg beam sweeping a 20 deg sector at 36
            % deg/s passes the other radar once per 20/36 s sweep, so 180
            % times in 100 s, each pass about 2/36 s long. Random mode puts
            % the beam within +/-1 deg of the path 2/20 of the time, 18
            % times the full-rotation share.
            [s1, s2] = pencil_vs_flat();
            s1.scan.sector_deg = 20;
            o = struct('seed', 4, 'mode', 'time_series', 'duration_s', 100, ...
                'n_trials', 1, 'event_thresh_dBi', 40 - 3);
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(out.events12.n_events, 180, "AbsTol", 1);
            testCase.verifyEqual(out.events12.mean_duration_s, 2/36, "AbsTol", 2*out.opts.dt_s);
            r = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(3));
            [f1, f2] = pencil_vs_flat();
            f = radar_mutual_gain_cdf_rev1(f1, f2, quick_opts(3));
            p_sec = mean(r.g12_dBi >= max(r.g12_dBi) - 3);
            p_full = mean(f.g12_dBi >= max(f.g12_dBi) - 3);
            testCase.verifyEqual(p_sec/p_full, 18, "RelTol", 0.1);
            s1.scan.sector_deg = 0;
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(3)), "radar_coupling:badSector");
        end

        function rpmEqualsDegPerSecond(testCase)
            % rate_rpm = 12.5 is 75 deg/s.
            [s1, s2] = two_radars(75, 90.3);
            a = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(6));
            s1.scan = rmfield(s1.scan, 'rate_deg_s');
            s1.scan.rate_rpm = 12.5;
            b = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(6));
            testCase.verifyEqual(b.g12_dBi, a.g12_dBi);
        end

        function beatAndJointPeriods(testCase)
            % Known values: 360/|75-90| = 24 s; 75:90 = 5:6, repeat 360*5/75 = 24 s.
            [s1, s2] = two_radars(75, 90);
            out = evalc_quiet(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)));
            testCase.verifyEqual(out.geom.beat_period_s, 24, "AbsTol", 1e-9);
            testCase.verifyEqual(out.geom.joint_period_s, 24, "AbsTol", 1e-9);
            testCase.verifyEqual(out.geom.ratio, [5 6]);
            [s1, s2] = two_radars(75, 90.37);
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1));
            testCase.verifyEqual(out.geom.joint_period_s, Inf);
        end

        function equalRatesRunInTimeSeries(testCase)
            % time_series handles the case random mode rejects.
            [s1, s2] = two_radars(75, 75);
            o = struct('seed', 2, 'mode', 'time_series', 'duration_s', 30, 'n_trials', 2);
            out = radar_mutual_gain_cdf_rev1(s1, s2, o);
            testCase.verifyEqual(numel(out.g12_dBi), 2*numel(out.ts.t_s));
            testCase.verifyLessThanOrEqual(max(out.g12_dBi), 33.5 + 34 + 1e-9);
        end

        function separateTxRxPatternsBreakReciprocity(testCase)
            % Transmit and receive tapers differ (M.1464-2 Radar F1), so
            % g12 and g21 differ, but both stay below the sum of peaks.
            s = two_array_radars();
            s{2}.rx = s{2}.tx;
            s{2}.tx.taylor_sll_dB = -17;
            o = quick_opts(9);
            out = evalc_quiet(@() radar_mutual_gain_cdf_rev1(s{1}, s{2}, o));
            gmax2 = antenna_prepare_rev1(s{2}.rx).g_max_dBi;
            testCase.verifyNotEqual(out.g12_dBi, out.g21_dBi);
            testCase.verifyLessThanOrEqual(max([out.g12_dBi; out.g21_dBi]), 33.5 + gmax2 + 0.5);
        end

        function multifaceRandomAndSequentialRun(testCase)
            % Four-face staring array in both beam modes. In sequential
            % mode the beam sweeps 360 deg every 4 s, so against a flat
            % antenna it must reach the peak gain less the face scan loss
            % (steer 45 deg: cos^2 element, -3 dB) and the dwell step
            % (4.5 deg steps of a ~4 deg beam).
            [s1, ~] = two_radars(75, 90.3);
            mf = struct('name', 'MF', 'tx', array_antenna(24, "taylor"), ...
                'scan', struct('type', 'array_multiface', 'n_faces', 4, 'face_az0_deg', 45, ...
                'beam_mode', 'random', 'beam_el_lim_deg', [0 10], 'dwell_s', 0.05));
            gmax = antenna_prepare_rev1(mf.tx).g_max_dBi;
            r = radar_mutual_gain_cdf_rev1(s1, mf, quick_opts(3));
            testCase.verifyLessThanOrEqual(max(r.g12_dBi), 33.5 + gmax + 1e-9);
            [~, flat] = pencil_vs_flat();
            mf.scan.beam_mode = 'sequential';
            mf.scan.rate_deg_s = 90;
            mf.scan.beam_el_lim_deg = [0 0];
            o = struct('seed', 3, 'mode', 'time_series', 'duration_s', 20, 'n_trials', 1);
            t = radar_mutual_gain_cdf_rev1(flat, mf, o);
            testCase.verifyLessThanOrEqual(max(t.g12_dBi), gmax + 1e-9);
            testCase.verifyGreaterThan(max(t.g12_dBi), gmax - 10);
        end

        function timeSeriesTrialsStartAnywhereInScanCycle(testCase)
            % Red team F2, fault 10: each time-series trial must start at
            % a random point of the scan cycle, not at its first step.
            % Radar 1 steps through six tilts (0 to 15 deg), one 5 s
            % revolution each; the other radar sits 9 deg above it. Trials
            % are 5 s long, so a trial sees the tilt-9 main beam only if
            % it starts near that step: with random starts, P(none of 40)
            % = (5/6)^40, about 0.07%. Starting every trial at tilt 0
            % caps radar 1's gain at its 9 deg off-axis level.
            s1 = struct('name', 'A', 'tx', struct('g_max_dBi', 46), 'scan', struct( ...
                'type', 'mech_rotating', 'rate_deg_s', 72 + zeros(6, 1), 'tilt_deg', (0:3:15)'));
            s2 = struct('name', 'B', 'tx', struct('g_max_dBi', 15), ...
                'scan', struct('type', 'mech_rotating', 'rate_deg_s', 75));
            o = struct('seed', 21, 'mode', 'time_series', 'duration_s', 5, 'n_trials', 40, ...
                'el12_deg', 9, 'el21_deg', -9);
            r = radar_mutual_gain_cdf_rev1(s1, s2, o);
            p1 = antenna_prepare_rev1(s1.tx);
            g1_at_9 = antenna_gain_rev1(p1, 0, 9, struct('mech_az', 0, 'tilt', 0, 'steer_az', 0, 'steer_el', 0));
            best_if_stuck = g1_at_9 + 15;                 % radar 1 at tilt 0, radar 2 at its peak
            testCase.verifyGreaterThan(max(r.g12_dBi), best_if_stuck + 10);
        end

        function multifaceUsesNearestFace(testCase)
            % Independent check of face selection: with faces at 0, 90,
            % 180, 270 deg and random beam azimuths, the fraction of time
            % the array's gain toward the other radar (at 180 deg) is
            % within 3 dB of peak must match a direct sweep of beam
            % azimuth with the nearest face chosen in the test.
            [~, flat] = pencil_vs_flat();
            arr = array_antenna(24, "taylor");
            mf = struct('tx', arr, 'scan', struct('type', 'array_multiface', 'n_faces', 4, ...
                'face_az0_deg', 0, 'beam_mode', 'random', 'beam_el_lim_deg', [0 0], 'dwell_s', 0.05));
            o = quick_opts(13);
            o.n_samples = 4e5;
            out = radar_mutual_gain_cdf_rev1(flat, mf, o);
            p = antenna_prepare_rev1(arr);
            beam = (0:0.005:360-0.005)';
            face = mod(90*round(beam/90), 360);
            st = struct('mech_az', face, 'tilt', zeros(size(beam)), ...
                'steer_az', wrapTo180(beam - face), 'steer_el', zeros(size(beam)));
            g = antenna_gain_rev1(p, 180, 0, st);
            expected = mean(g >= p.g_max_dBi - 3);
            sim = mean(out.g12_dBi >= p.g_max_dBi - 3);
            testCase.verifyEqual(sim, expected, "AbsTol", 4*sqrt(expected/o.n_samples) + 5e-4);
        end

        function commensurateRatesCanMissAlignment(testCase)
            % M.2069-0 rotation-mismatch effect: radar A at 75 deg/s and a
            % beam sweeping at 90 deg/s (5:6) revisit only 5 relative
            % phases, so a short time series can miss main-beam alignment
            % that random mode samples.
            [s1, ~] = two_radars(75, 90.3);
            mf = struct('name', 'MF', 'tx', array_antenna(24, "taylor"), ...
                'scan', struct('type', 'array_multiface', 'n_faces', 4, 'face_az0_deg', 45, ...
                'beam_mode', 'sequential', 'rate_deg_s', 90, 'beam_el_lim_deg', [0 0], 'dwell_s', 0.05));
            o = struct('seed', 3, 'mode', 'time_series', 'duration_s', 60, 'n_trials', 1);
            t = radar_mutual_gain_cdf_rev1(s1, mf, o);
            gmax = antenna_prepare_rev1(mf.tx).g_max_dBi;
            testCase.verifyLessThan(max(t.g12_dBi), 33.5 + gmax - 20);
        end

        function multifaceNeedsArray(testCase)
            [s1, s2] = two_radars(75, 90.3);
            s2.scan = struct('type', 'array_multiface', 'n_faces', 4, 'face_az0_deg', 0, ...
                'beam_el_lim_deg', [0 10], 'dwell_s', 0.05);
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)), ...
                "radar_coupling:scanAntennaMismatch");
        end

        function badJitterErrors(testCase)
            [s1, s2] = two_radars(75, 90.3);
            s1.scan.rate_jitter_deg_s = 80;
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1)), ...
                "radar_coupling:badJitter");
        end

        function badModeErrors(testCase)
            [s1, s2] = two_radars(75, 90.3);
            o = quick_opts(1);
            o.mode = "monte";
            testCase.verifyError(@() radar_mutual_gain_cdf_rev1(s1, s2, o), "radar_coupling:badMode");
        end

        function tailCountMatchesSampleSize(testCase)
            % Output contract: n_tail = n*(1 - p/100).
            [s1, s2] = two_radars(75, 90.3);
            out = radar_mutual_gain_cdf_rev1(s1, s2, quick_opts(1));
            testCase.verifyEqual(out.pct12.n_tail, 2e5*(1 - out.pct12.pct/100), "RelTol", 1e-12);
        end
    end

    methods (Test, TestTags = {'Parallel'})
        function parallelMatchesSerial(testCase)
            % parfor over chunks with RandStream substreams must reproduce
            % the serial result exactly.
            testCase.assumeNotEmpty(ver('parallel'));
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool(2);
            end
            parfevalOnAll(pool, @addpath, 0, fileparts(which('radar_mutual_gain_cdf_rev1')));
            s = two_array_radars();
            o = quick_opts(5);
            o.parallel_flag = 0;
            a = radar_mutual_gain_cdf_rev1(s{1}, s{2}, o);
            o.parallel_flag = 1;
            b = radar_mutual_gain_cdf_rev1(s{1}, s{2}, o);
            testCase.verifyEqual(b.g12_dBi, a.g12_dBi);
            o.mode = "time_series";
            o.duration_s = 120;
            o.n_trials = 3;
            o.parallel_flag = 0;
            c = radar_mutual_gain_cdf_rev1(s{1}, s{2}, o);
            o.parallel_flag = 1;
            d = radar_mutual_gain_cdf_rev1(s{1}, s{2}, o);
            testCase.verifyEqual(d.g12_dBi, c.g12_dBi);
        end
    end
end

% ---- Fixtures ------------------------------------------------------------
function st = one_state()
st = struct('mech_az', 0, 'tilt', 0, 'steer_az', 0, 'steer_el', 0);
end

function st = states(n)
z = zeros(n, 1);
st = struct('mech_az', z, 'tilt', z, 'steer_az', z, 'steer_el', z);
end

function st = two_states()
st = struct('mech_az', [0; 0], 'tilt', [0; 0], 'steer_az', [0; 0], 'steer_el', [0; 0]);
end

function a = csc2_antenna()
a = struct('type', 'rect', 'g_max_dBi', 33.5, 'az_bw_deg', 1.35, 'az_sll_dB', -26.2, ...
    'el_bw_deg', 4.8, 'el_pattern', 'csc2', 'csc2_end_deg', 30, ...
    'csc2_floor_dB', -55, 'envelope', 'peak');
end

function a = array_antenna(n, taper)
lam = physconst('LightSpeed')/3e9;
a = struct('type', 'array_planar', 'nx', n, 'ny', n, 'dx_m', lam/2, 'dy_m', lam/2, ...
    'freq_MHz', 3000, 'elem_gain_dBi', 5, 'elem_cos_exp', 2, 'taper', taper, ...
    'taylor_sll_dB', -25, 'floor_dB', -60, 'front_back_dB', -60);
end

function [s1, s2] = two_radars(rate1, rate2)
a = csc2_antenna();
b = a;
b.g_max_dBi = 34;
b.az_bw_deg = 1.3;
b.az_sll_dB = -24.5;
b.el_bw_deg = 4;
s1 = struct('name', 'A', 'tx', a, ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', rate1, 'tilt_deg', 2));
s2 = struct('name', 'B', 'tx', b, ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', rate2, 'tilt_deg', 2));
end

function s = two_array_radars()
[s1, ~] = two_radars(75, 90);
arr = array_antenna(24, "taylor");
s2 = struct('name', 'AESA', 'tx', arr, ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', 60, 'beam_mode', 'random', ...
    'beam_az_lim_deg', [-10 10], 'beam_el_lim_deg', [0 20], 'dwell_s', 0.05));
s = {s1, s2};
end

function [s1, s2] = pencil_vs_flat()
% 40 dBi, 2 deg pencil beam rotating at 36 deg/s against a flat 0 dBi
% antenna, so the coupling equals the pencil beam's gain.
pencil = struct('type', 'parabolic', 'g_max_dBi', 40, 'az_bw_deg', 2, ...
    'az_sll_dB', -16, 'envelope', 'peak');
flat = struct('type', 'measured', 'g_max_dBi', 0, 'az_deg', [-180 180], ...
    'az_gain_dB', [0 0], 'floor_dB', -100, 'az_bw_deg', 360);
s1 = struct('name', 'pencil', 'tx', pencil, ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', 36, 'tilt_deg', 0));
s2 = struct('name', 'flat', 'tx', flat, ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', 7.31));
end

function out = evalc_quiet(fcn)
% Run fcn with the commensurate-rate warning silenced.
w = warning('off', 'radar_coupling:commensurateRates');
restore = onCleanup(@() warning(w));
out = fcn();
end

function o = quick_opts(seed)
o = struct('seed', seed, 'mode', 'random', 'n_samples', 2e5, 'chunk_size', 5e4);
end

function o = sync_opts(seed, offset_deg)
o = quick_opts(seed);
o.scan_sync = "synchronized";
o.sync_offset_deg = offset_deg;
end

function [s1, s2] = fg_radars(g1, g2, rate)
% Two F.1245 radars at one rotation rate, tilt 0.
[s1, s2] = sg_radars(g1, g2, rate);
s1.tx.type = 'f1245';
s2.tx.type = 'f1245';
end

function [s1, s2] = sg_radars(g1, g2, rate)
% Two Statgain radars at one rotation rate, tilt 0.
s1 = struct('name', 'A', 'tx', struct('g_max_dBi', g1), ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', rate));
s2 = struct('name', 'B', 'tx', struct('g_max_dBi', g2), ...
    'scan', struct('type', 'mech_rotating', 'rate_deg_s', rate));
end
