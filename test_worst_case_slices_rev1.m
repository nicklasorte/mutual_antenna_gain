function tests = test_worst_case_slices_rev1
%TEST_WORST_CASE_SLICES_REV1 Tests for the worst-case horizontal slice functions.
%   Checks nexrad_, asr11_, asr9_ and statgain_worst_case_slice_rev1:
%   output shape, hand-worked values, and agreement with the repository
%   antenna model (antenna_prepare_rev1 and antenna_gain_rev1) at the same
%   settings.
%
%   DRAFT//PRE-DECISIONAL
tests = functiontests(localfunctions);
end

function testShapeAndAzimuthGrid(tc)
fns = {@nexrad_worst_case_slice_rev1, @asr11_worst_case_slice_rev1, @asr9_worst_case_slice_rev1, ...
    @() statgain_worst_case_slice_rev1(45.5)};
for k = 1:numel(fns)
    s = fns{k}();
    verifySize(tc, s, [3601 2]);
    verifyEqual(tc, s([1 end], 1), [0; 360]);
    verifyEqual(tc, s(1, 2), s(end, 2), "AbsTol", 1e-12, "0 and 360 deg are the same pointing");
    verifyEqual(tc, s(2:1800, 2), flip(s(1802:3600, 2)), "AbsTol", 1e-9, "symmetric about 180 deg");
end
s = nexrad_worst_case_slice_rev1(step_deg=0.25);
verifySize(tc, s, [1441 2]);
verifyError(tc, @() asr11_worst_case_slice_rev1(step_deg=0.7), "asr_worst_case_slice:badStep");
verifyError(tc, @() nexrad_worst_case_slice_rev1(step_deg=0.7), "nexrad_worst_case_slice:badStep");
end

function testNexradHandValues(tc)
[s, info] = nexrad_worst_case_slice_rev1(step_deg=0.005);
g = @(a) s(abs(s(:, 1) - a) < 1e-9, 2);
verifyEqual(tc, g(0), 45.5, "AbsTol", 1e-12, "peak on the horizon");
s3 = nexrad_worst_case_slice_rev1(step_deg=0.0125);
verifyEqual(tc, s3(38, 2), 42.5, "AbsTol", 1e-9, "3 dB down at half the 0.925 deg beamwidth (0.4625 deg)");
verifyEqual(tc, info.x_deg, 0.925*sqrt(29/12), "AbsTol", 1e-12);
verifyEqual(tc, g(5), 16.5, "AbsTol", 1e-12, "-29 dB to 10 deg");
verifyEqual(tc, g(10), 16.5, "AbsTol", 1e-12);
verifyEqual(tc, g(10.005), 5.5, "AbsTol", 1e-12, "-40 dB beyond 10 deg");
verifyEqual(tc, g(180), 5.5, "AbsTol", 1e-12);
verifyEqual(tc, info.far_field_m, 1384, "AbsTol", 2);
st = nexrad_worst_case_slice_rev1(tilt_deg=0.5);
verifyEqual(tc, st(1, 2), 45.5 - 12*(0.5/0.925)^2, "AbsTol", 1e-9, "0.5 deg cut puts 42.0 dBi on the horizon");
verifyLessThanOrEqual(tc, st(:, 2), s(1:20:end, 2) + 1e-12, "tilt 0 bounds the 0.5 deg tilt");
[sn, in] = nexrad_worst_case_slice_rev1(pattern="nssl_measured", step_deg=0.005);
gn = @(a) sn(abs(sn(:, 1) - a) < 1e-9, 2);
verifyEqual(tc, gn(2), 10.5, "AbsTol", 1e-9);
verifyEqual(tc, gn(6), -2.5, "AbsTol", 1e-9);
verifyEqual(tc, gn(20), -9.5, "AbsTol", 1e-9);
verifyEqual(tc, gn(4), 45.5 - 41.5, "AbsTol", 1e-9, "linear between 2 and 6 deg");
verifyEqual(tc, in.floor_dBi, -9.5);
verifyLessThanOrEqual(tc, sn(:, 2), s(:, 2) + 1e-12, "the ROC specification bounds the measurement");
end

function testNexradMatchesRepoModel(tc)
for tilt = [0 0.5]
    for pat = ["roc_spec", "nssl_measured"]
        [s, info] = nexrad_worst_case_slice_rev1(pattern=pat, tilt_deg=tilt, step_deg=0.05);
        % Repository pencil_table (nexrad_asr11_mag_rev1), rebuilt here.
        if pat == "roc_spec"
            env_deg = [0 10 10.001 180]; env_dB = [-29 -29 -40 -40];
        else
            env_deg = [0 2 6 20 90 90.001 180]; env_dB = [-35 -35 -48 -55 -55 -55 -55];
        end
        th = unique([(0:0.005:5)'; (5:0.05:180)'; env_deg(:)]);
        env = interp1(env_deg, env_dB, th, "linear");
        main = -12*(th/0.925).^2;
        x = find(main < env, 1);
        gt = main; gt(x:end) = env(x:end);
        ant_p = antenna_prepare_rev1(struct("type", "offaxis_table", "g_max_dBi", 45.5, ...
            "offaxis_deg", th, "gain_dB", gt, "az_bw_deg", 0.925));
        st = struct("mech_az", s(:, 1), "tilt", tilt*ones(height(s), 1));
        g_repo = antenna_gain_rev1(ant_p, 0, 0, st);
        theta = acosd(cosd(tilt)*cosd(s(:, 1)));
        keep = abs(theta - info.x_deg) > 0.006;     % the table grid smooths the corner
        verifyEqual(tc, s(keep, 2), g_repo(keep), "AbsTol", 1e-3, pat + " tilt " + tilt);
        verifyEqual(tc, s(~keep, 2), g_repo(~keep), "AbsTol", 0.25, pat + " tilt " + tilt + " corner");
    end
end
end

function testAsrHandValues(tc)
[s, info] = asr11_worst_case_slice_rev1(step_deg=0.01);
g = @(a) s(abs(s(:, 1) - a) < 1e-9, 2);
verifyEqual(tc, g(0), 34, "AbsTol", 1e-12, "beam nose on the horizon");
verifyEqual(tc, info.g_el_dB, 0, "AbsTol", 1e-12);
verifyEqual(tc, info.n, 1, "cosine taper for the -24.5 dB first sidelobe");
c = fzero(@(a) 20*log10(abs(cos_field(68.8, 1.45, a))) + 3, [0.1 1.4]);
verifyEqual(tc, 2*c, 1.4334, "AbsTol", 1e-4, "M.1851-2 K = 68.8 gives 1.433 deg, 1.1% under the 1.45 deg input");
verifyEqual(tc, g(2), 34 - 17.51*log(2.33*2/1.45), "AbsTol", 1e-9, "peak envelope");
verifyEqual(tc, g(180), 3.5, "AbsTol", 1e-12, "remote sidelobes");
verifyEqual(tc, info.floor_az_deg, 1.45/2.33*exp(30.5/17.51), "AbsTol", 1e-12);
verifyEqual(tc, 20*log10(abs(cos_field(68.8, 1.45, info.brk_deg))), -14.4, "AbsTol", 1e-9);
[s9, i9] = asr9_worst_case_slice_rev1(step_deg=0.01);
g9 = @(a) s9(abs(s9(:, 1) - a) < 1e-9, 2);
verifyEqual(tc, g9(0), 33.5, "AbsTol", 1e-12);
verifyEqual(tc, g9(90), 5.5, "AbsTol", 1e-12, "WEBER_1986 -28 dB bound");
verifyEqual(tc, g9(2), 33.5 - 17.51*log(2.33*2/1.3), "AbsTol", 1e-9);
verifyEqual(tc, i9.floor_az_deg, 1.3/2.33*exp(28/17.51), "AbsTol", 1e-12);
st = asr11_worst_case_slice_rev1(tilt_deg=2.5, step_deg=0.01);
verifyLessThan(tc, st(1, 2), 34, "tilting up lowers the horizon gain");
verifyLessThanOrEqual(tc, st(:, 2), s(:, 2) + 1e-12, "tilt 0 bounds every tilt");
verifyError(tc, @() asr11_worst_case_slice_rev1(tilt_deg=-2.5), "asr_worst_case_slice:badTilt");
end

function testAsrMatchesRepoModel(tc)
cases = struct("fn", {@asr11_worst_case_slice_rev1, @asr9_worst_case_slice_rev1}, ...
    "g", {34, 33.5}, "bw", {1.45, 1.3}, "sl", {9.5, 7.3}, "rsl", {3.5, 5.5});
for c = cases
    for tilt = [0 1 2.5]
        s = c.fn(tilt_deg=tilt, step_deg=0.05);
        fl = c.rsl - c.g;
        ant_p = antenna_prepare_rev1(struct("type", "rect", "g_max_dBi", c.g, ...
            "az_bw_deg", c.bw, "az_sll_dB", c.sl - c.g, "el_bw_deg", 4.8, "el_pattern", "csc2", ...
            "csc2_end_deg", 30, "csc2_floor_dB", fl, "envelope", "peak", "floor_dB", fl));
        st = struct("mech_az", s(:, 1), "tilt", tilt*ones(height(s), 1));
        g_repo = antenna_gain_rev1(ant_p, 0, 0, st);
        verifyEqual(tc, s(:, 2), g_repo, "AbsTol", 1e-9, func2str(c.fn) + " tilt " + tilt);
    end
end
end

function testAsrElevationPattern(tc)
[~, info] = asr11_worst_case_slice_rev1();
g_uni = @(x) 20*log10(abs(sinc(50.8*sind(x)/4.8)));
verifyEqual(tc, info.el_fn(0), 0, "AbsTol", 1e-12, "elevation peak at the tilt (0 deg)");
verifyEqual(tc, info.el_fn(2.4), g_uni(2.4), "AbsTol", 1e-12, "skirt starts at tilt + 2.4 deg");
verifyEqual(tc, info.el_fn(2.4), -3.0, "AbsTol", 0.05, "about 3 dB down there");
verifyEqual(tc, info.el_fn(10), 20*log10(sind(2.4)/sind(10)) + g_uni(2.4), "AbsTol", 1e-9, "csc^2 law");
verifyEqual(tc, info.el_fn(-1), g_uni(-1), "AbsTol", 1e-12, "uniform lower edge");
verifyEqual(tc, info.el_fn(40), -30.5, "AbsTol", 1e-12, "floor above 30 deg");
end

function testStatgainMatchesRepoModel(tc)
for G = [12 33.5 34 45.5 50]
    [s, info] = statgain_worst_case_slice_rev1(G, step_deg=0.05);
    ant_p = antenna_prepare_rev1(struct("type", "statgain", "g_max_dBi", G));
    st = struct("mech_az", s(:, 1), "tilt", zeros(height(s), 1));
    verifyEqual(tc, s(:, 2), antenna_gain_rev1(ant_p, 0, 0, st), "AbsTol", 1e-9, "G = " + G);
    verifyEqual(tc, s(1, 2), G, "AbsTol", 1e-12);
    verifyEqual(tc, info.peak_dBi, G, "AbsTol", 1e-12);
end
[s, info] = statgain_worst_case_slice_rev1(45.5);
verifyEqual(tc, info.category, 2);
verifyEqual(tc, info.back_dBi, 11 - 45.5/2, "AbsTol", 1e-12);
verifyEqual(tc, s(end - 1, 2), s(2, 2), "AbsTol", 1e-12);
verifyEqual(tc, info.bw_implied_deg, 0.919, "AbsTol", 1e-3, "matches the WSR-88D beamwidth");
[~, i34] = statgain_worst_case_slice_rev1(34);
verifyGreaterThan(tc, i34.bw_implied_deg, 3, "much wider than the ASR-11 1.45 deg");
verifyError(tc, @() statgain_worst_case_slice_rev1(9), "MATLAB:validators:mustBeGreaterThanOrEqual");
end

function F = cos_field(K, bw, a)
mu = pi*K*sind(a)/bw;
F = (sinc((mu + pi/2)/pi) + sinc((mu - pi/2)/pi))/(2*sinc(0.5));
end
