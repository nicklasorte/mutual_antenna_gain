function tests = test_validation_rev1
%TEST_VALIDATION_REV1 Input checks and less common branches of the radar
%   coupling code: every error identifier that the other test files do not
%   reach, measured elevation cuts, an explicit CSC^2 start angle, the
%   "theoretical" sidelobe envelope, Statgain back levels used as array
%   floors, fixed and sequential array beams, and synchronized arrays.
%   Found with MATLAB code coverage over all test files.
%   Run with runtests("test_validation_rev1").
tests = functiontests(localfunctions);
end

% =========================================================================
% antenna_prepare_rev1
% =========================================================================
function testAntennaPrepareErrors(tc)
verifyError(tc, @() antenna_prepare_rev1(struct("type", "dish", "g_max_dBi", 30)), "radar_coupling:badAntennaType");
r = rect();
r.el_pattern = "fan";
verifyError(tc, @() antenna_prepare_rev1(r), "radar_coupling:badElPattern");
r = rect();
r.envelope = "smooth";
verifyError(tc, @() antenna_prepare_rev1(r), "radar_coupling:badEnvelope");
m = measured();
m.el_deg = [-90 0 90];
verifyError(tc, @() antenna_prepare_rev1(m), "radar_coupling:missingField");
m = measured();
m.az_gain_dB = [0 0];
verifyError(tc, @() antenna_prepare_rev1(m), "radar_coupling:badMeasuredCut");
m = measured();
m.az_deg = [-180 10 0];
verifyError(tc, @() antenna_prepare_rev1(m), "radar_coupling:badMeasuredCut");
verifyError(tc, @() antenna_prepare_rev1(struct("type", "f1245", "g_max_dBi", 34, "d_over_lambda", -1)), ...
    "radar_coupling:f1245Range");
verifyError(tc, @() antenna_prepare_rev1(struct("type", "array_planar", "eff_dB", 1)), "radar_coupling:badEfficiency");
verifyError(tc, @() antenna_prepare_rev1(struct("type", "array_planar", "taper", "hann")), "radar_coupling:badTaper");
% One element of 3 dBi: below the 10 dBi Statgain range, so no floor default.
verifyError(tc, @() antenna_prepare_rev1(struct("type", "array_planar", "nx", 1, "ny", 1, "elem_gain_dBi", 3)), ...
    "radar_coupling:missingField");
end

function testArrayFloorFromStatgainBackLevel(tc)
% NTIA TM-13-489 Table 4-1 back level, used as the array floor default:
% -13 dBi for G >= 48 dBi, 0 dBi for G < 22 dBi, 11 - G/2 between.
for G = [15 30 50]
    p = antenna_prepare_rev1(struct("type", "array_planar", "g_max_dBi", G));
    back = (G >= 48)*(-13) + (G >= 22 && G < 48)*(11 - p.g_max_dBi/2);
    verifyEqual(tc, p.floor_dB + p.g_max_dBi, back, "AbsTol", 1e-9, sprintf("G = %g", G));
end
end

% =========================================================================
% antenna_gain_rev1
% =========================================================================
function testMeasuredElevationCut(tc)
% Gain = peak + azimuth cut + elevation cut (relative to the tilt).
m = measured();
m.el_deg = [-90 0 10 90];
m.el_gain_dB = [-40 0 -6 -40];
p = antenna_prepare_rev1(m);
st = states(3, 0);
st.tilt = [0; 0; 2];
g = antenna_gain_rev1(p, 0, [0; 10; 12], st);
verifyEqual(tc, g, 20 + [0; -6; -6], "AbsTol", 1e-9);
end

function testUnpreparedTypeError(tc)
verifyError(tc, @() antenna_gain_rev1(struct("type", "dish"), 0, 0, states(1, 0)), "radar_coupling:badAntennaType");
end

function testTheoreticalEnvelopeFloorBeyond90(tc)
% "theoretical" keeps the M.1851-2 pattern itself and drops to the floor
% beyond 90 deg.
r = rect();
r.envelope = "theoretical";
p = antenna_prepare_rev1(r);
st = states(2, 0);
st.mech_az = [120; 180];
g = antenna_gain_rev1(p, 0, 0, st);
verifyEqual(tc, g, p.g_max_dBi + p.az_cut.floor_dB + [0; 0], "AbsTol", 1e-9);
end

function testParabolicTheoreticalEnvelope(tc)
% Parabolic "theoretical" envelope: Eq. (34) pattern, floor beyond 90 deg.
p = antenna_prepare_rev1(struct("type", "parabolic", "g_max_dBi", 40, "az_bw_deg", 2, ...
    "az_sll_dB", -16, "envelope", "theoretical"));
st = states(2, 0);
st.mech_az = [120; 0];
g = antenna_gain_rev1(p, 0, 0, st);
verifyEqual(tc, g(1), p.g_max_dBi + p.floor_dB, "AbsTol", 1e-9);
verifyEqual(tc, g(2), p.g_max_dBi, "AbsTol", 1e-9);
end

function testExplicitCsc2Start(tc)
% csc2_start_deg fixes where the cosecant-squared shape starts:
% G(th) = G(start) + 20 log10(sin(start)/sin(th)) above it.
r = rect();
r.el_pattern = "csc2";
r.csc2_end_deg = 30;
r.csc2_floor_dB = -30;
r.csc2_start_deg = 3;
p = antenna_prepare_rev1(r);
st = states(2, 0);
g = antenna_gain_rev1(p, 0, [3; 10], st);
verifyEqual(tc, g(2) - g(1), 20*log10(sind(3)/sind(10)), "AbsTol", 1e-9);
r.csc2_start_deg = -1;
verifyError(tc, @() antenna_gain_rev1(antenna_prepare_rev1(r), 0, 5, st), "radar_coupling:csc2Start");
end

% =========================================================================
% radar_mutual_gain_cdf_rev1
% =========================================================================
function testSystemAndScanErrors(tc)
s = flat_sys();
o = struct("seed", 1, "n_samples", 1e3);
verifyError(tc, @() radar_mutual_gain_cdf_rev1(rmfield(s, "tx"), s, o), "radar_coupling:missingField");
a = s;
a.tx = struct("type", "array_planar");
a.rx = s.tx;                                % array tx with an aperture rx
verifyError(tc, @() radar_mutual_gain_cdf_rev1(a, s, o), "radar_coupling:txRxMismatch");
b = s;
b.scan = rmfield(b.scan, "type");
verifyError(tc, @() radar_mutual_gain_cdf_rev1(b, s, o), "radar_coupling:missingField");
b = s;
b.scan = rmfield(b.scan, "rate_deg_s");
verifyError(tc, @() radar_mutual_gain_cdf_rev1(b, s, o), "radar_coupling:missingField");
b = s;
b.scan.rate_deg_s = 0;
verifyError(tc, @() radar_mutual_gain_cdf_rev1(b, s, o), "radar_coupling:zeroRate");
b = s;
b.scan.rate_deg_s = [10 20];
b.scan.n_rev = 1.5;
verifyError(tc, @() radar_mutual_gain_cdf_rev1(b, s, o), "radar_coupling:badNRev");
b.scan.n_rev = 1;
b.scan.tilt_deg = [0 1 2];
verifyError(tc, @() radar_mutual_gain_cdf_rev1(b, s, o), "radar_coupling:stepMismatch");
r = arr_sys();
r.scan.beam_mode = "sweep";
verifyError(tc, @() radar_mutual_gain_cdf_rev1(r, s, o), "radar_coupling:badBeamMode");
f = face_sys();
f.scan.beam_mode = "fixed";
verifyError(tc, @() radar_mutual_gain_cdf_rev1(f, s, o), "radar_coupling:badBeamMode");
end

function testOptionErrors(tc)
s = flat_sys();
verifyError(tc, @() radar_mutual_gain_cdf_rev1(s, s, struct("seed", 1, "parallel_flag", 2)), ...
    "radar_coupling:badParallelFlag");
o = struct("seed", 1, "scan_sync", "synchronized", "sync_offset_deg", NaN);
verifyError(tc, @() radar_mutual_gain_cdf_rev1(s, s, o), "radar_coupling:badSyncOffset");
o.sync_offset_deg = 0;
c = s;
c.scan.sector_deg = 90;
verifyError(tc, @() radar_mutual_gain_cdf_rev1(c, s, o), "radar_coupling:syncSector");
a = arr_sys();
a.scan.beam_mode = "random";
verifyError(tc, @() radar_mutual_gain_cdf_rev1(a, arr_sys(), o), "radar_coupling:syncRandomBeam");
verifyError(tc, @() radar_mutual_gain_cdf_rev1(s, s, struct("seed", 1, "mode", "time_series", ...
    "duration_s", 1e7, "n_trials", 100)), "radar_coupling:tooManySteps");
end

function testSynchronizedFixedBeamArrays(tc)
% Two rotating arrays with fixed beams (broadside), same rate, offset 180:
% they face each other once a turn, so the maximum is the sum of the two
% broadside gains.
a = arr_sys();
o = struct("seed", 1, "scan_sync", "synchronized", "sync_offset_deg", 180);
out = radar_mutual_gain_cdf_rev1(a, a, o);
p = antenna_prepare_rev1(a.tx);
verifyEqual(tc, out.sync.max_g12_dBi, 2*p.g_max_dBi, "AbsTol", 0.05);
end

function testFixedAndSequentialBeamsTimeSeries(tc)
% Fixed-beam rotating array and sequential multi-face array in the time
% series: both run, and every MAG stays at or below the sum of the peaks.
s = flat_sys();
o = struct("seed", 2, "mode", "time_series", "duration_s", 20, "n_trials", 1);
a = arr_sys();
r1 = radar_mutual_gain_cdf_rev1(a, s, o);
f = face_sys();
f.scan.beam_mode = "sequential";
f.scan.rate_deg_s = 7.3;                    % 0.365 deg per 0.05 s dwell
o.duration_s = 60;                          % one full turn of the beam
r2 = radar_mutual_gain_cdf_rev1(f, s, o);
pa = antenna_prepare_rev1(a.tx);
pf = antenna_prepare_rev1(f.tx);
verifyLessThanOrEqual(tc, max(r1.g12_dBi), pa.g_max_dBi + 1e-9);
verifyLessThanOrEqual(tc, max(r2.g12_dBi), pf.g_max_dBi + 1e-9);
% A sequential beam steps around the faces, so in one turn it points at
% the flat radar within 0.4 deg, from a face within 45 deg of broadside.
verifyGreaterThan(tc, max(r2.g12_dBi), pf.g_max_dBi - 6);
end

function testRandomMultifaceTimeSeries(tc)
% Random beams on a multi-face array: every dwell draws a new azimuth.
o = struct("seed", 5, "mode", "time_series", "duration_s", 10, "n_trials", 1);
r = radar_mutual_gain_cdf_rev1(face_sys(), flat_sys(), o);
pf = antenna_prepare_rev1(face_sys().tx);
verifyLessThanOrEqual(tc, max(r.g12_dBi), pf.g_max_dBi + 1e-9);
verifyGreaterThan(tc, numel(unique(round(r.g12_dBi, 6))), 20, "new beam azimuth each dwell (many dwells sit on the floor)");
end

function testMeasuredWithoutBeamwidthWarns(tc)
m = measured();
s1 = struct("name", "m", "tx", m, "scan", struct("type", "mech_rotating", "rate_deg_s", 30));
o = struct("seed", 3, "mode", "time_series", "duration_s", 5, "n_trials", 1);
verifyWarning(tc, @() radar_mutual_gain_cdf_rev1(s1, flat_sys(), o), "radar_coupling:measuredBw");
end

% =========================================================================
% Other functions
% =========================================================================
function testOtherErrors(tc)
verifyError(tc, @() alignment_dwell_rev1(1000, 1, 1, 1, max_ratio=2, jitter_deg_s=[0 0]), "radar_coupling:noRatio");
verifyError(tc, @() mag_mode_envelope_rev1({1:10}, "A", 100), "radar_coupling:envelopeInput");
verifyError(tc, @() mag_mode_envelope_rev1({1:10, 1:5}, "A", 99), "radar_coupling:envelopeInput");
verifyError(tc, @() nexrad_asr11_mag_rev1(struct("seed", 1, "duration_s", 1e6, "n_trials", 1000)), ...
    "radar_coupling:tooManySteps");
end

function testSensitivityAddsThreshold45(tc)
% A base grid without 45 dBi gets it added, so max_run_45_revs exists.
b = struct("seed", 3, "vcp", 35, "n_samples", 1e4, "duration_s", 60, "n_trials", 1, "persist_thr_dBi", [50 60]);
c = struct("group", "ASR-11 tilt", "label", "+1 deg", "opts", struct("asr_tilt_deg", 1));
T = mag_sensitivity_rev1(b, c);
verifyEqual(tc, height(T), 2);
verifyTrue(tc, all(isfinite(T.max_run_45_revs)));
end

% =========================================================================
function a = rect()
a = struct("type", "rect", "g_max_dBi", 30, "az_bw_deg", 2, "az_sll_dB", -25, ...
    "el_bw_deg", 5, "el_sll_dB", -20, "envelope", "peak");
end

function m = measured()
m = struct("type", "measured", "g_max_dBi", 20, "az_deg", [-180 0 180], ...
    "az_gain_dB", [-30 0 -30], "floor_dB", -40);
end

function st = states(n, tilt)
z = zeros(n, 1);
st = struct("mech_az", z, "tilt", tilt + z, "steer_az", z, "steer_el", z);
end

function s = flat_sys()
flat = struct("type", "measured", "g_max_dBi", 0, "az_deg", [-180 180], ...
    "az_gain_dB", [0 0], "floor_dB", -100, "az_bw_deg", 360);
s = struct("name", "flat", "tx", flat, "scan", struct("type", "mech_rotating", "rate_deg_s", 7.31));
end

function s = arr_sys()
lam = physconst("LightSpeed")/3e9;
arr = struct("type", "array_planar", "nx", 16, "ny", 16, "dx_m", lam/2, "dy_m", lam/2, ...
    "freq_MHz", 3000, "elem_gain_dBi", 5, "elem_cos_exp", 2, "taper", "taylor", ...
    "taylor_sll_dB", -25, "floor_dB", -60, "front_back_dB", -60);
s = struct("name", "array", "tx", arr, "scan", struct("type", "mech_rotating", "rate_deg_s", 60));
end

function s = face_sys()
s = arr_sys();
s.scan = struct("type", "array_multiface", "dwell_s", 0.05);
end
