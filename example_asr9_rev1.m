%EXAMPLE_ASR9_REV1 WSR-88D vs ASR-9 mutual antenna gain, compared with the ASR-11.
%   Runs nexrad_asr11_mag_rev1 with atc_radar = "asr9" for every WSR-88D
%   VCP mode (vcp_mode_cases_rev1, 74 modes), finds the envelope over the
%   modes, and compares the ASR-9 with the ASR-11 at the ASR-9's worst mode.
%   ASR-9 values come from Rec. ITU-R M.1464-2 Table 1, Radar B (33.5 dBi,
%   1.3 deg azimuth, 4.8 deg elevation, +7.3 dBi first sidelobe, 75 deg/s);
%   see nexrad_asr11_mag_rev1 (assumptions, citations). M.1464-2 gives no
%   ASR-9 remote sidelobe level, so the ASR-11 value (+3.5 dBi) is used
%   (TODO). The paper baseline (F.1245, 34 dBi) is the same for both radars
%   and isn't rerun.
%
%   Set asr9_base before running to change the run size; the default
%   matches example_nexrad_asr11_rev1 (seed 20261105, 5e6 samples, four 1 h
%   trials per mode). The comparison runs asr9_detail_trials 1 h trials
%   (default 24) for both radars with the same seed, so the two runs share
%   their random pointing and rotation-rate draws and differ only in the
%   ATC radar antenna. This script doesn't clear the workspace, so
%   report_figures_rev5 can call it.
%
%   Defines:
%     T9       one row per VCP mode: p99 to p99_999 (time-weighted MAG),
%              rev_pct_ge_40/50/60 (% of ASR-9 revolutions whose maximum
%              reaches the level), vol_max (largest per-volume maximum)
%     env9     mag_mode_envelope_rev1 output over the 74 modes
%     w9, w11  ASR-9 and ASR-11 runs at the comparison mode (cmp9_label)
%     Tcmp9    ASR-11 and ASR-9 side by side at the comparison mode, with
%              the difference (ASR-9 minus ASR-11)
%
%   See also nexrad_asr11_mag_rev1, example_nexrad_asr11_rev1,
%   vcp_mode_cases_rev1, report_figures_rev5.

if ~exist("asr9_base", "var")
    asr9_base = struct("seed", 20261105, "n_samples", 5e6, "duration_s", 3600, "n_trials", 4);
end
asr9_detail_trials = 24;
if exist("example_overrides", "var")
    asr9_base = merge_overrides_rev1(asr9_base, example_overrides);
    if isfield(example_overrides, "n_trials")
        asr9_detail_trials = example_overrides.n_trials;
    end
end
asr9_base.atc_radar = "asr9";

% ---- Sweep the 74 VCP modes for the ASR-9 (c9{1} is the paper baseline)
c9 = vcp_mode_cases_rev1();
c9 = c9(2:end);
env_pct9 = [99 99.9 99.99 99.999];
n9 = numel(c9);
summ9 = cell(1, n9);
rows9 = cell(n9, 1);
for i9 = 1:n9
    o9 = asr9_base;
    f9 = setdiff(fieldnames(c9{i9}), "label");
    for k9 = 1:numel(f9)
        o9.(f9{k9}) = c9{i9}.(f9{k9});
    end
    r9 = nexrad_asr11_mag_rev1(o9);
    summ9{i9} = mag_mode_summary_rev1(r9.cdf.g21_dBi, env_pct9);
    q9 = r9.cdf.pct21;
    R9 = r9.scan_max.asr_rev(r9.scan_max.asr_rev.complete, :);
    rows9{i9} = {string(c9{i9}.label), q9.gain_dBi(q9.pct == 99), q9.gain_dBi(q9.pct == 99.9), ...
        q9.gain_dBi(q9.pct == 99.99), q9.gain_dBi(q9.pct == 99.999), ...
        100*mean(R9.max_g21_dBi >= 40), 100*mean(R9.max_g21_dBi >= 50), ...
        100*mean(R9.max_g21_dBi >= 60), max(r9.scan_max.nexrad_volume.max_g21_dBi)};
    fprintf("ASR-9 %-34s done\n", c9{i9}.label);
end
clear r9
T9 = cell2table(vertcat(rows9{:}), "VariableNames", ["case", "p99", "p99_9", "p99_99", ...
    "p99_999", "rev_pct_ge_40", "rev_pct_ge_50", "rev_pct_ge_60", "vol_max"]);
env9 = mag_mode_envelope_rev1(summ9, T9.case', env_pct9);
[~, j9] = max(T9.p99_9);
cmp9_label = T9.case(j9);
fprintf("\nASR-9 modes ranked by 99.9%% MAG (worst first):\n");
disp(sortrows(T9, ["p99_9", "p99"], "descend"))
fprintf("ASR-9 mode envelope (bound for any mix of the %d VCP modes):\n", n9);
disp(env9.pct(:, ["pct", "envelope_dBi", "worst_mode"]))

% ---- Paired comparison at the ASR-9 worst mode: same seed and draws,
%      only the ATC radar antenna differs
o9 = asr9_base;
f9 = setdiff(fieldnames(c9{j9}), "label");
for k9 = 1:numel(f9)
    o9.(f9{k9}) = c9{j9}.(f9{k9});
end
o9.n_trials = asr9_detail_trials;
w9 = nexrad_asr11_mag_rev1(o9);
o9.atc_radar = "asr11";
w11 = nexrad_asr11_mag_rev1(o9);
w9.cdf.g12_dBi = [];
w11.cdf.g12_dBi = [];

% ---- Side-by-side table
pick9 = @(w, p) w.cdf.pct21.gain_dBi(abs(w.cdf.pct21.pct - p) < 1e-9);
rev9 = @(w, v) 100*mean(w.scan_max.asr_rev.max_g21_dBi(w.scan_max.asr_rev.complete) >= v - 0.2);
vol9 = @(w, v) w.scan_max.persistence_volume.max_run_vols(w.scan_max.persistence_volume.threshold_dBi == v);
mbmb9 = @(w) w.lobe_pairs.time_pct(w.lobe_pairs.nexrad_region == "main beam" & w.lobe_pairs.asr11_region == "main beam");
qty9 = ["Peak gains, main beam to main beam (dBi)"
    "Main beam to main beam on the path (dBi)"
    "NEXRAD main beam + ATC near sidelobe (dBi)"
    "NEXRAD main beam + ATC back lobe (dBi)"
    "ATC main beam + NEXRAD near sidelobe (dBi)"
    "ATC main beam + NEXRAD back lobe (dBi)"
    "99% time-weighted MAG (dBi)"
    "99.9% time-weighted MAG (dBi)"
    "99.99% time-weighted MAG (dBi)"
    "99.999% time-weighted MAG (dBi)"
    "Revolutions reaching 45 dBi (%)"
    "Revolutions reaching 60 dBi (%)"
    "Longest run of volumes reaching 60 dBi"
    "Time with both main beams on the path (%)"
    "ATC main-beam region edge (deg)"];
val9 = @(w) [w.sensitivity.mag_dBi(1); w.ref_levels.mag_dBi; arrayfun(@(p) pick9(w, p), env_pct9(:)); ...
    rev9(w, 45); rev9(w, 60); vol9(w, 60); mbmb9(w); w.lobe_pairs.Properties.UserData.asr11(1)];
v11 = val9(w11);
v9 = val9(w9);
Tcmp9 = table(qty9, v11, v9, v9 - v11, 'VariableNames', ["quantity", "asr11", "asr9", "asr9_minus_asr11"]);
Tcmp9.Properties.Description = sprintf("%s, seed %d, %d trials of %g s, %g samples", cmp9_label, ...
    o9.seed, o9.n_trials, o9.duration_s, o9.n_samples);
fprintf("\nASR-11 and ASR-9 at %s (same seed and draws):\n", cmp9_label);
disp(Tcmp9)
disp(w9.assumptions(startsWith(w9.assumptions.parameter, "ASR-9"), :))
