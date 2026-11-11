%REPORT_NUMBERS_REV6 Headline tables for methodology report rev6.
%   Runs the revision 6 defaults (ROC specification WSR-88D pattern, all
%   216 scan modes, ASR-9 remote sidelobe +5.5 dBi, exact CDFs, event time
%   series) and writes to docs/figures:
%     headline_rev6.csv     exact MAG percentiles at a horizontal path (0/0
%                           deg) for each case: envelope over the 216 modes
%                           and over the 74 Rev 5 modes, the worst mode,
%                           and the largest MAG on the path
%     worst_modes_rev6.csv  mode that sets each envelope percentile, each case
%     noaa_tilt_rev6.csv    NOAA case: 99.99% envelope against ATC tilt
%     timeseries_rev6.csv   120 h event time series at the worst mode, each
%                           case: share of ATC revolutions and NEXRAD
%                           volumes reaching each level, number of runs,
%                           longest run
%     harm_rev6.csv         harm metric: share of 0.5 deg cuts with at
%                           least one radial at or above each level, mean
%                           and largest number of such radials per cut, for
%                           1 deg and 0.5 deg radials
%     lookup_pair_rev6.csv  two-part model, worked pair (2.0/1.7 deg),
%                           envelope over the 216 modes, both ATC radars
%     ref_levels_rev6.csv   main-beam and sidelobe reference levels at the
%                           worst mode (ROC specification, ASR-11)
%   Cases: "roc_asr11" (baseline), "roc_asr9", "nssl_asr11" and
%   "nssl_asr9" (Rev 5 pattern, sensitivity), "noaa_asr11" (ROC
%   specification, ASR remote sidelobe +6.5 dBi, worst ATC tilt). The ATC
%   remote sidelobe is each radar's default (ASR-11 +3.5 dBi, ASR-9
%   +5.5 dBi) except in the NOAA case.
%
%   Set rev6_overrides (struct: out_dir, n_trials, duration_s) before
%   running for a small test run.
%
%   See also mag_pair_cdf_rev1, gain_dist_table_rev1,
%   mag_pair_timeseries_rev1, scan_event_library_rev1, redteam_runs_rev1.

if ~exist("rev6_overrides", "var"), rev6_overrides = struct(); end
here = fileparts(which("report_numbers_rev6"));
out_dir = fullfile(here, "docs", "figures");
n_trials = 120;
duration_s = 3600;
if isfield(rev6_overrides, "out_dir"), out_dir = rev6_overrides.out_dir; end
if isfield(rev6_overrides, "n_trials"), n_trials = rev6_overrides.n_trials; end
if isfield(rev6_overrides, "duration_s"), duration_s = rev6_overrides.duration_s; end
if ~isfolder(out_dir), mkdir(out_dir); end

pct = [50 90 99 99.9 99.99 99.999];
seed_ts = 20261105;                              % Rev 5 long-run seed
grid0 = -0.1:0.05:0.1;                           % tables around the horizontal path
[scans, names] = nexrad_mode_scans_rev1();       % 216 modes
[~, names74] = nexrad_mode_scans_rev1("lowest");
cases = vcp_mode_cases_rev1();
labels = cellfun(@(c) string(c.label), cases(2:end));

% ---- NOAA case tilt: the ATC tilt with the largest 99.99% envelope
quick = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
o_roc = nexrad_asr11_mag_rev1(quick);
T_roc = gain_dist_table_rev1(o_roc.sys_nexrad.tx, scans, struct("name", "WSR-88D roc_spec", ...
    "mode_names", names, "el_grid_deg", grid0));
tilt = -1:0.5:2.5;
p9999 = zeros(size(tilt));
for i = 1:numel(tilt)
    q = quick;
    q.asr_tilt_deg = tilt(i);
    q.asr_remote_sl_dBi = 6.5;
    o = nexrad_asr11_mag_rev1(q);
    Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("el_grid_deg", grid0));
    R = mag_pair_cdf_rev1(Ta, T_roc, 0, 0, struct("pct", 99.99));
    p9999(i) = R.pct{1, 3};
end
writetable(table(tilt', p9999', 'VariableNames', ["asr_tilt_deg", "p99_99_dBi"]), ...
    fullfile(out_dir, "noaa_tilt_rev6.csv"));
[~, iw] = max(p9999);
noaa_tilt = tilt(iw);

% ---- Cases
cs = struct("name", ["roc_asr11" "roc_asr9" "nssl_asr11" "nssl_asr9" "noaa_asr11"], ...
    "pattern", ["roc_spec" "roc_spec" "nssl_measured" "nssl_measured" "roc_spec"], ...
    "atc", ["asr11" "asr9" "asr11" "asr9" "asr11"], ...
    "remote", [NaN NaN NaN NaN 6.5], "tilt", [0 0 0 0 noaa_tilt]);
n_c = numel(cs.name);
H = cell(n_c, 1);
W = cell(n_c, 1);
TS = cell(n_c, 1);
HM = cell(n_c, 1);
worst_label = strings(n_c, 1);
for c = 1:n_c
    q = quick;
    q.nexrad_pattern = cs.pattern(c);
    q.atc_radar = cs.atc(c);
    q.asr_tilt_deg = cs.tilt(c);
    q.asr_remote_sl_dBi = cs.remote(c);
    o = nexrad_asr11_mag_rev1(q);
    Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", grid0));
    Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("el_grid_deg", grid0));
    % Path ATC -> WSR-88D (g21, the WSR-88D as victim); tx = rx patterns,
    % so the other direction is the same distribution.
    R = mag_pair_cdf_rev1(Ta, Tn, 0, 0, struct("pct", pct, "per_mode", true));
    R74 = mag_pair_cdf_rev1(Ta, Tn, 0, 0, struct("pct", pct, "rx_modes", names74));
    env = R.pct{1, 3:2 + numel(pct)};
    env74 = R74.pct{1, 3:2 + numel(pct)};
    wm = extractAfter(R.worst_mode(1, :), " | ");
    H{c} = table(repmat(cs.name(c), numel(pct), 1), pct', env', env74', wm', ...
        repmat(R.pct.max_on_path_dBi(1), numel(pct), 1), repmat(R.bias_dB, numel(pct), 1), ...
        'VariableNames', ["case", "pct", "env216_dBi", "env74_dBi", "worst_mode", "max_on_path_dBi", "bias_dB"]);
    % Worst mode for the time series: the mode that sets the 99.9% envelope.
    worst_label(c) = wm(pct == 99.9);
    k = find(labels == worst_label(c), 1);
    if isempty(k)
        error("report:noCase", "Scan mode %s not found in vcp_mode_cases_rev1.", worst_label(c));
    end
    ow = rmfield(cases{k + 1}, "label");
    f = fieldnames(q);
    for i = 1:numel(f)
        ow.(f{i}) = q.(f{i});
    end
    ow.cdf_method = "exact";
    w = nexrad_asr11_mag_rev1(ow);
    if c == 1
        writetable(w.ref_levels, fullfile(out_dir, "ref_levels_rev6.csv"));
    end
    lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", seed_ts, ...
        "n_trials", n_trials, "duration_s", duration_s));
    S1 = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0, ...
        struct("radial_thr_dBi", [40 50 60 70], "radial_deg", 1));
    S05 = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0, ...
        struct("radial_thr_dBi", [40 50 60 70], "radial_deg", 0.5));
    P = S1.persistence;
    PV = S1.persistence_volume;
    P.Properties.VariableNames(2:end) = "rev_" + erase(string(P.Properties.VariableNames(2:end)), "rev_");
    PV.Properties.VariableNames(2:end) = "vol_" + erase(string(PV.Properties.VariableNames(2:end)), "vol_");
    TS{c} = [addvars(P, repmat(cs.name(c), height(P), 1), repmat(worst_label(c), height(P), 1), ...
        'Before', 1, 'NewVariableNames', ["case", "mode"]), PV(:, 2:end)];
    h1 = addvars(S1.radial_summary, repmat(1, height(S1.radial_summary), 1), 'Before', 1, ...
        'NewVariableNames', "radial_deg"); %#ok<REPMAT>
    h05 = addvars(S05.radial_summary, repmat(0.5, height(S05.radial_summary), 1), 'Before', 1, ...
        'NewVariableNames', "radial_deg");
    hm = [h1; h05];
    HM{c} = addvars(hm, repmat(cs.name(c), height(hm), 1), 'Before', 1, 'NewVariableNames', "case");
    fprintf("%s: worst mode %s, 99.9%% %.2f dBi, floor %.1f dBi\n", cs.name(c), worst_label(c), ...
        env(pct == 99.9), S1.floor21_dBi);
end
writetable(vertcat(H{:}), fullfile(out_dir, "headline_rev6.csv"));
writetable(vertcat(TS{:}), fullfile(out_dir, "timeseries_rev6.csv"));
writetable(vertcat(HM{:}), fullfile(out_dir, "harm_rev6.csv"));

% ---- Two-part model, worked pair (Rev 5 lookup example angles)
gridp = 1.5:0.05:2.5;
L = cell(2, 1);
atc = ["asr11" "asr9"];
for a = 1:2
    q = quick;
    q.atc_radar = atc(a);
    o = nexrad_asr11_mag_rev1(q);
    Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", gridp));
    Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("el_grid_deg", gridp));
    R = mag_pair_cdf_rev1(Tn, Ta, 2.0, 1.7, struct("pct", pct));
    L{a} = addvars(R.pct, atc(a), 'Before', 1, 'NewVariableNames', "atc");
end
writetable(vertcat(L{:}), fullfile(out_dir, "lookup_pair_rev6.csv"));
