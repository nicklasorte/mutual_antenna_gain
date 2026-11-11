%REPORT_NUMBERS_REV7 Tables and figures for methodology report rev7.
%   Runs the analyses of the revision 6 red team's way forward
%   (docs/rev6_redteam_plan_rev1.md, Section 4.2) under the revision 6
%   defaults and writes to docs/figures:
%     takeoff_surface_rev7.csv  envelope percentiles over 216 modes for
%                               takeoff angles -1 to 3 deg at each end
%                               (0.25 deg steps), each case (finding P1)
%     takeoff_worst_rev7.csv    largest value of each percentile over that
%                               range, and the angles where it occurs
%     bracket_cdf_rev7.csv      percentiles at 0/0 and at 0.5/0.5 deg for
%                               the pattern bracket (finding P4)
%     bracket_ts_rev7.csv       120 h event time series at the worst mode,
%                               same cases and angles: share of volumes
%                               reaching each level, longest run, runs cut
%                               by a trial end
%     bracket_harm_rev7.csv     harm metric for the same runs
%     vol_exceedance_rev7.csv   per-volume and per-revolution exceedance
%                               curves (finding P3)
%     drift_runs_rev7.csv       run length of volumes at 55, 60 and 70 dBi
%                               against the ATC rate-drift model, 10 trials
%                               of 24 h each (finding P2)
%     fig10_takeoff_curve.png ... fig14_harm.png
%   Cases: roc (ROC specification, ASR-11), roc_m3 and roc_m6 (the
%   specification lowered by 3 and 6 dB beyond the main beam), nssl (NSSL
%   measured), roc_asr9 (ASR-9, +5.5 dBi remote sidelobe).
%
%   Set rev7_overrides (struct: out_dir, n_trials, duration_s,
%   drift_trials, drift_hours, grid_step) before running for a small run.
%
%   Style: NTIA palette (mid gray context, light blue highlight, purple
%   second highlight, dark blue text), Times New Roman, 6.5 in wide.
%   Titles, captions and sources are in the report text, not the images.
%
%   See also report_numbers_rev6, mag_pair_cdf_rev1, mag_pair_timeseries_rev1.

if ~exist("rev7_overrides", "var"), rev7_overrides = struct(); end
here = fileparts(which("report_numbers_rev7"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "n_trials", 120, "duration_s", 3600, ...
    "drift_trials", 10, "drift_hours", 24, "grid_step", 0.25);
f = fieldnames(rev7_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev7_overrides.(f{i});
end
if ~isfolder(ov.out_dir), mkdir(ov.out_dir); end
out_dir = ov.out_dir;
pct = [50 90 99 99.9 99.99 99.999];
c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "purple", [110 79 146]/255, "gray", [138 141 145]/255, "grid", [230 230 230]/255);

% ---- Systems and tables for each case
cases = ["roc", "roc_m3", "roc_m6", "nssl", "roc_asr9"];
copt = {struct(), struct("nexrad_sl_margin_dB", 3), struct("nexrad_sl_margin_dB", 6), ...
    struct("nexrad_pattern", "nssl_measured"), struct("atc_radar", "asr9")};
[scans, names] = nexrad_mode_scans_rev1();
g_el = -1:0.05:3;
sys = cell(numel(cases), 1);
Tn = cell(numel(cases), 1);
Ta = cell(numel(cases), 1);
for k = 1:numel(cases)
    q = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1);
    fk = fieldnames(copt{k});
    for i = 1:numel(fk)
        q.(fk{i}) = copt{k}.(fk{i});
    end
    sys{k} = nexrad_asr11_mag_rev1(q);
    Tn{k} = gain_dist_table_rev1(sys{k}.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", g_el));
    Ta{k} = gain_dist_table_rev1(sys{k}.sys_asr11.tx, sys{k}.sys_asr11.scan, struct("el_grid_deg", g_el));
end

% ---- A. Takeoff-angle surface (finding P1)
e = (-1:ov.grid_step:3)';
[EN, EA] = ndgrid(e, e);
S = cell(numel(cases), 1);
W = cell(numel(cases), 1);
for k = 1:numel(cases)
    R = mag_pair_cdf_rev1(Ta{k}, Tn{k}, EA(:), EN(:), struct("pct", pct));   % path ATC -> WSR-88D
    P = R.pct{:, 3:2 + numel(pct)};
    S{k} = [table(repmat(cases(k), numel(EN), 1), EN(:), EA(:), 'VariableNames', ["case", "el_nexrad_deg", "el_atc_deg"]), ...
        array2table(P, "VariableNames", "p" + replace(string(pct), ".", "_"))];
    [mx, im] = max(P, [], 1);
    z = P(EN(:) == 0 & EA(:) == 0, :);
    W{k} = table(repmat(cases(k), numel(pct), 1), pct', z', mx', EN(im(:)), EA(im(:)), mx' - z', ...
        'VariableNames', ["case", "pct", "at_0_0_dBi", "worst_dBi", "worst_el_nexrad_deg", "worst_el_atc_deg", "increase_dB"]);
end
T_surf = vertcat(S{:});
T_worst = vertcat(W{:});
writetable(T_surf, fullfile(out_dir, "takeoff_surface_rev7.csv"));
writetable(T_worst, fullfile(out_dir, "takeoff_worst_rev7.csv"));

% ---- B. Pattern bracket: CDFs and time series at 0/0 and 0.5/0.5 deg
ang = [0 0; 0.5 0.5];
B = cell(numel(cases), 1);
for k = 1:numel(cases)
    R = mag_pair_cdf_rev1(Ta{k}, Tn{k}, ang(:, 2), ang(:, 1), struct("pct", pct));
    B{k} = [table(repmat(cases(k), 2, 1), ang(:, 1), ang(:, 2), 'VariableNames', ["case", "el_nexrad_deg", "el_atc_deg"]), ...
        R.pct(:, 3:end), table(extractAfter(R.worst_mode(:, pct == 99.9), " | "), 'VariableNames', "worst_mode_99_9")];
end
T_bcdf = vertcat(B{:});
writetable(T_bcdf, fullfile(out_dir, "bracket_cdf_rev7.csv"));

% Worst mode: VCP 212 SAILSx3 AVSET 6.4 sets the 99.9% envelope at both
% angles in every case (checked here); one event library serves them all.
mode_ts = "VCP 212 SAILSx3 AVSET 6.4";
if any(T_bcdf.worst_mode_99_9 ~= mode_ts)
    warning("report:worstMode", "The 99.9%% worst mode differs from %s in some cases; the time series uses %s.", mode_ts, mode_ts);
end
mc = vcp_mode_cases_rev1();
lab = cellfun(@(x) string(x.label), mc);
mo = rmfield(mc{find(lab == mode_ts, 1)}, "label");
thr = [40 45 50 55 60 65 70];
TS = cell(numel(cases), 1);
HM = cell(numel(cases), 1);
EX = cell(numel(cases), 1);
for k = 1:numel(cases)
    q = mo;
    q.seed = 1; q.n_samples = 1e3; q.duration_s = 60; q.n_trials = 1;
    fk = fieldnames(copt{k});
    for i = 1:numel(fk)
        q.(fk{i}) = copt{k}.(fk{i});
    end
    w = nexrad_asr11_mag_rev1(q);
    lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, ...
        "n_trials", ov.n_trials, "duration_s", ov.duration_s));
    Rp = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, ang(:, 1), ang(:, 2), ...
        struct("persist_thr_dBi", thr, "radial_thr_dBi", [40 50 60 70]));
    for a = 1:2
        PV = Rp(a).persistence_volume;
        PR = Rp(a).persistence;
        TS{k}{a} = table(repmat(cases(k), numel(thr), 1), repmat(ang(a, 1), numel(thr), 1), thr', PR.rev_pct, ...
            PV.vol_pct, PV.max_run_vols, PV.max_run_s/60, PV.censored_runs, PV.n_runs, ...
            'VariableNames', ["case", "el_deg", "threshold_dBi", "rev_pct", "vol_pct", "max_run_vols", ...
            "max_run_min", "censored_runs", "n_runs"]);
        H = Rp(a).radial_summary;
        HM{k}{a} = addvars(H, repmat(cases(k), height(H), 1), repmat(ang(a, 1), height(H), 1), ...
            'Before', 1, 'NewVariableNames', ["case", "el_deg"]);
        X = Rp(a).exceedance;
        EX{k}{a} = addvars(X, repmat(cases(k), height(X), 1), repmat(ang(a, 1), height(X), 1), ...
            'Before', 1, 'NewVariableNames', ["case", "el_deg"]);
    end
end
T_bts = vertcat(TS{:}); T_bts = vertcat(T_bts{:});
T_harm = vertcat(HM{:}); T_harm = vertcat(T_harm{:});
T_ex = vertcat(EX{:}); T_ex = vertcat(T_ex{:});
writetable(T_bts, fullfile(out_dir, "bracket_ts_rev7.csv"));
writetable(T_harm, fullfile(out_dir, "bracket_harm_rev7.csv"));
writetable(T_ex, fullfile(out_dir, "vol_exceedance_rev7.csv"));

% ---- C. Run length against the drift model (finding P2), ROC, ASR-11
q = mo;
q.seed = 1; q.n_samples = 1e3; q.duration_s = 60; q.n_trials = 1;
w = nexrad_asr11_mag_rev1(q);
models = ["per_trial", "random_walk", "random_walk", "random_walk", "random_walk"];
steps = [NaN 0.005 0.01 0.02 0.05];
thr_d = [55 60 70];
D = cell(numel(models), 1);
for m = 1:numel(models)
    spec = w.scan_max.clock_spec;
    spec.jitter_mode = models(m);
    if models(m) == "random_walk"
        spec.drift_step_a = steps(m);
    end
    lib = scan_event_library_rev1(spec, struct("seed", 20261111, "n_trials", ov.drift_trials, ...
        "duration_s", 3600*ov.drift_hours));
    Rd = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, ang(:, 1), ang(:, 2), ...
        struct("persist_thr_dBi", thr_d, "radial_thr_dBi", 60, "keep_scans", true));
    rows = {};
    for a = 1:2
        V = Rd(a).nexrad_volume(Rd(a).nexrad_volume.complete, :);
        for i = 1:numel(thr_d)
            len = run_lengths(V, thr_d(i));
            PV = Rd(a).persistence_volume;
            rows(end + 1, :) = {models(m), steps(m), ang(a, 1), thr_d(i), PV.vol_pct(i), numel(len), ...
                median([len; NaN], "omitnan"), prctile([len; NaN], 95), max([len; 0]), ...
                max([len; 0])*w.scan_max.t_volume_s/60, PV.censored_runs(i)}; %#ok<SAGROW>
        end
    end
    D{m} = cell2table(rows, 'VariableNames', ["drift_model", "step_deg_s_per_rev", "el_deg", "threshold_dBi", ...
        "vol_pct", "n_runs", "median_run_vols", "p95_run_vols", "max_run_vols", "max_run_min", "censored_runs"]);
end
T_drift = vertcat(D{:});
writetable(T_drift, fullfile(out_dir, "drift_runs_rev7.csv"));

% ---- Figures
% Figure 10: 99.9% and 99.999% envelope against takeoff angle (same angle
% at both ends).
f = new_fig(6.5, 2.6);
tl = tiledlayout(f, 1, 2, "TileSpacing", "compact", "Padding", "compact");
pk = ["p99_9", "p99_999"];
ttl = ["99.9% of the time", "99.999% of the time"];
for j = 1:2
    ax = nexttile(tl);
    hold(ax, "on");
    for k = [4 3 1]
        Sk = S{k}(S{k}.el_nexrad_deg == S{k}.el_atc_deg, :);
        col = c.gray;
        lw = 1.2;
        ls = "-";
        if cases(k) == "roc", col = c.blue; lw = 2; end
        if cases(k) == "roc_m6", ls = "--"; end
        plot(ax, Sk.el_nexrad_deg, Sk.(pk(j)), ls, "Color", col, "LineWidth", lw);
    end
    xline(ax, 0, ":", "Color", c.dark);
    style(ax, c, "Takeoff angle at both ends (deg)", "MAG (dBi)", ttl(j));
end
put_label(ax, 0.8, 74, "All three patterns", c.blue, true);
ax = tl.Children(end);                          % first tile
put_label(ax, 1.55, 47.6, "ROC specification", c.blue, true);
put_label(ax, 0.9, 38.6, "Specification - 6 dB, dashed", c.gray, false);
put_label(ax, 0.9, 36.6, "NSSL measured, solid", c.gray, false);
save_fig(f, out_dir, "fig10_takeoff_curve");

% Figure 11: exceedance of the MAG at 0.5/0.5 deg for the pattern bracket.
f = new_fig(6.5, 3.0);
ax = axes(f);
hold(ax, "on");
for k = [4 3 2 1]
    R = mag_pair_cdf_rev1(Ta{k}, Tn{k}, 0.5, 0.5, struct("return_curve", true));
    cv = R.curve{1};
    col = c.gray;
    lw = 1.2;
    if cases(k) == "roc", col = c.blue; lw = 2; end
    semilogy(ax, cv.gain_dBi, cv.exceed_pct, "Color", col, "LineWidth", lw);
end
ax.YScale = "log";
ylim(ax, [1e-4 100]);
style(ax, c, "MAG (dBi)", "Time exceeded (%, log scale)", "");
put_label(ax, 51.5, 0.03, ["Step at 50.5 dBi: ATC main beam", "+ flat specification sidelobe"], c.blue, true);
put_label(ax, 0, 2e-3, "Gray: specification lowered 3 dB and 6 dB, and NSSL measured", c.gray, false);
save_fig(f, out_dir, "fig11_pattern_bracket");

% Figure 12: share of WSR-88D volumes whose maximum reaches each level.
f = new_fig(6.5, 2.8);
ax = axes(f);
hold(ax, "on");
for k = [4 3 2 5 1]
    X = T_ex(T_ex.case == cases(k) & T_ex.el_deg == 0, :);
    col = c.gray;
    lw = 1.2;
    if cases(k) == "roc", col = c.blue; lw = 2; end
    if cases(k) == "roc_asr9", col = c.purple; end
    plot(ax, X.gain_dBi, X.vol_pct, "Color", col, "LineWidth", lw);
end
xlim(ax, [40 80]);
style(ax, c, "Largest MAG in the volume (dBi)", "Volumes reaching it (%)", "");
put_label(ax, 50.7, 90, "ASR-11, specification", c.blue, true);
put_label(ax, 49.8, 60, "ASR-9", c.purple, true, "right");
put_label(ax, 45.3, 60, "Specification - 6 dB, NSSL", c.gray, false, "right");
save_fig(f, out_dir, "fig12_volume_exceedance");

% Figure 13: run lengths at 60 dBi (path at 0.5 deg) for each drift model:
% median, 95th percentile and longest run.
f = new_fig(6.5, 2.6);
ax = axes(f);
hold(ax, "on");
Dd = T_drift(T_drift.threshold_dBi == 60 & T_drift.el_deg == 0.5, :);
xl = ["Fixed offset" , "Walk " + compose("%g", Dd.step_deg_s_per_rev(2:end)')];
x = (1:height(Dd))';
plot(ax, x, Dd.max_run_vols, "o", "Color", c.blue, "MarkerFaceColor", c.blue, "MarkerSize", 6);
plot(ax, x, Dd.p95_run_vols, "s", "Color", c.gray, "MarkerFaceColor", c.gray, "MarkerSize", 5);
plot(ax, x, Dd.median_run_vols, "d", "Color", c.gray, "MarkerSize", 5);
ax.XTick = x;
ax.XTickLabel = xl;
xlim(ax, [0.5 height(Dd) + 0.5]);
ylim(ax, [0 max(Dd.max_run_vols) + 2]);
style(ax, c, "ATC rate drift model (random-walk step in deg/s per revolution)", "Run of volumes at 60 dBi", "");
put_label(ax, 1.15, Dd.max_run_vols(1), "Longest", c.blue, true);
put_label(ax, 1.15, Dd.p95_run_vols(1), "95th percentile", c.gray, false);
put_label(ax, 1.15, Dd.median_run_vols(1) + 0.3, "Median", c.gray, false);
save_fig(f, out_dir, "fig13_runs_vs_drift");

% Figure 14: share of 0.5 deg cuts with a radial at or above each level.
f = new_fig(6.5, 2.6);
ax = axes(f);
hold(ax, "on");
for k = [4 3 2 5 1]
    H = T_harm(T_harm.case == cases(k) & T_harm.el_deg == 0, :);
    col = c.gray;
    lw = 1.2;
    if cases(k) == "roc", col = c.blue; lw = 2; end
    if cases(k) == "roc_asr9", col = c.purple; end
    plot(ax, H.threshold_dBi, H.cuts_with_hit_pct, "-o", "Color", col, "MarkerFaceColor", col, "LineWidth", lw);
end
ax.YScale = "log";
style(ax, c, "MAG level (dBi)", "0.5° cuts with a radial at the level (%)", "");
put_label(ax, 51, 30, "ROC specification, ASR-11", c.blue, true);
put_label(ax, 50.5, 1.5, ["Specification - 3 dB and - 6 dB,", "NSSL measured, and ASR-9"], c.gray, false);
save_fig(f, out_dir, "fig14_harm");

% =========================================================================
function len = run_lengths(V, thr)
% Lengths of runs of consecutive complete volumes reaching thr, per trial.
len = [];
for j = unique(V.trial)'
    hit = V.max_g21_dBi(V.trial == j) >= thr;
    d = diff([0; hit(:); 0]);
    len = [len; find(d == -1) - find(d == 1)]; %#ok<AGROW>
end
end

function f = new_fig(w_in, h_in)
f = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(f, "light");
end
f.Position(3:4) = [w_in h_in];
end

function style(ax, c, xl, yl, ttl)
box(ax, "off");
ax.Color = "w";
ax.FontName = "Times New Roman";
ax.FontSize = 9;
ax.TickDir = "out";
ax.XColor = c.dark;
ax.YColor = c.dark;
ax.YGrid = "on";
ax.GridColor = c.grid;
ax.GridAlpha = 1;
xlabel(ax, xl);
ylabel(ax, yl);
if ttl ~= ""
    title(ax, ttl, "FontWeight", "bold", "Color", c.dark, "FontSize", 9);
    ax.TitleHorizontalAlignment = "left";
end
end

function put_label(ax, x, y, s, col, bold, halign)
if nargin < 7
    halign = "left";
end
w = "normal";
if bold
    w = "bold";
end
text(ax, x, y, s, "Color", col, "FontWeight", w, "FontName", "Times New Roman", ...
    "FontSize", 8, "HorizontalAlignment", halign);
end

function save_fig(f, d, name)
exportgraphics(f, fullfile(d, name + ".png"), "Resolution", 300);
close(f);
end
