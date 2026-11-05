%REPORT_FIGURES_REV5 Figures and tables for methodology report rev5.
%   Runs example_nexrad_asr11_rev1 (WSR-88D vs ASR-11, all VCP modes) and
%   writes to docs/figures:
%     fig1_antenna_patterns.png   NEXRAD and ASR-11 antenna patterns
%     fig2_vcp_schedule.png       elevation versus time in one volume
%     fig3_low_tilt_driver.png    99.9% MAG versus time at 0.5 deg
%     fig4_exceedance.png         MAG exceedance with reference levels
%     fig5_max_per_rev.png        maximum MAG per ASR-11 revolution
%     fig7_mbmb_vs_statistical.png  main beam to main beam versus statistical MAG
%     fig6_peak_anatomy.png       gains around the two largest peaks (120 h run, the
%                                 trial with the most peaks above the 45.5 dBi plateau)
%     fig8_sensitivity.png        MAG change for each antenna and scan assumption
%     fig9_asr9_comparison.png    ASR-9 and ASR-11 azimuth patterns and exceedance
%     asr9_*_rev5.csv             ASR-9 sweep, envelope, comparison with the ASR-11,
%                                 reference levels, lobe pairs, persistence, MAG
%                                 choices and assumptions (example_asr9_rev1)
%     plateaus_rev5.csv, peaks_above_plateau_120h_rev5.csv, rev_exceedance_120h_rev5.csv,
%     asr_rev_pct_120h_rev5.csv, convergence_rev5.csv
%     lobe_pairs_rev5.csv         time share and largest MAG per pattern-region pair
%     persistence_120h_rev5.csv   runs of consecutive high ASR-11 revolutions
%     persistence_volume_120h_rev5.csv  runs of consecutive high NEXRAD volumes
%     sensitivity_sweep_rev5.csv  one-at-a-time sensitivity (example_sensitivity_rev1)
%     alignment_dwell_rev5.csv    equal and multiple scan rates (example_asr11_drift_rev1)
%     vcp_sweep_rev5.csv          sweep table (all modes)
%     ref_levels_rev5.csv, by_elevation_rev5.csv, asr_rev_pct_rev5.csv,
%     sensitivity_rev5.csv, mode_envelope_rev5.csv, mbmb_vs_statistical_rev5.csv,
%     cdf_array_rev5.csv (CDF arrays: percent 0-100 and MAG dBi)
%     Tables typed in the report, exported here so every number has a file:
%     case1_percentiles_rev5.csv  Case 1 of example_radar_mutual_gain_rev1 (Section 12)
%     vcp_table_rev5.csv          VCP use, cuts, volume time, time at 0.5 deg (Step 5)
%     worst_mode_weights_rev5.csv time weights of the worst mode by elevation (Step 5)
%     nexrad_gain_by_cut_rev5.csv NEXRAD gain toward the ASR-11, beam azimuth on the path (Step 2)
%     asr11_gain_cut_rev5.csv     ASR-11 gain toward the NEXRAD by azimuth offset (Step 3)
%     fig_a1_call_tree.png        call tree of the repository functions (Appendix A),
%                                 from matlab.codetools.requiredFilesAndProducts
%   Style: NTIA palette (mid gray context, light blue highlight, purple
%   second highlight, dark blue text), Times New Roman, 6.5 in wide.
%   Titles, captions and sources are in the report text, not the images.
%   A workspace struct example_overrides makes a small run (tests): its
%   n_samples, duration_s, n_trials apply to every run, and out_dir, n_h
%   (trials of the long run), long_duration_s (seconds per long-run trial)
%   and n_list (convergence sample sizes) replace the defaults.

% ---- Case 1 of example_radar_mutual_gain_rev1 (Section 12 table). It runs
%      first, because each example script clears the workspace except
%      example_overrides.
example_radar_mutual_gain_rev1;
if ~exist("example_overrides", "var")
    example_overrides = struct();
end
writetable(res{1}.pct21, fullfile(report_dir(example_overrides), "case1_percentiles_rev5.csv"));

example_nexrad_asr11_rev1;          % defines T, V, res, w, i999, i99
close all;

out_dir = report_dir(example_overrides);
n_h = 120;
n_list = [5e5 5e6 2e7];
if isfield(example_overrides, "n_h"), n_h = example_overrides.n_h; end
if isfield(example_overrides, "n_list"), n_list = example_overrides.n_list; end
writetable(T, fullfile(out_dir, "vcp_sweep_rev5.csv"));
writetable(w.ref_levels, fullfile(out_dir, "ref_levels_rev5.csv"));
writetable(w.scan_max.by_elevation, fullfile(out_dir, "by_elevation_rev5.csv"));
writetable(w.scan_max.asr_rev_pct, fullfile(out_dir, "asr_rev_pct_rev5.csv"));
writetable(w.sensitivity, fullfile(out_dir, "sensitivity_rev5.csv"));

% ---- Worst mode (99.9% MAG) options, from the sweep case list
i_w = 1 + i999;                     % worst mode in res, T and c
wopt = rmfield(c{i_w}, "label");
wopt.seed = base.seed;
wopt.n_samples = base.n_samples;

% ---- Worst mode, long time series: 120 trials of 1 h (5 days, about
%      90,000 ASR-11 revolutions; a run of two to five days was asked for at
%      the November 4, 2026 meeting). Same seed as the example, so trial 1
%      and the random-mode CDF match the sweep.
o_long = wopt;
o_long.duration_s = 3600;
if exist("example_overrides", "var")
    o_long = merge_overrides_rev1(o_long, example_overrides);
    if isfield(example_overrides, "long_duration_s")
        o_long.duration_s = example_overrides.long_duration_s;
    end
end
o_long.n_trials = n_h;
wl = nexrad_asr11_mag_rev1(o_long);
Rl = wl.scan_max.asr_rev(wl.scan_max.asr_rev.complete, :);
% Trial shown in Figures 5 and 6: the hour with the most peaks above the
% 45.5 dBi plateau (ties: the largest peak). Rerun only that trial to keep
% its samples; trial j always uses substream 1000 + j, so it repeats.
lev_hi = wl.ref_levels.mag_dBi(3) + 1;
[G, tr] = findgroups(Rl.trial);
n_pk = splitapply(@(g) sum(g >= lev_hi), Rl.max_g21_dBi, G);
m_pk = splitapply(@max, Rl.max_g21_dBi, G);
[~, k_tr] = sortrows([n_pk, m_pk], [-1 -2]);
j_fig = tr(k_tr(1));
o_fig = o_long;
o_fig.first_trial = j_fig;
o_fig.n_trials = 1;
o_fig.n_samples = 1e4;                          % the random-mode CDF is not used
o_fig.keep_ts1 = true;
wf = nexrad_asr11_mag_rev1(o_fig);
writetable(wl.scan_max.plateaus, fullfile(out_dir, "plateaus_rev5.csv"));
writetable(wl.scan_max.asr_rev_pct, fullfile(out_dir, "asr_rev_pct_120h_rev5.csv"));
writetable(wl.lobe_pairs, fullfile(out_dir, "lobe_pairs_rev5.csv"));
writetable(wl.scan_max.persistence, fullfile(out_dir, "persistence_120h_rev5.csv"));
writetable(wl.scan_max.persistence_volume, fullfile(out_dir, "persistence_volume_120h_rev5.csv"));
pk = Rl(Rl.max_g21_dBi >= lev_hi, ["trial", "rev", "t_at_max_s", "max_g21_dBi", ...
    "nexrad_el_at_max_deg", "nexrad_gain_at_max_dBi", "nexrad_az_off_at_max_deg", ...
    "asr11_gain_at_max_dBi", "asr11_az_off_at_max_deg"]);
writetable(sortrows(pk, "max_g21_dBi", "descend"), fullfile(out_dir, "peaks_above_plateau_120h_rev5.csv"));
lv = [wl.ref_levels.mag_dBi(5), 30, 40, wl.ref_levels.mag_dBi(3), 50, 60, 70, wl.ref_levels.mag_dBi(1) - 0.2];
exc = table(lv(:), arrayfun(@(v) sum(Rl.max_g21_dBi >= v - 0.2), lv(:)), ...
    arrayfun(@(v) 100*mean(Rl.max_g21_dBi >= v - 0.2), lv(:)), ...
    'VariableNames', ["level_dBi", "revolutions", "pct_of_revolutions"]);
exc.per_hour = exc.revolutions/(n_h*o_long.duration_s/3600);
writetable(exc, fullfile(out_dir, "rev_exceedance_120h_rev5.csv"));

% ---- Sample-size convergence of the time-weighted percentiles (worst mode)
cv = cell(numel(n_list), 1);
for k = 1:numel(n_list)
    o = wopt;
    o.n_samples = n_list(k);
    o.duration_s = min(60, base.duration_s);
    o.n_trials = 1;
    r = nexrad_asr11_mag_rev1(o);
    q = r.cdf.pct21(r.cdf.pct21.pct >= 99, :);
    cv{k} = [table(n_list(k) + zeros(height(q), 1), 'VariableNames', "n_samples"), q];
end
writetable(vertcat(cv{:}), fullfile(out_dir, "convergence_rev5.csv"));
writetable(env.pct, fullfile(out_dir, "mode_envelope_rev5.csv"));

% ---- Sensitivity (one assumption at a time, from the worst mode) and
%      alignment of equal or multiple scan rates
sens_base = wopt;
sens_base.duration_s = base.duration_s;
sens_base.n_trials = base.n_trials;
example_sensitivity_rev1;           % defines Tsens
writetable(Tsens, fullfile(out_dir, "sensitivity_sweep_rev5.csv"));
example_asr11_drift_rev1;           % defines Tdrift
writetable(Tdrift, fullfile(out_dir, "alignment_dwell_rev5.csv"));

% ---- ASR-9 (M.1464-2 Radar B): every VCP mode, envelope, and the paired
%      comparison with the ASR-11 at the ASR-9's worst mode
asr9_base = struct("seed", base.seed, "n_samples", base.n_samples, ...
    "duration_s", base.duration_s, "n_trials", base.n_trials);
example_asr9_rev1;                  % defines T9, env9, w9, w11, Tcmp9, cmp9_label
writetable(T9, fullfile(out_dir, "asr9_vcp_sweep_rev5.csv"));
writetable(env9.pct, fullfile(out_dir, "asr9_mode_envelope_rev5.csv"));
writetable(Tcmp9, fullfile(out_dir, "asr9_comparison_rev5.csv"));
writetable(w9.ref_levels, fullfile(out_dir, "asr9_ref_levels_rev5.csv"));
writetable(w9.lobe_pairs, fullfile(out_dir, "asr9_lobe_pairs_rev5.csv"));
writetable(w9.scan_max.persistence, fullfile(out_dir, "asr9_persistence_rev5.csv"));
writetable(w9.scan_max.persistence_volume, fullfile(out_dir, "asr9_persistence_volume_rev5.csv"));
writetable(w9.sensitivity, fullfile(out_dir, "asr9_sensitivity_rev5.csv"));
writetable(w9.assumptions, fullfile(out_dir, "asr9_assumptions_rev5.csv"));

c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "purple", [110 79 146]/255, "gray", [138 141 145]/255, "grid", [230 230 230]/255);
i_p = 1;                            % paper baseline
i_215 = find(T.case == "VCP 215");
worst = T.case(i_w);

% ---- Figure 1: antenna patterns ------------------------------------------
A = res{i_w};
pn = antenna_prepare_rev1(A.sys_nexrad.tx);
pr = antenna_prepare_rev1(nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, ...
    "duration_s", 60, "n_trials", 1, "nexrad_pattern", "roc_spec")).sys_nexrad.tx);
pf = antenna_prepare_rev1(struct("type", "f1245", "g_max_dBi", 46));
pa = antenna_prepare_rev1(A.sys_asr11.tx);
paf = antenna_prepare_rev1(struct("type", "f1245", "g_max_dBi", 34));
x = logspace(-2, log10(180), 2000)';
st = state_az(x, 0);
f = new_fig(6.5, 2.9);
tl = tiledlayout(f, 1, 3, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(tl);
semilogx(ax, x, antenna_gain_rev1(pf, 0, 0, st), ":", "Color", c.gray, "LineWidth", 1.2); hold(ax, "on");
semilogx(ax, x, antenna_gain_rev1(pr, 0, 0, st), "--", "Color", c.gray, "LineWidth", 1.2);
semilogx(ax, x, antenna_gain_rev1(pn, 0, 0, st), "-", "Color", c.blue, "LineWidth", 2);
put_label(ax, 0.012, 2, "NSSL measured (default)", c.blue, true);
put_label(ax, 20, 10, "ROC spec", c.gray, false);
put_label(ax, 150, -16, "F.1245 46 dBi", c.gray, false, "right");
style(ax, c, "Off-axis angle (deg, log scale)", "Gain (dBi)", "a. WSR-88D");
ylim(ax, [-20 50]); xlim(ax, [0.01 180]);
ax = nexttile(tl);
semilogx(ax, x, antenna_gain_rev1(paf, 0, 0, st), ":", "Color", c.gray, "LineWidth", 1.2); hold(ax, "on");
semilogx(ax, x, antenna_gain_rev1(pa, 0, 0, st), "-", "Color", c.blue, "LineWidth", 2);
put_label(ax, 0.012, 2, "M.1851-2 cosine taper", c.blue, true);
put_label(ax, 150, -14, "F.1245 34 dBi", c.gray, false, "right");
style(ax, c, "Azimuth offset (deg, log scale)", "", "b. ASR-11 azimuth");
ylim(ax, [-20 50]); xlim(ax, [0.01 180]);
ax = nexttile(tl);
e = (-10:0.05:40)';
plot(ax, e, antenna_gain_rev1(pa, 0, e, state_az(zeros(size(e)), 0)), "-", "Color", c.blue, "LineWidth", 2);
hold(ax, "on");
put_label(ax, -9, 45, "Cosecant squared to +30 deg", c.blue, true);
style(ax, c, "Elevation (deg), 0 deg tilt", "", "c. ASR-11 elevation");
ylim(ax, [-20 50]);
save_fig(f, out_dir, "fig1_antenna_patterns");

% ---- Figure 2: VCP elevation schedule ------------------------------------
f = new_fig(6.5, 2.8);
ax = axes(f);
sw = res{i_w}.sys_nexrad.scan;
s215 = res{i_215}.sys_nexrad.scan;
[tx, ey] = steps(s215);
stairs(ax, tx, ey, "-", "Color", c.gray, "LineWidth", 1.2); hold(ax, "on");
[tx, ey] = steps(sw);
stairs(ax, tx, ey, "-", "Color", c.blue, "LineWidth", 2);
put_label(ax, 5, 16.3, "VCP 215", c.gray, false);
low = sw.cuts.el_deg == 0.5 & sw.cuts.waveform == "CS";
t0 = cumsum([0; sw.cuts.dur_s(1:end-1)]);
for k = find(low)'
    text(ax, t0(k) + 3, 2.2, "0.5 deg", "Color", c.blue, "FontWeight", "bold", ...
        "FontName", "Times New Roman", "FontSize", 8);
end
put_label(ax, 5, 18.3, string(worst), c.blue, true);
style(ax, c, "Time in volume (s)", "Elevation (deg)", "");
ylim(ax, [0 20]); xlim(ax, [0 max(sw.volume_s, s215.volume_s)]);
save_fig(f, out_dir, "fig2_vcp_schedule");

% ---- Figure 3: MAG versus time at 0.5 deg ---------------------------------
f = new_fig(6.5, 3.0);
ax = axes(f);
[tl, o] = sort(V.t_low_pct);
plot(ax, tl, V.p99_99(o), "-", "Color", c.gray, "LineWidth", 1.2); hold(ax, "on");
plot(ax, tl, V.p99(o), "--", "Color", c.gray, "LineWidth", 1.2);
plot(ax, tl, V.p99_9(o), "-", "Color", c.blue, "LineWidth", 2.2);
yline(ax, T.p99_9(i_p), ":", "Color", c.purple, "LineWidth", 1.2);
put_label(ax, tl(end), V.p99_99(o(end)) + 1.8, "99.99%", c.gray, false, "right");
put_label(ax, tl(end), V.p99_9(o(end)) - 2.2, "99.9%", c.blue, true, "right");
put_label(ax, tl(end), V.p99(o(end)) + 1.8, "99%", c.gray, false, "right");
put_label(ax, tl(1), T.p99_9(i_p) + 1.6, sprintf("Paper baseline 99.9%%, %.1f dBi", T.p99_9(i_p)), c.purple, true);
style(ax, c, sprintf("NEXRAD time at 0.5 deg (%% of volume), %d VCP modes", height(V)), "Mutual antenna gain (dBi)", "");
ylim(ax, [10 50]); xlim(ax, [tl(1) - 1, tl(end) + 1]);
save_fig(f, out_dir, "fig3_low_tilt_driver");

% ---- Figure 4: exceedance with reference levels ----------------------------
f = new_fig(6.5, 3.8);
ax = axes(f);
cw = res{i_w}.cdf.cdf21; c215 = res{i_215}.cdf.cdf21; cp = res{i_p}.cdf.cdf21;
semilogy(ax, cp.gain_dBi, cp.exceed_pct, "--", "Color", c.purple, "LineWidth", 1.4); hold(ax, "on");
semilogy(ax, c215.gain_dBi, c215.exceed_pct, "-", "Color", c.gray, "LineWidth", 1.4);
semilogy(ax, cw.gain_dBi, cw.exceed_pct, "-", "Color", c.blue, "LineWidth", 2.2);
semilogy(ax, env.exceed.gain_dBi, max(env.exceed.envelope_exceed_pct, 1e-6), "--", "Color", c.dark, "LineWidth", 1.0);
R = w.ref_levels;
tag = ["a", "b", "c", "d", "e"];           % defined in the figure caption
y_tag = [40, 40, 40, 12, 40];              % stagger the 44.5/45.5 dBi pair
for k = 1:height(R)
    xline(ax, R.mag_dBi(k), ":", "Color", c.dark, "LineWidth", 0.8);
    text(ax, R.mag_dBi(k), y_tag(k), tag(k), "Color", c.dark, "FontName", "Times New Roman", ...
        "FontSize", 9, "FontWeight", "bold", "HorizontalAlignment", "center", "BackgroundColor", "w");
end
put_label(ax, -18, 30, "Paper baseline", c.purple, true);
put_label(ax, 2, 0.6, "VCP 215", c.gray, false);
put_label(ax, 27, 0.6, string(worst), c.blue, true);
put_label(ax, 60, 2.5e-4, "Dashed: envelope of all " + height(V) + " VCP modes", c.dark, false, "right");
style(ax, c, "Mutual antenna gain (dBi)", "Time exceeded (%, log scale)", "");
ylim(ax, [1e-4 100]); xlim(ax, [-25 82]);
ax.YGrid = "on";
save_fig(f, out_dir, "fig4_exceedance");

% ---- Figure 5: maximum MAG per ASR-11 revolution, trial 1, peaks labeled ----
f = new_fig(6.5, 3.4);
ax = axes(f);
Rv = wf.scan_max.asr_rev(wf.scan_max.asr_rev.complete, :);   % trial j_fig
R = wl.ref_levels;
P = wl.scan_max.plateaus;
plot(ax, Rv.t_start_s, Rv.max_g21_dBi, "-", "Color", c.blue, "LineWidth", 0.9); hold(ax, "on");
lev09 = P.level_dBi(P.el_deg == 0.9);
for v = [R.mag_dBi(5), lev09, R.mag_dBi(3), R.mag_dBi(1)]
    yline(ax, v, ":", "Color", c.dark);
end
put_label(ax, 20, R.mag_dBi(5) - 3, sprintf("ASR-11 main beam + NEXRAD back lobe, %.1f dBi (every revolution)", R.mag_dBi(5)), c.dark, false, "left", "w");
put_label(ax, 3580, lev09 - 2.2, sprintf("NEXRAD main beam at 0.9 deg + ASR-11 back lobe, %.1f dBi", lev09), c.dark, false, "right", "w");
put_label(ax, 3580, R.mag_dBi(3) - 2.2, sprintf("NEXRAD main beam at 0.5 deg + ASR-11 back lobe, %.1f dBi (%.1f%% of revolutions)", ...
    R.mag_dBi(3), P.simulated_rev_pct(P.el_deg == 0.5)), c.dark, false, "right", "w");
put_label(ax, 3580, R.mag_dBi(1) + 1.6, sprintf("Main beam to main beam on the path, %.1f dBi", R.mag_dBi(1)), c.dark, false, "right", "w");
% Label every peak between the 45.5 dBi plateau and main beam to main beam.
% Main beam = within 3 dB of that antenna's peak gain on the path (NEXRAD at
% the cut's tilt).
pk1 = Rv(Rv.max_g21_dBi >= lev_hi, :);
% Mark every peak above the plateau; letter the three largest and explain
% them in a key (the full list is in peaks_above_plateau_120h_rev5.csv).
plot(ax, pk1.t_start_s, pk1.max_g21_dBi, "o", "Color", c.purple, "MarkerFaceColor", c.purple, "MarkerSize", 3.5);
[~, o3] = sort(pk1.max_g21_dBi, "descend");
o3 = o3(1:min(3, numel(o3)));
for k = 1:numel(o3)
    text(ax, pk1.t_start_s(o3(k)), pk1.max_g21_dBi(o3(k)) + 1.6, char(96 + k), "Color", c.purple, ...
        "FontName", "Times New Roman", "FontSize", 9, "FontWeight", "bold", "HorizontalAlignment", "center");
    txt = peak_label_rev1(pk1(o3(k), :), P, wl.ref_levels.asr11_gain_dBi(1));
    put_label(ax, 20, 93 - 3.4*k, sprintf("%c  %s, t = %.0f s", 96 + k, txt, pk1.t_at_max_s(o3(k))), c.purple, false, "left", "w");
end
style(ax, c, sprintf("Time (s), trial %d", j_fig), "Maximum MAG per revolution (dBi)", "");
ylim(ax, [15 95]); xlim(ax, [0 max(Rv.t_start_s) + 5]);
save_fig(f, out_dir, "fig5_max_per_rev");

% ---- Figure 6: anatomy of the two largest peaks in trial j_fig ---------------
T1 = wf.scan_max.ts1;
[~, o] = sort(pk1.max_g21_dBi, "descend");
f = new_fig(6.5, 3.0);
tl = tiledlayout(f, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for k = 1:min(2, numel(o))
    ax = nexttile(tl);
    t0 = pk1.t_at_max_s(o(k));
    i = abs(T1.t_s - t0) <= 0.4;
    tt = T1.t_s(i) - t0;
    plot(ax, tt, T1.nexrad_gain_dBi(i), "-", "Color", c.gray, "LineWidth", 1.2); hold(ax, "on");
    plot(ax, tt, T1.asr11_gain_dBi(i), "--", "Color", c.purple, "LineWidth", 1.2);
    plot(ax, tt, T1.g21_dBi(i), "-", "Color", c.blue, "LineWidth", 2.0);
    put_label(ax, 0.38, 80, sprintf("MAG, peak %.1f dBi", pk1.max_g21_dBi(o(k))), c.blue, true, "right");
    put_label(ax, 0.38, 72, "NEXRAD gain toward ASR-11", c.gray, false, "right");
    put_label(ax, 0.38, 64, "ASR-11 gain toward NEXRAD", c.purple, true, "right");
    style(ax, c, "Time from peak (s)", "Gain (dBi)", sprintf("%s. Peak at t = %.1f s", char(96 + k), t0));
    ylim(ax, [-15 85]); xlim(ax, [-0.4 0.4]);
end
save_fig(f, out_dir, "fig6_peak_anatomy");

% ---- Figure 7: main beam to main beam versus the statistical MAG -------------
%      MAG against percent of time exceeded, from the CDF arrays
f = new_fig(6.5, 3.2);
ax = axes(f);
Sw = res{i_w}.sensitivity; S215 = res{i_215}.sensitivity; Sp = res{i_p}.sensitivity;
p_fine = 100 - logspace(0, -3, 200)';         % 99% to 99.999%, finer than opts.cdf_pct
cases = cellfun(@(r) [p_fine, prctile(r.cdf.g21_dBi, p_fine)], {res{i_p}, res{i_215}, res{i_w}}, ...
    "UniformOutput", false);
sty = {{"--", c.purple, 1.4}, {"-", c.gray, 1.4}, {"-", c.blue, 2.2}};
for k = 1:3
    A = cases{k};
    A = A(A(:, 1) >= 99 & A(:, 1) < 100, :);
    semilogx(ax, 100 - A(:, 1), A(:, 2), sty{k}{1}, "Color", sty{k}{2}, "LineWidth", sty{k}{3});
    hold(ax, "on");
end
ax.XDir = "reverse";
yline(ax, Sw.mag_dBi(1), "-", "Color", c.dark, "LineWidth", 1.0);
yline(ax, Sw.mag_dBi(2), ":", "Color", c.dark, "LineWidth", 1.0);
put_label(ax, 1, Sw.mag_dBi(1) + 2.5, sprintf("Main beam to main beam, peak gains, %.1f dBi", Sw.mag_dBi(1)), c.dark, false);
put_label(ax, 1, Sw.mag_dBi(2) - 3, sprintf("Main beam to main beam on the path, %.1f dBi", Sw.mag_dBi(2)), c.dark, false);
g01 = Sw.delta_vs_peak_mbmb_dB(Sw.mag_choice == "99.9% time-weighted MAG");
put_label(ax, 0.0011, 71, string(worst), c.blue, true, "right");
put_label(ax, 0.1, 31, sprintf("Worst mode at 0.1%%: %.0f dB below peak main beam to main beam", -g01), c.blue, false);
put_label(ax, 0.01, 65, "Paper baseline", c.purple, true, "right");
put_label(ax, 0.01, 37, "VCP 215", c.gray, false, "right");
style(ax, c, "Time exceeded (%, log scale)", "Mutual antenna gain (dBi)", "");
xlim(ax, [1e-3 1]); ylim(ax, [10 85]);
save_fig(f, out_dir, "fig7_mbmb_vs_statistical");
A3 = {res{i_p}.cdf.cdf_array21, res{i_215}.cdf.cdf_array21, res{i_w}.cdf.cdf_array21};   % opts.cdf_pct grid
cdf_tab = array2table([A3{1}(:, 1), A3{1}(:, 2), A3{2}(:, 2), A3{3}(:, 2)], ...
    "VariableNames", ["pct", "paper_baseline_dBi", "vcp215_dBi", "worst_mode_dBi"]);
writetable(cdf_tab, fullfile(out_dir, "cdf_array_rev5.csv"));
writetable([table(repmat(T.case(i_p), height(Sp), 1), 'VariableNames', "case"), Sp; ...
    table(repmat(T.case(i_215), height(S215), 1), 'VariableNames', "case"), S215; ...
    table(repmat(worst, height(Sw), 1), 'VariableNames', "case"), Sw], ...
    fullfile(out_dir, "mbmb_vs_statistical_rev5.csv"));

% ---- Figure 8: sensitivity, change in MAG at 99.9% and 99.999% -------------
f = new_fig(6.5, 4.2);
ax = axes(f);
Ts = Tsens(2:end, :);
n_s = height(Ts);
yy = (n_s:-1:1)';
for k = 1:n_s
    plot(ax, [0 Ts.d99_9_dB(k)], yy(k) + 0.15*[1 1], "-", "Color", c.blue, "LineWidth", 2.2); hold(ax, "on");
    plot(ax, [0 Ts.d99_999_dB(k)], yy(k) - 0.15*[1 1], "-", "Color", c.gray, "LineWidth", 2.2);
end
plot(ax, Ts.d99_9_dB, yy + 0.15, "o", "Color", c.blue, "MarkerFaceColor", c.blue, "MarkerSize", 4);
plot(ax, Ts.d99_999_dB, yy - 0.15, "o", "Color", c.gray, "MarkerFaceColor", c.gray, "MarkerSize", 4);
xline(ax, 0, "-", "Color", c.dark, "LineWidth", 0.8);
ax.YTick = flip(yy);
ax.YTickLabel = flip(Ts.group + ": " + Ts.label);
xr = max(abs([Ts.d99_9_dB; Ts.d99_999_dB]));
put_label(ax, -xr*0.98, n_s + 0.8, "Upper line of each pair: 99.9%", c.blue, true, "left");
put_label(ax, -xr*0.98, n_s + 0.35, "Lower line: 99.999%", c.gray, false, "left");
style(ax, c, "Change in MAG from the baseline (dB)", "", "");
ylim(ax, [0.4 n_s + 1.1]); xlim(ax, [-xr - 1, xr + 1]);
ax.YGrid = "off";
save_fig(f, out_dir, "fig8_sensitivity");

% ---- Figure 9: ASR-9 and ASR-11, azimuth patterns and exceedance ------------
f = new_fig(6.5, 3.0);
tl = tiledlayout(f, 1, 2, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(tl);
p11 = antenna_prepare_rev1(w11.sys_asr11.tx);
p9 = antenna_prepare_rev1(w9.sys_asr11.tx);
semilogx(ax, x, antenna_gain_rev1(p11, 0, 0, st), "-", "Color", c.gray, "LineWidth", 1.4); hold(ax, "on");
semilogx(ax, x, antenna_gain_rev1(p9, 0, 0, st), "-", "Color", c.blue, "LineWidth", 2);
put_label(ax, 0.012, 12, "ASR-11, M.1464-2 Radar C", c.gray, false);
put_label(ax, 0.012, 7, "ASR-9, M.1464-2 Radar B", c.blue, true);
style(ax, c, "Azimuth offset (deg, log scale)", "Gain (dBi)", "a. Azimuth patterns");
ylim(ax, [-5 40]); xlim(ax, [0.01 180]);
ax = nexttile(tl);
semilogy(ax, w11.cdf.cdf21.gain_dBi, w11.cdf.cdf21.exceed_pct, "-", "Color", c.gray, "LineWidth", 1.4); hold(ax, "on");
semilogy(ax, w9.cdf.cdf21.gain_dBi, w9.cdf.cdf21.exceed_pct, "-", "Color", c.blue, "LineWidth", 2);
semilogy(ax, env9.exceed.gain_dBi, max(env9.exceed.envelope_exceed_pct, 1e-6), "--", "Color", c.dark, "LineWidth", 1.0);
put_label(ax, -20, 3e-4, "Dashed: ASR-9 envelope of all " + height(T9) + " VCP modes", c.dark, false);
put_label(ax, 78, 30, "ASR-9 (blue), ASR-11 (gray)", c.dark, false, "right");
put_label(ax, 78, 10, string(cmp9_label), c.dark, false, "right");
style(ax, c, "Mutual antenna gain (dBi)", "Time exceeded (%, log scale)", "b. Exceedance");
ylim(ax, [1e-4 100]); xlim(ax, [-25 82]);
ax.YGrid = "on";
save_fig(f, out_dir, "fig9_asr9_comparison");

% ---- Tables typed in the report (Steps 2, 3 and 5) -------------------------
vl = [35 31 32 34 215 212 12 112];
vt = cell(numel(vl), 1);
for k = 1:numel(vl)
    sv = nexrad_vcp_rev1(vl(k));
    vt{k} = {vl(k), sv.use, height(sv.cuts), sv.volume_s, ...
        100*sum(sv.cuts.dur_s(sv.cuts.el_deg == 0.5))/sv.volume_s, sv.limits.sails, string(sv.pulse_width)};
end
writetable(cell2table(vertcat(vt{:}), "VariableNames", ["vcp", "use", "cuts", "volume_s", ...
    "t_low_pct", "max_sails", "pulse_width"]), fullfile(out_dir, "vcp_table_rev5.csv"));
sw = res{i_w}.sys_nexrad.scan;
[Gw, elw] = findgroups(sw.cuts.el_deg);
rates = splitapply(@(r) strjoin(string(unique(r)'), " "), sw.cuts.rate_deg_s, Gw);
tw = splitapply(@sum, sw.cuts.dur_s, Gw);
writetable(table(elw, splitapply(@numel, sw.cuts.el_deg, Gw), rates, tw, 100*tw/sw.volume_s, ...
    'VariableNames', ["el_deg", "cuts_per_volume", "rates_deg_s", "time_s", "share_pct"]), ...
    fullfile(out_dir, "worst_mode_weights_rev5.csv"));
pn = antenna_prepare_rev1(res{i_w}.sys_nexrad.tx);
el_all = unique(nexrad_vcp_rev1(212).tilt_deg);
writetable(table(el_all, antenna_gain_rev1(pn, 0, 0, state_az(zeros(size(el_all)), el_all)), ...
    'VariableNames', ["el_deg", "nexrad_gain_on_path_dBi"]), fullfile(out_dir, "nexrad_gain_by_cut_rev5.csv"));
pa = antenna_prepare_rev1(res{i_w}.sys_asr11.tx);
off = [0; 0.725; 2.9; 10; 90; 180];
writetable(table(off, antenna_gain_rev1(pa, 0, 0, state_az(off, 0)), ...
    'VariableNames', ["az_off_deg", "asr11_gain_dBi"]), fullfile(out_dir, "asr11_gain_cut_rev5.csv"));

% ---- Figure A-1: call tree of the repository functions ---------------------
%      Direct calls only (requiredFilesAndProducts "toponly"); tests and the
%      small-run helper merge_overrides_rev1 left out for clarity.
src = fileparts(mfilename("fullpath"));
fl = string({dir(fullfile(src, "*.m")).name});
fl = fl(~startsWith(fl, "test_") & fl ~= "merge_overrides_rev1.m");
stem = erase(fl, ".m");
e1 = strings(0, 1);
e2 = strings(0, 1);
for k = 1:numel(fl)
    dep = matlab.codetools.requiredFilesAndProducts(fullfile(src, fl(k)), "toponly");
    [~, dn] = fileparts(string(dep));
    dn = unique(dn(ismember(dn, stem) & dn ~= stem(k)));
    e1 = [e1; repmat(stem(k), numel(dn), 1)]; %#ok<AGROW>
    e2 = [e2; dn(:)]; %#ok<AGROW>
end
G = digraph(e1, e2);
f = new_fig(6.5, 4.6);
ax = axes(f);
h = plot(ax, G, "Layout", "layered", "Direction", "down", "NodeColor", c.blue, "EdgeColor", c.gray, ...
    "MarkerSize", 5, "ArrowSize", 7, "LineWidth", 0.9, "NodeFontName", "Times New Roman", "NodeFontSize", 7, ...
    "Interpreter", "none");
h.NodeLabel = replace(G.Nodes.Name, "_rev1", "");
axis(ax, "off");
save_fig(f, out_dir, "fig_a1_call_tree");
writetable(table(e1, e2, 'VariableNames', ["caller", "callee"]), fullfile(out_dir, "call_tree_rev5.csv"));

fprintf("Figures and tables written to %s\n", out_dir);

% =========================================================================
function d = report_dir(ov)
% Output folder: docs/figures next to this script, or ov.out_dir (tests).
d = fullfile(fileparts(mfilename("fullpath")), "docs", "figures");
if isfield(ov, "out_dir")
    d = ov.out_dir;
end
if ~isfolder(d)
    mkdir(d);
end
end

function st = state_az(az, tilt)
z = zeros(numel(az), 1);
st = struct("mech_az", az(:), "tilt", tilt + z, "steer_az", z, "steer_el", z);
end

function [tx, ey] = steps(scan)
tx = [0; cumsum(scan.cuts.dur_s)];
ey = [scan.cuts.el_deg; scan.cuts.el_deg(end)];
end

function f = new_fig(w_in, h_in)
f = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(f, "light");              % R2025a+ desktop dark theme; R2024b figures are light already
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

function put_label(ax, x, y, s, col, bold, halign, bg)
if nargin < 7
    halign = "left";
end
if nargin < 8
    bg = "none";
end
w = "normal";
if bold
    w = "bold";
end
text(ax, x, y, s, "Color", col, "FontWeight", w, "FontName", "Times New Roman", ...
    "FontSize", 8, "HorizontalAlignment", halign, "BackgroundColor", bg, "Margin", 0.5);
end

function save_fig(f, d, name)
exportgraphics(f, fullfile(d, name + ".png"), "Resolution", 300);
close(f);
end
