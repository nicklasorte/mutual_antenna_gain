%REPORT_FIGURES_REV4 Figures and tables for methodology report rev4.
%   Runs example_nexrad_asr11_rev1 (WSR-88D vs ASR-11, all VCP modes) and
%   writes to docs/figures:
%     fig1_antenna_patterns.png   NEXRAD and ASR-11 antenna patterns
%     fig2_vcp_schedule.png       elevation versus time in one volume
%     fig3_low_tilt_driver.png    99.9% MAG versus time at 0.5 deg
%     fig4_exceedance.png         MAG exceedance with reference levels
%     fig5_max_per_rev.png        maximum MAG per ASR-11 revolution
%     fig6_mbmb_vs_statistical.png  main beam to main beam versus statistical MAG
%     vcp_sweep_rev4.csv          sweep table (all modes)
%     ref_levels_rev4.csv, by_elevation_rev4.csv, asr_rev_pct_rev4.csv,
%     sensitivity_rev4.csv, mode_envelope_rev4.csv
%   Style: NTIA palette (mid gray context, light blue highlight, purple
%   second highlight, dark blue text), Times New Roman, 6.5 in wide.
%   Titles, captions and sources are in the report text, not the images.

example_nexrad_asr11_rev1;          % defines T, V, res, w, i999, i99
close all;

out_dir = fullfile(fileparts(mfilename("fullpath")), "docs", "figures");
if ~isfolder(out_dir)
    mkdir(out_dir);
end
writetable(T, fullfile(out_dir, "vcp_sweep_rev4.csv"));
writetable(w.ref_levels, fullfile(out_dir, "ref_levels_rev4.csv"));
writetable(w.scan_max.by_elevation, fullfile(out_dir, "by_elevation_rev4.csv"));
writetable(w.scan_max.asr_rev_pct, fullfile(out_dir, "asr_rev_pct_rev4.csv"));
writetable(w.sensitivity, fullfile(out_dir, "sensitivity_rev4.csv"));
writetable(env.pct, fullfile(out_dir, "mode_envelope_rev4.csv"));

c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "purple", [110 79 146]/255, "gray", [138 141 145]/255, "grid", [230 230 230]/255);
i_w = 1 + i999;                     % worst mode (99.9%) in res and T
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

% ---- Figure 3: 99.9% MAG versus time at 0.5 deg ----------------------------
f = new_fig(6.5, 3.0);
ax = axes(f);
scatter(ax, V.t_low_pct, V.p99_9, 28, c.gray, "filled"); hold(ax, "on");
scatter(ax, V.t_low_pct(i999), V.p99_9(i999), 60, c.blue, "filled");
put_label(ax, V.t_low_pct(i999) - 1, V.p99_9(i999) + 0.4, string(worst) + " (worst)", c.blue, true, "right");
[~, i_lo] = min(V.p99_9);
put_label(ax, V.t_low_pct(i_lo) + 0.6, V.p99_9(i_lo) + 0.2, string(V.case(i_lo)), c.dark, false);
yline(ax, T.p99_9(i_p), "--", "Color", c.purple, "LineWidth", 1.2);
put_label(ax, 13.5, T.p99_9(i_p) + 0.35, sprintf("Paper baseline, %.1f dBi", T.p99_9(i_p)), c.purple, true);
style(ax, c, "NEXRAD time at 0.5 deg (% of volume)", "99.9% mutual gain (dBi)", "");
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

% ---- Figure 5: maximum MAG per ASR-11 revolution ---------------------------
f = new_fig(6.5, 3.0);
ax = axes(f);
Rv = w.scan_max.asr_rev(w.scan_max.asr_rev.trial == 1 & w.scan_max.asr_rev.complete, :);
hi = Rv.max_g21_dBi >= R.mag_dBi(3) - 0.2;
scatter(ax, Rv.t_start_s(~hi), Rv.max_g21_dBi(~hi), 8, c.gray, "filled"); hold(ax, "on");
scatter(ax, Rv.t_start_s(hi), Rv.max_g21_dBi(hi), 22, c.blue, "filled");
yline(ax, R.mag_dBi(5), ":", "Color", c.dark);
yline(ax, R.mag_dBi(1), ":", "Color", c.dark);
put_label(ax, 20, R.mag_dBi(5) - 3, sprintf("ASR-11 main beam + NEXRAD back lobe, %.1f dBi", R.mag_dBi(5)), c.dark, false);
put_label(ax, 20, R.mag_dBi(1) + 3, sprintf("Main beam to main beam, %.1f dBi", R.mag_dBi(1)), c.dark, false);
put_label(ax, 1900, 62, sprintf("%.1f%% of revolutions reach %.1f dBi or more", ...
    R.asr_rev_reached_pct(3), R.mag_dBi(3)), c.blue, true);
style(ax, c, "Time (s), trial 1", "Maximum mutual gain per revolution (dBi)", "");
ylim(ax, [15 85]); xlim(ax, [0 max(Rv.t_start_s) + 5]);
save_fig(f, out_dir, "fig5_max_per_rev");

% ---- Figure 6: main beam to main beam versus the statistical MAG -------------
f = new_fig(6.5, 3.2);
ax = axes(f);
Sw = res{i_w}.sensitivity; S215 = res{i_215}.sensitivity; Sp = res{i_p}.sensitivity;
keep = [1 2 7 6 5 4];                          % peak MB-MB, path MB-MB, 99.999 ... 99%
yl = ["Main beam to main beam (peak gains)", "Main beam to main beam (on path)", ...
    "99.999% MAG", "99.99% MAG", "99.9% MAG", "99% MAG"];
y = numel(keep):-1:1;
plot(ax, Sp.mag_dBi(keep), y + 0.18, "d", "MarkerFaceColor", c.purple, "MarkerEdgeColor", "none", "MarkerSize", 6);
hold(ax, "on");
plot(ax, S215.mag_dBi(keep), y - 0.18, "o", "MarkerFaceColor", c.gray, "MarkerEdgeColor", "none", "MarkerSize", 6);
plot(ax, Sw.mag_dBi(keep), y, "o", "MarkerFaceColor", c.blue, "MarkerEdgeColor", "none", "MarkerSize", 8);
xline(ax, Sw.mag_dBi(1), "-", "Color", c.dark, "LineWidth", 0.8);
for k = 3:numel(keep)
    text(ax, Sw.mag_dBi(keep(k)) + 1.5, y(k), sprintf("%.0f dB below", -Sw.delta_vs_peak_mbmb_dB(keep(k))), ...
        "Color", c.blue, "FontWeight", "bold", "FontName", "Times New Roman", "FontSize", 8);
end
put_label(ax, 2, 6.75, string(worst), c.blue, true);
put_label(ax, 2, 6.35, "VCP 215", c.gray, false);
put_label(ax, 2, 5.95, "Paper baseline (F.1245, 3 rpm)", c.purple, true);
ax.YTick = 1:numel(keep);
ax.YTickLabel = fliplr(yl);
style(ax, c, "Mutual antenna gain (dBi)", "", "");
ax.YGrid = "off"; ax.XGrid = "on";
xlim(ax, [0 85]); ylim(ax, [0.4 numel(keep) + 0.9]);
save_fig(f, out_dir, "fig6_mbmb_vs_statistical");
writetable([table(repmat(T.case(i_p), height(Sp), 1), 'VariableNames', "case"), Sp; ...
    table(repmat(T.case(i_215), height(S215), 1), 'VariableNames', "case"), S215; ...
    table(repmat(worst, height(Sw), 1), 'VariableNames', "case"), Sw], ...
    fullfile(out_dir, "mbmb_vs_statistical_rev4.csv"));

fprintf("Figures and tables written to %s\n", out_dir);

% =========================================================================
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
