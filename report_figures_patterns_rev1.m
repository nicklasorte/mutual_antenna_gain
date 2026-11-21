%REPORT_FIGURES_PATTERNS_REV1 Figures and tables for the antenna patterns working paper.
%   Runs the worst-case horizontal slice functions for the WSR-88D, ASR-11
%   and ASR-9, compares them with the generic statgain (NTIA TM-13-489)
%   and Rec. ITU-R F.1245 patterns, and writes to docs/figures:
%
%     fig_patterns_slices_rev1.png     0-360 deg worst-case slices, one panel per radar
%     fig_patterns_zoom_rev1.png       the main beams, -10 to +10 deg
%     fig_patterns_elevation_rev1.png  ASR elevation pattern; horizon gain versus tilt
%     fig_patterns_compare_rev1.png    each radar against statgain and F.1245
%     nexrad_worst_case_slice_rev1.csv, asr11_worst_case_slice_rev1.csv,
%     asr9_worst_case_slice_rev1.csv   the N-by-2 slices at 0.1 deg
%     patterns_compare_rev1.csv        gain at selected angles, every pattern
%     patterns_params_rev1.csv         every parameter with its source
%
%   Run from the repository root. Set out_dir before running to write
%   elsewhere.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also nexrad_worst_case_slice_rev1, asr11_worst_case_slice_rev1,
%   asr9_worst_case_slice_rev1, statgain_worst_case_slice_rev1.

if ~exist("out_dir", "var")
    out_dir = fullfile("docs", "figures");
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---- Slices (0.01 deg for the figures, 0.1 deg for the CSV files)
[s_nex, i_nex] = nexrad_worst_case_slice_rev1(step_deg=0.01);
[s_nssl, i_nssl] = nexrad_worst_case_slice_rev1(step_deg=0.01, pattern="nssl_measured");
[s_a11, i_a11] = asr11_worst_case_slice_rev1(step_deg=0.01);
[s_a9, i_a9] = asr9_worst_case_slice_rev1(step_deg=0.01);
names = ["nexrad", "asr11", "asr9"];
fns = {@nexrad_worst_case_slice_rev1, @asr11_worst_case_slice_rev1, @asr9_worst_case_slice_rev1};
for k = 1:3
    s = fns{k}(step_deg=0.1);
    writetable(array2table(s, 'VariableNames', ["azimuth_deg", "gain_dBi"]), ...
        fullfile(out_dir, names(k) + "_worst_case_slice_rev1.csv"));
end

% ---- Generic patterns at each radar's peak gain
[s_sg_nex, i_sg_nex] = statgain_worst_case_slice_rev1(i_nex.g_max_dBi, step_deg=0.01);
[s_sg_a11, i_sg_a11] = statgain_worst_case_slice_rev1(i_a11.g_max_dBi, step_deg=0.01);
[s_sg_a9, i_sg_a9] = statgain_worst_case_slice_rev1(i_a9.g_max_dBi, step_deg=0.01);
% F.1245 as the 2.7 GHz link budget working paper uses it (WP_REV3 Table
% 4.4-1: 46 dBi WSR-88D, 34 dBi ASR), from the repository antenna model.
f1245 = @(g, az) antenna_gain_rev1(antenna_prepare_rev1(struct("type", "f1245", "g_max_dBi", g)), ...
    0, 0, struct("mech_az", az, "tilt", zeros(size(az))));
az = s_nex(:, 1);
g_f_nex = f1245(46, az);
g_f_asr = f1245(34, az);

% ---- Comparison table at selected angles
ang = [0 0.25 0.5 0.75 1 1.5 2 3 5 10 20 45 90 180]';
at = @(s) interp1(s(:, 1), s(:, 2), ang);
T = table(ang, at(s_nex), at(s_nssl), at(s_sg_nex), interp1(az, g_f_nex, ang), ...
    at(s_a11), at(s_a9), at(s_sg_a11), interp1(az, g_f_asr, ang), ...
    'VariableNames', ["azimuth_deg", "wsr88d_roc_spec_dBi", "wsr88d_nssl_measured_dBi", ...
    "statgain_45p5_dBi", "f1245_46_dBi", "asr11_dBi", "asr9_dBi", "statgain_34_dBi", "f1245_34_dBi"]);
writetable(T, fullfile(out_dir, "patterns_compare_rev1.csv"));
P = [addvars(i_nex.params, repmat("WSR-88D", height(i_nex.params), 1), 'Before', 1, 'NewVariableNames', "radar"); ...
    addvars(i_a11.params, repmat("ASR-11", height(i_a11.params), 1), 'Before', 1, 'NewVariableNames', "radar"); ...
    addvars(i_a9.params, repmat("ASR-9", height(i_a9.params), 1), 'Before', 1, 'NewVariableNames', "radar")];
writetable(P, fullfile(out_dir, "patterns_params_rev1.csv"));

% ---- Style (NTIA palette)
C = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, "gray", [138 141 145]/255, ...
    "purple", [112 48 160]/255, "grid", [0.9 0.9 0.9]);
fs = {"FontName", "Times New Roman", "FontSize", 8};
deg = char(176);

% ---- Figure 1: 0-360 deg slices
fh = new_fig(6.5, 6.2);
tl = tiledlayout(fh, 3, 1, "TileSpacing", "compact", "Padding", "compact");
S = {s_nex, s_a11, s_a9};
I = {i_nex, i_a11, i_a9};
for k = 1:3
    ax = nexttile(tl);
    plot(ax, S{k}(:, 1), S{k}(:, 2), "Color", C.blue, "LineWidth", 1.1);
    style_ax(ax, C, fs);
    xlim(ax, [0 360]); xticks(ax, 0:45:360); ylim(ax, [-5 50]);
    ylabel(ax, "Gain (dBi)", fs{:});
    fl = min(S{k}(:, 2));
    yline(ax, fl, ":", "Color", C.gray, "LineWidth", 0.8);
    text(ax, 180, 44, sprintf("%s: peak %.1f dBi at 0%s and 360%s; floor %.1f dBi", ...
        I{k}.name, I{k}.peak_dBi, deg, deg, fl), fs{:}, "Color", C.dark, "HorizontalAlignment", "center", ...
        "FontWeight", "bold");
    if k == 3
        xlabel(ax, "Beam azimuth relative to the bearing of the station (degrees)", fs{:});
    end
end
exportgraphics(fh, fullfile(out_dir, "fig_patterns_slices_rev1.png"), "Resolution", 300); close(fh);

% ---- Figure 2: main beams, -10 to +10 deg
fh = new_fig(6.5, 3.4);
ax = axes(fh);
sg = @(s) [s(s(:, 1) > 180, 1) - 360; s(s(:, 1) <= 180, 1)];
gv = @(s) [s(s(:, 1) > 180, 2); s(s(:, 1) <= 180, 2)];
plot(ax, sg(s_nex), gv(s_nex), "Color", C.blue, "LineWidth", 1.3); hold(ax, "on");
plot(ax, sg(s_a11), gv(s_a11), "Color", C.dark, "LineWidth", 1.1);
plot(ax, sg(s_a9), gv(s_a9), "--", "Color", C.gray, "LineWidth", 1.1);
yline(ax, i_a11.first_sl_dBi, ":", "Color", C.purple, "LineWidth", 0.8);
style_ax(ax, C, fs);
xlim(ax, [-10 10]); ylim(ax, [0 50]);
xlabel(ax, "Beam azimuth relative to the bearing of the station (degrees)", fs{:});
ylabel(ax, "Gain (dBi)", fs{:});
text(ax, 1.0, 44, "WSR-88D (ROC specification)", fs{:}, "Color", C.blue, "FontWeight", "bold");
text(ax, 5, 1.8, "ASR-11 (floor 3.5 dBi)", fs{:}, "Color", C.dark, "FontWeight", "bold");
text(ax, 5, 7.3, "ASR-9 (floor 5.5 dBi)", fs{:}, "Color", C.gray, "FontWeight", "bold");
text(ax, -9.7, i_a11.first_sl_dBi + 1.3, "ASR-11 first sidelobe in M.1464-2 (+9.5 dBi)", fs{:}, "Color", C.purple, "FontWeight", "bold");
exportgraphics(fh, fullfile(out_dir, "fig_patterns_zoom_rev1.png"), "Resolution", 300); close(fh);

% ---- Figure 3: ASR elevation pattern; horizon gain versus tilt
fh = new_fig(6.5, 3.2);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(tl);
el = (-6:0.01:40)';
[~, i_t] = asr11_worst_case_slice_rev1(tilt_deg=2.5);
plot(ax, el, i_a11.g_max_dBi + i_a11.el_fn(el), "Color", C.blue, "LineWidth", 1.2); hold(ax, "on");
plot(ax, el, i_a11.g_max_dBi + i_t.el_fn(el), "--", "Color", C.gray, "LineWidth", 1.0);
xline(ax, 0, ":", "Color", C.dark);
style_ax(ax, C, fs); xlim(ax, [-6 40]); ylim(ax, [0 37]);
xlabel(ax, "Elevation (degrees)", fs{:}); ylabel(ax, "ASR-11 gain at azimuth 0 (dBi)", fs{:});
text(ax, 6, 33, "Tilt 0" + deg + ": peak on the horizon", fs{:}, "Color", C.blue, "FontWeight", "bold");
text(ax, 17, 24, "Tilt 2.5" + deg, fs{:}, "Color", C.gray, "FontWeight", "bold");
text(ax, 0.5, 2, "Horizon", fs{:}, "Color", C.dark);
title(ax, "a. Cosecant-squared elevation pattern", fs{:}, "FontWeight", "normal", "Color", C.dark);
ax = nexttile(tl);
tilt = (0:0.05:2.5)';
h_nex = arrayfun(@(t) peak_at(@nexrad_worst_case_slice_rev1, t), tilt);
h_a11 = arrayfun(@(t) peak_at(@asr11_worst_case_slice_rev1, t), tilt);
h_a9 = arrayfun(@(t) peak_at(@asr9_worst_case_slice_rev1, t), tilt);
plot(ax, tilt, h_nex, "Color", C.blue, "LineWidth", 1.2); hold(ax, "on");
plot(ax, tilt, h_a11, "Color", C.dark, "LineWidth", 1.1);
plot(ax, tilt, h_a9, "--", "Color", C.gray, "LineWidth", 1.1);
plot(ax, 0.5, peak_at(@nexrad_worst_case_slice_rev1, 0.5), "o", "Color", C.blue, "MarkerFaceColor", C.blue, "MarkerSize", 4);
style_ax(ax, C, fs); xlim(ax, [0 2.5]); ylim(ax, [0 50]);
xlabel(ax, "Beam elevation (degrees)", fs{:}); ylabel(ax, "Peak gain toward the horizon (dBi)", fs{:});
text(ax, 0.6, 43, sprintf("WSR-88D 0.5%s cut: %.1f dBi", deg, peak_at(@nexrad_worst_case_slice_rev1, 0.5)), ...
    fs{:}, "Color", C.blue, "FontWeight", "bold");
text(ax, 1.2, 36.5, "ASR-11", fs{:}, "Color", C.dark, "FontWeight", "bold");
text(ax, 1.2, 28.5, "ASR-9", fs{:}, "Color", C.gray, "FontWeight", "bold");
title(ax, "b. Highest at 0" + deg + " tilt for every radar", fs{:}, "FontWeight", "normal", "Color", C.dark);
exportgraphics(fh, fullfile(out_dir, "fig_patterns_elevation_rev1.png"), "Resolution", 300); close(fh);

% ---- Figure 4: each radar against statgain and F.1245 (log angle)
fh = new_fig(6.5, 7.0);
tl = tiledlayout(fh, 3, 1, "TileSpacing", "compact", "Padding", "compact");
half = az <= 180 & az > 0;
rows = {
    "WSR-88D", s_nex, s_sg_nex, g_f_nex, "46 dBi", i_sg_nex
    "ASR-11", s_a11, s_sg_a11, g_f_asr, "34 dBi", i_sg_a11
    "ASR-9", s_a9, s_sg_a9, g_f_asr, "34 dBi", i_sg_a9};
for k = 1:3
    ax = nexttile(tl);
    s = rows{k, 2};
    semilogx(ax, az(half), s(half, 2), "Color", C.blue, "LineWidth", 1.4); hold(ax, "on");
    if k == 1
        semilogx(ax, az(half), s_nssl(half, 2), "-.", "Color", C.blue, "LineWidth", 0.9);
    end
    semilogx(ax, az(half), rows{k, 3}(half, 2), "Color", C.gray, "LineWidth", 1.1);
    semilogx(ax, az(half), rows{k, 4}(half), "--", "Color", C.purple, "LineWidth", 1.0);
    style_ax(ax, C, fs);
    xlim(ax, [0.05 180]); ylim(ax, [-15 50]);
    xticks(ax, [0.1 0.3 1 3 10 30 100 180]); xticklabels(ax, ["0.1" "0.3" "1" "3" "10" "30" "100" "180"]);
    ylabel(ax, "Gain (dBi)", fs{:});
    title(ax, sprintf("%s: radar-specific pattern (blue), statgain at %.4g dBi (gray), F.1245 at %s (purple, dashed)", ...
        rows{k, 1}, rows{k, 6}.g_max_dBi, rows{k, 5}), fs{:}, "FontWeight", "normal", "Color", C.dark);
    if k == 1
        text(ax, 25, 8.5, "ROC specification", fs{:}, "Color", C.blue, "FontWeight", "bold");
        text(ax, 25, -4.5, "NSSL measured", fs{:}, "Color", C.blue);
    end
    if k == 3
        xlabel(ax, "Beam azimuth from the station (degrees, log scale)", fs{:});
    end
end
exportgraphics(fh, fullfile(out_dir, "fig_patterns_compare_rev1.png"), "Resolution", 300); close(fh);

disp(T);
fprintf("WSR-88D: x = %.4f deg; far field %.0f m\n", i_nex.x_deg, i_nex.far_field_m);
fprintf("ASR-11: brk %.4f deg, floor at %.4f deg; ASR-9: brk %.4f deg, floor at %.4f deg\n", ...
    i_a11.brk_deg, i_a11.floor_az_deg, i_a9.brk_deg, i_a9.floor_az_deg);
fprintf("Statgain implied beamwidths: %.3f (45.5), %.3f (34), %.3f (33.5) deg\n", ...
    i_sg_nex.bw_implied_deg, i_sg_a11.bw_implied_deg, i_sg_a9.bw_implied_deg);
fprintf("Share of the revolution above the floor: WSR-88D %.2f%%, ASR-11 %.2f%%, ASR-9 %.2f%%\n", ...
    100*mean(s_nex(1:end-1, 2) > min(s_nex(:, 2)) + 1e-9), 100*mean(s_a11(1:end-1, 2) > min(s_a11(:, 2)) + 1e-9), ...
    100*mean(s_a9(1:end-1, 2) > min(s_a9(:, 2)) + 1e-9));

% -------------------------------------------------------------------------
function g = peak_at(fn, tilt)
s = fn(step_deg=0.5, tilt_deg=tilt);
g = s(1, 2);
end

function fh = new_fig(w, h)
fh = figure("Visible", "off", "Units", "inches", "Position", [1 1 w h], "Color", "w");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(fh, "light");
end
end

function style_ax(ax, C, fs)
set(ax, fs{:}, "Box", "off", "TickDir", "out", "XColor", C.dark, "YColor", C.dark, ...
    "GridColor", C.grid, "GridAlpha", 1, "YGrid", "on", "Layer", "top");
hold(ax, "on");
end
