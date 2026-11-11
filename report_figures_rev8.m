function report_figures_rev8(out_dir)
%REPORT_FIGURES_REV8 Draws Figures 12 to 14 of revision 8 from the CSV files.
%   report_figures_rev8(out_dir) reads vol_exceed_modes_rev8.csv,
%   seeds_runs_rev8.csv, and ts_all_modes_rev8.csv in out_dir, written by
%   report_numbers_rev8, with ts_envelope_rev8.csv for the persistence-worst
%   mode of Figure 13, and saves fig12_volume_exceedance.png,
%   fig13_runs_vs_drift.png, and fig14_harm.png there. Reading the tables
%   instead of the workspace lets the figures be redrawn without the
%   long runs.
%
%   DRAFT//PRE-DECISIONAL

arguments
    out_dir (1,1) string = fullfile(fileparts(mfilename("fullpath")), "docs", "figures")
end

c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "gray", [138 141 145]/255, "grid", [0.9 0.9 0.9]);
T_ex = readtable(fullfile(out_dir, "vol_exceed_modes_rev8.csv"));
T_seed = readtable(fullfile(out_dir, "seeds_runs_rev8.csv"), "TextType", "string");
T_all = readtable(fullfile(out_dir, "ts_all_modes_rev8.csv"), "TextType", "string");
T_env = readtable(fullfile(out_dir, "ts_envelope_rev8.csv"), "TextType", "string");

% Figure 12: share of volumes reaching each level, envelope over the 216
% modes and the CDF-worst mode, path at 0.5 deg.
f = new_fig(6.5, 2.8);
ax = axes(f);
hold(ax, "on");
plot(ax, T_ex.gain_dBi, T_ex.cdf_worst_05_pct, "Color", c.gray, "LineWidth", 1.5);
plot(ax, T_ex.gain_dBi, T_ex.env_05_pct, "Color", c.blue, "LineWidth", 2);
xlim(ax, [45 80]);
ylim(ax, [0 100]);
style(ax, c, "Largest MAG in the volume (dBi)", "Volumes reaching it (%)");
put_label(ax, 62.5, interp1(T_ex.gain_dBi, T_ex.env_05_pct, 62) + 7, ...
    "Largest over all 216 modes", c.blue, true);
put_label(ax, 52, 12, "VCP 212 SAILS x3 AVSET 6.4 (worst mode for the CDF)", c.gray, false);
save_fig(f, out_dir, "fig12_volume_exceedance");

% Figure 13: longest run at 60 dBi over 20 seeds, by drift model, path at
% 0.5 deg, persistence-worst mode.
mode_p = T_env.worst_mode(T_env.el_deg == 0.5 & T_env.metric == "vol_pct_60");   % persistence-worst mode
models = ["per_trial", "random_walk", "random_walk"];
steps = [NaN 0.01 0.05];
lbl = ["Fixed offset", "Walk 0.01", "Walk 0.05"];
rst = RandStream("mrg32k3a", "Seed", 3);            % jitter of the dots only
f = new_fig(6.5, 2.8);
ax = axes(f);
hold(ax, "on");
y_top = 0;
for d = 1:numel(models)
    sel = T_seed.mode == mode_p & T_seed.drift_model == models(d) & ...
        (isnan(steps(d)) & isnan(T_seed.step_deg_s_per_rev) | T_seed.step_deg_s_per_rev == steps(d)) & ...
        T_seed.el_deg == 0.5 & T_seed.threshold_dBi == 60;
    if ~any(sel)
        error("report_figures:noSeeds", "No seed rows for %s, %s.", mode_p, lbl(d));
    end
    y = T_seed.max_run_vols(sel);
    y_top = max(y_top, max(y));
    x = d + 0.25*(rand(rst, numel(y), 1) - 0.5);
    plot(ax, x, y, "o", "Color", c.blue, "MarkerFaceColor", c.blue, "MarkerSize", 4);
    plot(ax, d + [-0.2 0.2], median(T_seed.p95_run_vols(sel))*[1 1], "-", "Color", c.gray, "LineWidth", 2);
end
ax.XTick = 1:numel(models);
ax.XTickLabel = lbl;
xlim(ax, [0.5 numel(models) + 0.5]);
ylim(ax, [0 10*ceil((y_top + 1)/10)]);
style(ax, c, "ATC rate drift model (random-walk step in deg/s per revolution)", "Run of volumes at 60 dBi");
yl = ylim(ax);
put_label(ax, 1.6, 0.92*yl(2), "Dots: longest run in each of 20 seeds", c.blue, true);
put_label(ax, 1.6, 0.80*yl(2), "Bars: median of the 95th-percentile runs", c.gray, false);
save_fig(f, out_dir, "fig13_runs_vs_drift");

% Figure 14: share of 0.5 deg cuts with a radial at 60 dBi for each mode,
% sorted, path at 0 and 0.5 deg.
f = new_fig(6.5, 2.6);
ax = axes(f);
hold(ax, "on");
n_m = 0;
for a = [0 0.5]
    v = sort(T_all.cuts_hit_60_pct(T_all.el_deg == a), "descend");
    n_m = max(n_m, numel(v));
    col = c.gray;
    if a == 0.5, col = c.blue; end
    plot(ax, 1:numel(v), v, "-", "Color", col, "LineWidth", 1.8);
end
xlim(ax, [0 n_m]);
style(ax, c, "WSR-88D scan modes, sorted from worst", "0.5° cuts with a radial at 60 dBi (%)");
v5 = T_all.cuts_hit_60_pct(T_all.el_deg == 0.5);
put_label(ax, 8, max(v5)*0.97, "Path at 0.5°", c.blue, true);
put_label(ax, 8, max(T_all.cuts_hit_60_pct(T_all.el_deg == 0))*0.85, "Path at 0°", c.gray, false);
save_fig(f, out_dir, "fig14_harm");
end

function f = new_fig(w_in, h_in)
f = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(f, "light");
end
f.Position(3:4) = [w_in h_in];
end

function style(ax, c, xl, yl)
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
end

function put_label(ax, x, y, s, col, bold)
w = "normal";
if bold
    w = "bold";
end
text(ax, x, y, s, "Color", col, "FontWeight", w, "FontName", "Times New Roman", "FontSize", 8);
end

function save_fig(f, d, name)
exportgraphics(f, fullfile(d, name + ".png"), "Resolution", 300);
close(f);
end
