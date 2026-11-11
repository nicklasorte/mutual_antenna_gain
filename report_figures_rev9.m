function report_figures_rev9(out_dir)
%REPORT_FIGURES_REV9 Draws Figure 13 and the lock figure of revision 9 from the CSV files.
%   report_figures_rev9(out_dir) reads seeds_runs_rev9.csv and
%   lock_trials_rev9.csv in out_dir, written by
%   report_numbers_rev9, and saves fig13_runs_vs_drift_rev9.png and
%   fig_lock_vs_slip_rev9.png there. Runs cut by a trial end are drawn
%   hollow: their length is a lower bound.
%
%   DRAFT//PRE-DECISIONAL

arguments
    out_dir (1,1) string = fullfile(fileparts(mfilename("fullpath")), "docs", "figures")
end

c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "gray", [138 141 145]/255, "grid", [0.9 0.9 0.9]);
T_seed = readtable(fullfile(out_dir, "seeds_runs_rev9.csv"), "TextType", "string");
T_lock = readtable(fullfile(out_dir, "lock_trials_rev9.csv"), "TextType", "string");
modes = unique(T_seed.mode, "stable");
mode_p = modes(1);                                  % persistence-worst mode, first in the run order

% Figure 13: longest run at 60 dBi in each seed, by drift model, path at
% 0.5 deg, persistence-worst mode; hollow when the trial end cut it.
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
    cen = T_seed.max_censored_run_vols(sel) >= y & y > 0;   % the longest run is a censored one
    y_top = max(y_top, max(y));
    x = d + 0.25*(rand(rst, numel(y), 1) - 0.5);
    plot(ax, x(~cen), y(~cen), "o", "Color", c.blue, "MarkerFaceColor", c.blue, "MarkerSize", 4);
    plot(ax, x(cen), y(cen), "o", "Color", c.blue, "MarkerFaceColor", "w", "MarkerSize", 5, "LineWidth", 1.2);
    plot(ax, d + [-0.2 0.2], median(T_seed.p95_run_vols(sel))*[1 1], "-", "Color", c.gray, "LineWidth", 2);
end
ax.XTick = 1:numel(models);
ax.XTickLabel = lbl;
xlim(ax, [0.5 numel(models) + 0.5]);
ylim(ax, [0 10*ceil((y_top + 1)/10)]);
style(ax, c, "ATC rate drift model (random-walk step standard deviation, deg/s per revolution)", ...
    "Run of volumes at 60 dBi");
yl = ylim(ax);
put_label(ax, 1.6, 0.92*yl(2), "Dots: longest run in each of 20 seeds (hollow: cut by the trial end)", c.blue, true);
put_label(ax, 1.6, 0.80*yl(2), "Bars: median of the 95th-percentile runs", c.gray, false);
save_fig(f, out_dir, "fig13_runs_vs_drift_rev9");

% Lock figure: longest run in each fixed-offset trial against the ATC
% phase slip per WSR-88D volume, path at 0.5 deg, 60 dBi, both modes.
f = new_fig(6.5, 3.0);
ax = axes(f);
hold(ax, "on");
mk = ["o", "s"];
cols = {c.blue, c.gray};
for m = 1:min(2, numel(modes))
    k = T_lock.mode == modes(m) & T_lock.el_deg == 0.5 & T_lock.threshold_dBi == 60 & T_lock.longest_run_vols > 0;
    L = T_lock(k, :);
    x = abs(L.slip_deg_per_vol);
    y = L.longest_run_vols;
    cen = L.longest_censored == 1;                  % read back from CSV as 0 or 1
    plot(ax, x(~cen), y(~cen), mk(m), "Color", cols{m}, "MarkerFaceColor", cols{m}, "MarkerSize", 3.5);
    plot(ax, x(cen), y(cen), mk(m), "Color", cols{m}, "MarkerFaceColor", "w", "MarkerSize", 5, "LineWidth", 1.1);
    put_label(ax, 0.15, 160/(1.8^(m - 1)), modes(m) + " (" + mk(m) + ")", cols{m}, m == 1);
end
ax.XScale = "log";
ax.YScale = "log";
xlim(ax, [0.1 200]);
ylim(ax, [1 300]);
style(ax, c, "ATC phase slip per WSR-88D volume (degrees, absolute)", "Longest run in the trial (volumes)");
put_label(ax, 8, 200, "Hollow: cut by the trial end (a lower bound)", c.dark, false);
save_fig(f, out_dir, "fig_lock_vs_slip_rev9");
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
