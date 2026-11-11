function report_figures_rev11(out_dir)
%REPORT_FIGURES_REV11 Draws the periodic-lock figure of revision 11 from the CSV files.
%   report_figures_rev11(out_dir) reads lock_trials_rev11.csv in out_dir,
%   written by report_numbers_rev11, and saves fig_lock_periodic_rev11.png:
%   the longest periodic run (high volumes recurring every k = 1 to 4
%   volumes) in each fixed-offset trial against m4. Runs cut by a trial end
%   are drawn hollow: their length is a lower bound.
%
%   DRAFT//PRE-DECISIONAL

arguments
    out_dir (1,1) string = fullfile(fileparts(mfilename("fullpath")), "docs", "figures")
end

c = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, ...
    "gray", [138 141 145]/255, "grid", [0.9 0.9 0.9]);
T_lock = readtable(fullfile(out_dir, "lock_trials_rev11.csv"), "TextType", "string");
modes = unique(T_lock.mode, "stable");
f = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(f, "light");
end
f.Position(3:4) = [6.5 3.0];
ax = axes(f);
hold(ax, "on");
mk = ["o", "s"];
cols = {c.blue, c.gray};
for m = 1:min(2, numel(modes))
    k = T_lock.mode == modes(m) & T_lock.el_deg == 0.5 & T_lock.threshold_dBi == 60 & T_lock.per_vols > 0;
    L = T_lock(k, :);
    x = max(L.m4_deg, 0.01);                        % keep exact locks on the log axis
    y = L.per_vols;
    cen = L.per_censored == 1;                      % read back from CSV as 0 or 1
    plot(ax, x(~cen), y(~cen), mk(m), "Color", cols{m}, "MarkerFaceColor", cols{m}, "MarkerSize", 3.5);
    plot(ax, x(cen), y(cen), mk(m), "Color", cols{m}, "MarkerFaceColor", "w", "MarkerSize", 5, "LineWidth", 1.1);
    fw = "normal";
    if m == 1, fw = "bold"; end
    text(ax, 0.012, 170/(1.5^(m - 1)), modes(m) + " (" + mk(m) + ")", "Color", cols{m}, ...
        "FontWeight", fw, "FontName", "Times New Roman", "FontSize", 8);
end
ax.XScale = "log";
ax.YScale = "log";
xlim(ax, [0.01 200]);
ylim(ax, [1 300]);
box(ax, "off");
ax.FontName = "Times New Roman";
ax.FontSize = 9;
ax.TickDir = "out";
ax.XColor = c.dark;
ax.YColor = c.dark;
ax.YGrid = "on";
ax.GridColor = c.grid;
ax.GridAlpha = 1;
xlabel(ax, "Smallest drift from whole ATC revolutions over 1 to 4 volumes, m4 (degrees)");
ylabel(ax, "Longest periodic run in the trial (volumes)");
text(ax, 0.012, 250, "Hollow: cut by the trial end (a lower bound)", "Color", c.dark, ...
    "FontName", "Times New Roman", "FontSize", 8);
exportgraphics(f, fullfile(out_dir, "fig_lock_periodic_rev11.png"), "Resolution", 300);
close(f);
end
