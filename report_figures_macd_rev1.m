%REPORT_FIGURES_MACD_REV1 Walkthrough figures and tables for the MACD section.
%   Runs macd_rev1, the MATLAB port of the C# MACD calculation, for the
%   WSR-88D and ASR-11 of Section 1.1 and writes the figures and tables of
%   the MACD section of the methodology report:
%
%     fig_macd_setup_rev1.png      the setup: two beams, two pointing angles
%     fig_macd_sampling_rev1.png   Step 1: pointings at a quarter beamwidth
%     fig_macd_hist_rev1.png       Step 2: each antenna's 1 dB gain histogram
%     fig_macd_modes_rev1.png      Step 3: WSR-88D cuts averaged by time share
%     fig_macd_convgrid_rev1.png   Step 4: every pair of bins and the MAG
%                                  histogram
%     fig_macd_lookup_rev1.png     Step 5: reading the 99% and 99.999% values
%     fig_macd_conv_rev1.png       result: MAG exceedance, MACD settings
%                                  against a fine reference
%     macd_pct_rev1.csv            percentiles: case, settings, both lookups
%     macd_resolution_rev1.csv     percentiles against step and bin width
%
%   Cases: "one cut" (WSR-88D at 0.5 deg only, as Section 1.1) and
%   "VCP 212" (VCP 212 SAILSx3 AVSET 6.4, each cut weighted by its share of
%   the volume time). Both on a horizontal path. The patterns are the
%   Section 1.1 equations (simple_mag_nexrad_asr11_rev1).
%
%   Set macd_fig_overrides (struct: out_dir, fine_step_div, fine_bin_dB) for
%   a quick run.
%
%   DRAFT//PRE-DECISIONAL
%
%   See also macd_rev1, simple_mag_nexrad_asr11_rev1.

if ~exist("macd_fig_overrides", "var"), macd_fig_overrides = struct(); end
here = fileparts(which("report_figures_macd_rev1"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "fine_step_div", 400, "fine_bin_dB", 0.001);
f = fieldnames(macd_fig_overrides);
for i = 1:numel(f)
    ov.(f{i}) = macd_fig_overrides.(f{i});
end
if ~isfolder(ov.out_dir), mkdir(ov.out_dir); end

% ---- The two antennas (Section 1.1 equations)
nex_gain = @(az, tilt) 45.5 + max(-12*(acosd(min(cosd(tilt).*cosd(az), 1))/0.925).^2, ...
    -29*(acosd(min(cosd(tilt).*cosd(az), 1)) <= 10) - 40*(acosd(min(cosd(tilt).*cosd(az), 1)) > 10));
asr_gain = @(az, ~) max(34 - 12*(abs(az)/1.45).^2, 9.5*(abs(az) <= 3.55) + 3.5*(abs(az) > 3.55));
nex = struct("name", "WSR-88D", "g_max_dBi", 45.5, "bw_deg", 0.925, "gain_fn", nex_gain, ...
    "backlobe_dBi", 5.5, "scan_width_deg", 360, "modes", table(0.5, 1, 'VariableNames', ["tilt_deg", "prob"]));
asr = struct("name", "ASR-11", "g_max_dBi", 34, "bw_deg", 1.45, "gain_fn", asr_gain, ...
    "backlobe_dBi", 3.5, "scan_width_deg", 360);
v = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
nexv = nex;
nexv.modes = table(v.cuts.el_deg, v.cuts.dur_s, 'VariableNames', ["tilt_deg", "prob"]);

pcts = [50 90 99 99.9 99.99 99.999]';
runs = {"one cut", nex; "VCP 212", nexv};
R = {};
res = struct();
for r = 1:2
    a1 = runs{r, 2};
    res(r).c = macd_rev1(a1, asr, pcts=pcts);                                   % C# settings
    res(r).f = macd_rev1(a1, asr, pcts=pcts, step_div=ov.fine_step_div, bin_dB=ov.fine_bin_dB);
    for s = ["C# (step BW/4, 1 dB bins)", "fine"]
        o = res(r).c;
        st = ov.fine_step_div; bw = ov.fine_bin_dB;
        if s ~= "fine"
            st = 4; bw = 1;
        else
            o = res(r).f;
        end
        for i = 1:numel(pcts)
            R(end + 1, :) = {runs{r, 1}, s, st, bw, pcts(i), o.pct.upper_edge_dBi(i), o.pct.interp_dBi(i)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(R, 'VariableNames', ["case", "settings", "step_div", "bin_dB", "pct", "upper_edge_dBi", "interp_dBi"]), ...
    fullfile(ov.out_dir, "macd_pct_rev1.csv"));

% ---- Resolution: one cut, step and bin width
S = {};
for sd = [4 8 40]
    for bw = [1 0.1]
        o = macd_rev1(nex, asr, pcts=pcts, step_div=sd, bin_dB=bw);
        for i = 3:numel(pcts)
            S(end + 1, :) = {sd, bw, pcts(i), o.pct.upper_edge_dBi(i), o.pct.interp_dBi(i), ...
                res(1).f.pct.interp_dBi(i)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(S, 'VariableNames', ["step_div", "bin_dB", "pct", "upper_edge_dBi", "interp_dBi", "fine_dBi"]), ...
    fullfile(ov.out_dir, "macd_resolution_rev1.csv"));

% ---- Figures
C = struct("dark", [22 51 83]/255, "blue", [86 116 196]/255, "gray", [138 141 145]/255, ...
    "purple", [112 48 160]/255, "grid", [0.9 0.9 0.9]);
fs = {"FontName", "Times New Roman", "FontSize", 8};

% Step 1: sampling
fh = new_fig(6.5, 3.6);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for a = 1:2
    ax = nexttile(tl);
    if a == 1, A = nex; tilt = 0.5; smp = res(1).c.samples{1}; else, A = asr; tilt = 0; smp = res(1).c.samples{2}; end
    az = (-3:0.001:3)';
    plot(ax, az, A.gain_fn(abs(az), tilt), "Color", C.gray, "LineWidth", 1.2); hold(ax, "on");
    k = abs(smp.rel_az_deg) <= 3;
    plot(ax, smp.rel_az_deg(k), smp.gain_dBi(k), "o", "Color", C.blue, "MarkerFaceColor", C.blue, "MarkerSize", 3);
    % The scan's first and last pointings are less than one step apart
    % (numBearings = (int)(360/step) + 1), so the scan start is sampled twice.
    if abs(smp.rel_az_deg(1)) <= 3
        gap = abs(mod(smp.rel_az_deg(end) - smp.rel_az_deg(1) + 180, 360) - 180);
        plot(ax, smp.rel_az_deg([1 end]), smp.gain_dBi([1 end]), "o", "Color", C.purple, "MarkerFaceColor", C.purple, "MarkerSize", 3);
        text(ax, 0.25, smp.gain_dBi(1) + 3, sprintf("Scan start and end, %.2f%s apart", gap, char(176)), fs{:}, ...
            "Color", C.purple, "FontWeight", "bold");
    end
    style_ax(ax, C);
    xlim(ax, [-3 3]); ylim(ax, [0 50]);
    xlabel(ax, "Pointing from the path (degrees)", fs{:});
    if a == 1, ylabel(ax, "Gain toward the other radar (dBi)", fs{:}); end
    title(ax, sprintf("%s: a sample every %.3g%s", A.name, A.bw_deg/4, char(176)), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_sampling_rev1.png"), "Resolution", 300); close(fh);

% Step 2: histograms
fh = new_fig(6.5, 3.4);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for a = 1:2
    ax = nexttile(tl);
    h = res(1).c.hist{a};
    bar(ax, h.bin_lo_dBi + 0.5, h.prob, 1, "FaceColor", C.blue, "EdgeColor", "w");
    ax.YScale = "log";
    style_ax(ax, C);
    ylim(ax, [1e-4 1]); xlim(ax, [0 48]);
    xlabel(ax, "Gain toward the other radar (dBi), 1 dB bins", fs{:});
    if a == 1, ylabel(ax, "Share of pointings", fs{:}); end
    nm = [string(nex.name) string(asr.name)];
    title(ax, sprintf("%s: %d pointings", nm(a), res(1).c.n_samples(a)), ...
        fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_hist_rev1.png"), "Resolution", 300); close(fh);

% Step 3: mode averaging
fh = new_fig(6.5, 3.2);
ax = axes(fh);
H = res(2).c.mode_hist{1};
h = res(2).c.hist{1};
els = nexv.modes.tilt_deg;
[ue, ia] = unique(els, "stable");
show = ue(ismember(ue, [0.5 0.9 1.3]));              % only these cuts reach above the 16.5 dBi sidelobe
for e = show'
    j = ia(ue == e);
    stairs(ax, h.bin_lo_dBi, H(:, j), "Color", C.gray, "LineWidth", 0.8); hold(ax, "on");
    kmax = find(H(:, j) > 0, 1, "last");
    text(ax, h.bin_lo_dBi(kmax) + 1.2, H(kmax, j)*1.6, sprintf("%.1f%s cut", e, char(176)), fs{:}, "Color", C.gray, ...
        "FontWeight", "bold");
end
stairs(ax, h.bin_lo_dBi, h.prob, "Color", C.blue, "LineWidth", 1.6);
text(ax, 7, 2e-1, "Time-weighted average of all cuts", fs{:}, "Color", C.blue, "FontWeight", "bold");
ax.YScale = "log";
style_ax(ax, C);
ylim(ax, [1e-4 1]); xlim(ax, [0 48]);
xlabel(ax, "WSR-88D gain toward the ASR-11 (dBi), 1 dB bins", fs{:});
ylabel(ax, "Share of pointings", fs{:});
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_modes_rev1.png"), "Resolution", 300); close(fh);

% Steps 4 and 5: convolution and percentiles
fh = new_fig(6.5, 3.6);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
for r = 1:2
    ax = nexttile(tl);
    mf = res(r).f.mag; mc = res(r).c.mag;
    k = mf.exceed > 1e-9;
    semilogy(ax, mf.bin_lo_dBi(k), mf.exceed(k), "Color", C.gray, "LineWidth", 1.2); hold(ax, "on");
    k = mc.exceed > 1e-9;
    stairs(ax, mc.bin_lo_dBi(k), mc.exceed(k), "Color", C.blue, "LineWidth", 1.4);
    q = res(r).c.pct(res(r).c.pct.pct >= 99, :);
    plot(ax, q.upper_edge_dBi, 1 - q.pct/100, "o", "Color", C.dark, "MarkerFaceColor", "w", "MarkerSize", 4);
    for i = 1:height(q)
        text(ax, q.upper_edge_dBi(i) - 1, (1 - q.pct(i)/100)*0.4, sprintf("%g%%: %.0f", q.pct(i), q.upper_edge_dBi(i)), ...
            fs{:}, "Color", C.dark, "HorizontalAlignment", "right", "FontSize", 7);
    end
    style_ax(ax, C);
    ax.YScale = "log";
    xlim(ax, [0 85]); ylim(ax, [1e-7 1]);
    xlabel(ax, "MAG (dBi)", fs{:});
    if r == 1
        ylabel(ax, "Share of time at or above", fs{:});
        text(ax, 52, 3e-2, "MACD settings", fs{:}, "Color", C.blue, "FontWeight", "bold");
        text(ax, 52, 1e-2, "Fine reference", fs{:}, "Color", C.gray, "FontWeight", "bold");
    end
    title(ax, ifelse_str(r == 1, "WSR-88D at 0.5° only", "VCP 212, cuts weighted by time"), fs{:}, ...
        "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_conv_rev1.png"), "Resolution", 300); close(fh);

% Setup: top view of the two radars and their beams
fh = new_fig(6.5, 2.6);
ax = axes(fh); hold(ax, "on"); axis(ax, "equal"); axis(ax, "off");
P1 = [0 0]; P2 = [10 0]; L = 3.2; ph = [28 140];      % pointing angles from the path, degrees
plot(ax, [P1(1) P2(1)], [0 0], "--", "Color", C.gray, "LineWidth", 1);
text(ax, 5, -0.45, "Path between the radars", fs{:}, "Color", C.gray, "HorizontalAlignment", "center");
cols = {C.blue, C.dark}; nm = ["WSR-88D", "ASR-11"]; base = [0 180]; P = {P1, P2};
for r = 1:2
    a0 = base(r) + ph(r);
    w = 6;                                             % drawn beam width, for visibility
    t = linspace(a0 - w/2, a0 + w/2, 20);
    patch(ax, [P{r}(1), P{r}(1) + L*cosd(t)], [P{r}(2), P{r}(2) + L*sind(t)], cols{r}, "FaceAlpha", 0.35, "EdgeColor", cols{r});
    plot(ax, P{r}(1), P{r}(2), "o", "MarkerSize", 7, "MarkerFaceColor", cols{r}, "MarkerEdgeColor", cols{r});
    text(ax, P{r}(1) + (r == 2)*0.2, P{r}(2) - 0.9 - (r == 2)*0.6, nm(r), fs{:}, "Color", cols{r}, "FontWeight", "bold", "HorizontalAlignment", "center");
    ta = linspace(base(r), a0, 30);
    plot(ax, P{r}(1) + 1.3*cosd(ta), P{r}(2) + 1.3*sind(ta), "Color", C.purple, "LineWidth", 1);
    am = (base(r) + a0)/2;
    text(ax, P{r}(1) + 1.75*cosd(am), P{r}(2) + 1.75*sind(am), sprintf("\\phi_%d", r), fs{:}, "Color", C.purple, ...
        "FontWeight", "bold", "HorizontalAlignment", "center");
end
text(ax, 5, 2.6, "MAG = G_{WSR-88D}(\phi_1) + G_{ASR-11}(\phi_2)", fs{:}, "Color", C.dark, "HorizontalAlignment", "center", "FontSize", 10);
xlim(ax, [-1.5 11.5]); ylim(ax, [-1.5 3.2]);
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_setup_rev1.png"), "Resolution", 300); close(fh);

% Step 4: every pair of bins, and the MAG histogram
o1 = res(1).c;
h1 = o1.hist{1}; h2 = o1.hist{2};
J = h2.prob*h1.prob.';                                 % rows: ASR-11 bins, columns: WSR-88D bins
fh = new_fig(6.5, 3.4);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(tl);
Jl = log10(J); Jl(J == 0) = NaN;
imagesc(ax, h1.bin_lo_dBi + 0.5, h2.bin_lo_dBi + 0.5, Jl, "AlphaData", ~isnan(Jl));
ax.YDir = "normal"; hold(ax, "on");
colormap(ax, flipud(gray(256))*0.85 + 0.15*[86 116 196]/255);
cb = colorbar(ax); cb.Label.String = "log_{10} share of time"; cb.FontName = "Times New Roman"; cb.FontSize = 7;
for g = [50 60 70]
    x = 5:46; plot(ax, x, g - x, "Color", C.purple, "LineWidth", 0.8);
    text(ax, g - 35.5, 35.5, sprintf("%d dBi", g), fs{:}, "Color", C.purple, "FontWeight", "bold", "HorizontalAlignment", "center");
end
text(ax, 6.5, 5.2, sprintf("Both far sidelobes: %.1f%%", 100*J(1, 1)), fs{:}, "Color", C.dark, "FontWeight", "bold");
style_ax(ax, C); ax.YGrid = "off";
xlim(ax, [4 46]); ylim(ax, [2 37]);
xlabel(ax, "WSR-88D gain bin (dBi)", fs{:}); ylabel(ax, "ASR-11 gain bin (dBi)", fs{:});
title(ax, "(a) Each square: one bin of each antenna", fs{:}, "FontWeight", "normal", "Color", C.dark);
ax = nexttile(tl);
m = o1.mag;
bar(ax, m.bin_lo_dBi + 0.5, m.prob, 1, "FaceColor", C.blue, "EdgeColor", "none");
ax.YScale = "log"; style_ax(ax, C);
ylim(ax, [1e-7 1]); xlim(ax, [0 80]);
xlabel(ax, "MAG bin (dBi): sum of the two bins", fs{:}); ylabel(ax, "Share of time", fs{:});
title(ax, "(b) Add up the squares on each diagonal", fs{:}, "FontWeight", "normal", "Color", C.dark);
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_convgrid_rev1.png"), "Resolution", 300); close(fh);

% Step 5: reading the 99% and 99.999% values
fh = new_fig(6.5, 3.2);
tl = tiledlayout(fh, 1, 2, "TileSpacing", "compact", "Padding", "compact");
zoom = {[30 40], [1e-3 1e-1], 99; [68 77], [1e-7 1e-4], 99.999};
for z = 1:2
    ax = nexttile(tl);
    q = zoom{z, 3};
    x = m.bin_lo_dBi; ex = m.exceed;
    k = ex > 0;
    stairs(ax, x(k), ex(k), "Color", C.blue, "LineWidth", 1.4); hold(ax, "on");
    xe = [x(k); x(find(k, 1, "last")) + 1]; ee = [ex(k); 0];
    plot(ax, xe, max(ee, 1e-12), "--", "Color", C.gray, "LineWidth", 0.9);
    yline(ax, 1 - q/100, "Color", C.purple, "LineWidth", 0.9);
    r = o1.pct(o1.pct.pct == q, :);
    plot(ax, r.upper_edge_dBi, 1 - q/100, "s", "Color", C.dark, "MarkerFaceColor", "w", "MarkerSize", 5);
    plot(ax, r.interp_dBi, 1 - q/100, "o", "Color", C.dark, "MarkerFaceColor", C.dark, "MarkerSize", 4);
    text(ax, r.upper_edge_dBi + 0.2, (1 - q/100)*1.5, sprintf("Bin upper edge: %.0f dBi", r.upper_edge_dBi), fs{:}, "Color", C.dark);
    text(ax, r.interp_dBi - 0.2, (1 - q/100)*0.6, sprintf("Interpolated: %.1f dBi", r.interp_dBi), fs{:}, "Color", C.dark, ...
        "HorizontalAlignment", "right");
    style_ax(ax, C); ax.YScale = "log";
    xlim(ax, zoom{z, 1}); ylim(ax, zoom{z, 2});
    xlabel(ax, "MAG (dBi)", fs{:});
    if z == 1, ylabel(ax, "Share of time at or above", fs{:}); end
    title(ax, sprintf("%g%%: the line at %.3g crosses one bin", q, 1 - q/100), fs{:}, "FontWeight", "normal", "Color", C.dark);
end
exportgraphics(fh, fullfile(ov.out_dir, "fig_macd_lookup_rev1.png"), "Resolution", 300); close(fh);

% -------------------------------------------------------------------------
function fh = new_fig(w, h)
fh = figure("Color", "w", "Units", "inches", "Visible", "off");
if ~isMATLABReleaseOlderThan("R2025a")
    theme(fh, "light");
end
fh.Position(3:4) = [w h];
end

function style_ax(ax, C)
box(ax, "off");
ax.FontName = "Times New Roman";
ax.FontSize = 8;
ax.TickDir = "out";
ax.XColor = C.dark;
ax.YColor = C.dark;
ax.YGrid = "on";
ax.GridColor = C.grid;
ax.GridAlpha = 1;
end

function s = ifelse_str(c, a, b)
if c, s = a; else, s = b; end
end
