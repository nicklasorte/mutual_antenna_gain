%EXAMPLE_NEXRAD_ASR11_REV1 WSR-88D vs ASR-11 mutual antenna gain, all VCPs.
%   Runs nexrad_asr11_mag_rev1 for the 2.7 GHz working paper baseline
%   (F.1245, NEXRAD at 3 rpm on the horizon) and for every operational VCP
%   with its scan options, then ranks the VCP modes by worst-case MAG.
%   All antenna and scan assumptions, with URL citations, are in
%   nexrad_asr11_mag_rev1 (assumptions, citations) and nexrad_vcp_rev1.
%
%   Worst-case metrics (path ASR-11 -> NEXRAD, g21; g12 is equal):
%     p99, p99_9, p99_99  time-weighted MAG percentiles (random mode)
%     rev_pct_ge_X        % of ASR-11 revolutions whose maximum MAG >= X dBi
%     vol_max             largest per-volume maximum in the time series
%     t_low_pct           % of NEXRAD volume time at the 0.5 deg cuts

clearvars -except example_overrides; close all;   % example_overrides: small runs (tests)

base = struct("seed", 20261105, "n_samples", 5e6, "duration_s", 3600, "n_trials", 4, "nexrad_pattern", "nssl_measured", "cdf_method", "monte_carlo");   % Rev 5 used the measured pattern   % M.2069-0: 5e6 samples
if exist("example_overrides", "var")
    base = merge_overrides_rev1(base, example_overrides);
end
thr = [40 50 60];                                  % per-revolution thresholds (dBi)

% ---- Case list: the paper baseline, then every VCP mode (74), with the
%      option combinations the RDA allows (vcp_mode_cases_rev1)
c = vcp_mode_cases_rev1("lowest");      % the 74 modes of the Rev 5 study
fprintf("%d cases (%d VCP modes and the paper baseline)\n", numel(c), numel(c) - 1);

% ---- Run
env_pct = [99 99.9 99.99 99.999];
summ = cell(1, numel(c));
n = numel(c);
rows = cell(n, 1);
res = cell(n, 1);
for i = 1:n
    o = base;
    f = setdiff(fieldnames(c{i}), "label");
    for k = 1:numel(f)
        o.(f{k}) = c{i}.(f{k});
    end
    r = nexrad_asr11_mag_rev1(o);
    % Keep memory bounded over many modes: keep a compact summary for the
    % envelope (mag_mode_summary_rev1) and drop the samples, except for
    % the paper baseline and VCP 215, which the figures plot. The worst
    % mode is rerun after the sweep with the same seed.
    summ{i} = mag_mode_summary_rev1(r.cdf, env_pct);
    r.cdf.g12_dBi = [];
    if i == 1 || c{i}.label == "VCP 215"
        r.cdf.g21_dBi = single(r.cdf.g21_dBi);
    else
        r.cdf.g21_dBi = [];
    end
    res{i} = r;
    q = r.cdf.pct21;
    R = r.scan_max.asr_rev(r.scan_max.asr_rev.complete, :);
    s = r.sys_nexrad.scan;
    if isfield(s, "cuts")
        t_low = 100*sum(s.cuts.dur_s(s.cuts.el_deg == 0.5))/s.volume_s;
        t_vol = s.volume_s;
    else
        t_low = NaN;
        t_vol = 20;
    end
    rows{i} = {string(c{i}.label), t_vol, t_low, q.gain_dBi(q.pct == 99), ...
        q.gain_dBi(q.pct == 99.9), q.gain_dBi(q.pct == 99.99), ...
        100*mean(R.max_g21_dBi >= thr(1)), 100*mean(R.max_g21_dBi >= thr(2)), ...
        100*mean(R.max_g21_dBi >= thr(3)), max(r.scan_max.nexrad_volume.max_g21_dBi)};
    fprintf("%-34s done\n", c{i}.label);
end
T = cell2table(vertcat(rows{:}), "VariableNames", ["case", "volume_s", "t_low_pct", ...
    "p99", "p99_9", "p99_99", "rev_pct_ge_40", "rev_pct_ge_50", "rev_pct_ge_60", "vol_max"]);

% ---- Rank the NOAA scan modes (paper baseline excluded) by worst case
V = T(2:end, :);
fprintf("\nVCP modes ranked by 99.9%% MAG (worst first):\n");
disp(sortrows(V, ["p99_9", "p99"], "descend"))
fprintf("\nPaper baseline:\n");
disp(T(1, :))
[~, i99] = max(V.p99);
[~, i999] = max(V.p99_9);
[~, irev] = max(V.rev_pct_ge_50);
fprintf("Worst case: p99 %s; p99.9 %s; revolutions >= 50 dBi %s\n", ...
    V.case(i99), V.case(i999), V.case(irev));

% ---- Envelope over all VCP modes: bounds any mix of modes a site uses
%      (mag_mode_envelope_rev1), so the unknown mode statistics drop out.
env = mag_mode_envelope_rev1(summ(2:end), V.case', env_pct);

% ---- Rerun the worst mode (same seed: identical samples) to keep its
%      samples for the figures.
o = base;
f = setdiff(fieldnames(c{1 + i999}), "label");
for k = 1:numel(f)
    o.(f{k}) = c{1 + i999}.(f{k});
end
rw = nexrad_asr11_mag_rev1(o);
rw.cdf.g12_dBi = [];
rw.cdf.g21_dBi = single(rw.cdf.g21_dBi);
res{1 + i999} = rw;
fprintf("\nMode envelope (bound for any mix of the %d VCP modes):\n", height(V));
disp(env.pct(:, ["pct", "envelope_dBi", "worst_mode"]))

% ---- Detail for the worst mode: assumptions, citations, per-scan maxima
w = res{1 + i999};
fprintf("\nSensitivity of the MAG choice, %s:\n", V.case(i999));
disp(w.sensitivity)
disp(w.assumptions)
disp(w.citations(:, ["key", "url"]))
disp(w.scan_max.by_elevation)
disp(w.scan_max.asr_rev_pct)

% ---- Plots
figure("Name", "MAG exceedance by mode");
hold on; grid on;
pick = unique([1, 1 + i999, 1 + i99, find(T.case == "VCP 215"), find(T.case == "VCP 212")]);
for i = pick
    semilogy(res{i}.cdf.cdf21.gain_dBi, res{i}.cdf.cdf21.exceed_pct, "LineWidth", 1.5);
end
set(gca, "YScale", "log");
ylim([1e-4 100]);
% Main-beam / sidelobe reference levels of the worst mode (both
% combinations), from nexrad_asr11_mag_rev1 out.ref_levels
R = w.ref_levels;
sty = ["-", "--", ":", "--", ":"];
for i = 1:height(R)
    xline(R.mag_dBi(i), sty(i), sprintf("%s (%.1f dBi)", R.combination(i), R.mag_dBi(i)), ...
        "LabelOrientation", "aligned", "LabelVerticalAlignment", "bottom", ...
        "FontSize", 8, "HandleVisibility", "off");
end
xlabel("Mutual antenna gain (dBi)");
ylabel("Time exceeded (%)");
legend(T.case(pick), "Location", "southwest");
title("WSR-88D and ASR-11 mutual antenna gain");
disp(R)

figure("Name", "Max MAG per ASR-11 revolution");
R = w.scan_max.asr_rev(w.scan_max.asr_rev.trial == 1, :);
plot(R.t_start_s, R.max_g21_dBi, ".");
grid on;
xlabel("Time (s)");
ylabel("Maximum MAG per revolution (dBi)");
title(sprintf("%s: maximum per 4.8 s ASR-11 scan", V.case(i999)));
