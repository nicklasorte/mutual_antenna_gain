%REPORT_NUMBERS_REV14 Tables and figure for methodology report rev14.
%   Answers the independent recheck of revision 13
%   (docs/rev13_independent_review_rev1.md, findings Y1 to Y12). Reads the
%   per-trial results of revision 13 (no rerun of the 20-seed study) and adds:
%     - the share of high volumes between the repeating volumes of each long
%       periodic run, against the trial's own share, and the span rate that
%       independent volumes would give (Y1);
%     - the chance comparison with the sampling of the 200 trial shares
%       included: the chance of a long run is found on a grid of shares,
%       then the 200 shares are resampled; also the chance without the
%       trials that hold long runs (Y2);
%     - the share of trials that reach a long run as a function of the
%       lock threshold, from the slip alone (Y4);
%     - hits and recurrence intervals of the long-period runs (Y7);
%     - the random-walk chance table with the per-trial chance in percent
%       (Y10);
%     - Figure 15 with the x axis to 25 degrees (Y5).
%   Reads lock_trials_rev13.csv, seeds_runs_rev13.csv, lock_fit_rev13.csv,
%   lock_long_k_rev13.csv, lock_walk_null_rev13.csv, and ts_all_modes_rev8.csv
%   from in_dir.
%
%     lock_offlattice_rev14.csv   long periodic runs: high volumes between
%                                 the repeating ones, rates
%     lock_share_drift_rev14.csv  trials within each lock threshold and how
%                                 many of them reach a long periodic run
%     lock_null_grid_rev14.csv    chance of a long run (periods up to 4) on
%                                 a grid of shares of high volumes
%     lock_null_boot_rev14.csv    chance of a long run with the trial shares
%                                 resampled, without the long-run trials,
%                                 and the ratio ranges
%     lock_long_k_rev14.csv       the long-period trials with hits and
%                                 recurrence intervals
%     lock_walk_null_rev14.csv    random-walk maxima against chance
%     fig_lock_any16_rev14.png    Figure 15
%
%   Set rev14_overrides (struct: out_dir, in_dir, n_grid_sims, p_step,
%   n_boot, parallel_flag, cells) for a small run. cells is an N-by-2 matrix
%   of [el_deg threshold_dBi] rows (default all four).
%
%   DRAFT//PRE-DECISIONAL
%
%   See also report_numbers_rev13, null_periodic_trial_rev1, periodic_runs_rev1.

if ~exist("rev14_overrides", "var"), rev14_overrides = struct(); end
here = fileparts(which("report_numbers_rev14"));
ov = struct("out_dir", fullfile(here, "docs", "figures"), "in_dir", fullfile(here, "docs", "figures"), ...
    "n_grid_sims", 20000, "p_step", 0.01, "n_boot", 10000, "parallel_flag", 1, ...
    "cells", [0 60; 0 70; 0.5 60; 0.5 70]);
f = fieldnames(rev14_overrides);
for i = 1:numel(f)
    ov.(f{i}) = rev14_overrides.(f{i});
end
out_dir = ov.out_dir;
if ~isfolder(out_dir), mkdir(out_dir); end
if ov.parallel_flag == 1 && isempty(gcp("nocreate"))
    parpool("Processes");
end
use_par = ov.parallel_flag == 1 && ~isempty(gcp("nocreate"));

T_lock = readtable(fullfile(ov.in_dir, "lock_trials_rev13.csv"), "TextType", "string");
T_fit = readtable(fullfile(ov.in_dir, "lock_fit_rev13.csv"), "TextType", "string");
T_all = readtable(fullfile(ov.in_dir, "ts_all_modes_rev8.csv"), "TextType", "string");
modes = unique(T_lock.mode, "stable");
t_vol = @(m) T_all.t_volume_s(find(T_all.mode == m, 1));

% ---- Y1: high volumes between the repeating volumes of each long run.
% The run spans per_vols volumes, per_hits of them on the period; the
% other per_vols - per_hits are between. Independent volumes with the
% trial's own share p would make p of those between high.
T_lock.between_vols = T_lock.per_vols - T_lock.per_hits;
T_lock.between_high = T_lock.per_span_high - T_lock.per_hits;
T_lock.between_share = T_lock.between_high./T_lock.between_vols;       % NaN when the run has no volumes between
T_lock.trial_share = T_lock.trial_vol_pct/100;
O = {};
for m = modes'
    tv = t_vol(m);
    for a = [0 0.5]
        for th = [60 70]
            L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
            span_h = L.per_span_high./max(L.per_vols, 1)*3600/tv;
            indep_h = (L.per_hits + L.trial_share.*L.between_vols)./max(L.per_vols, 1)*3600/tv;
            avg_h = mean(L.trial_share)*3600/tv;
            for H = [4 8]
                lg = L.per_vols*tv/3600 >= H;
                if ~any(lg)
                    continue
                end
                O(end + 1, :) = {m, a, th, H, nnz(lg), median(L.between_share(lg), "omitnan"), ...
                    median(L.trial_share(lg)), median(L.trial_share), ...
                    nnz(lg & L.between_share < L.trial_share), avg_h, median(span_h(lg)), ...
                    median(indep_h(lg)), nnz(lg & span_h > avg_h), nnz(lg & span_h < indep_h)}; %#ok<SAGROW>
            end
        end
    end
end
if isempty(O), O = cell(0, 14); end
writetable(cell2table(O, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_long", ...
    "median_between_share_long", "median_trial_share_long", "median_trial_share_all", ...
    "n_long_between_below_own_share", "mode_avg_high_per_h", "median_span_per_h_long", ...
    "median_indep_span_per_h_long", "n_long_span_above_avg", "n_long_span_below_indep"]), ...
    fullfile(out_dir, "lock_offlattice_rev14.csv"));

% ---- Y4: share of trials within each lock threshold (slip alone, any 1 to
% 16 volumes) that reach a long periodic run (periods up to 4). Every slip
% has a lock within 360/17 = 21.2 deg over 1 to 16 volumes.
S = {};
thr_deg = [0.5 1 2 3 5 8 10 15 360/17];
for m = modes'
    tv = t_vol(m);
    L = T_lock(T_lock.mode == m & T_lock.el_deg == 0.5 & T_lock.threshold_dBi == 60, :);
    for H = [4 8]
        lg = L.per_vols*tv/3600 >= H;
        for d = thr_deg
            tight = L.lock_any16_drift_deg <= d;
            S(end + 1, :) = {m, 0.5, 60, H, d, nnz(tight), nnz(tight & lg), nnz(tight & lg)/max(nnz(tight), 1), ...
                100*nnz(tight)/height(L), nnz(lg), height(L)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(S, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "drift_max_deg", ...
    "n_tight", "n_tight_long", "share_tight_long", "pct_trials_tight", "n_long", "n_trials"]), ...
    fullfile(out_dir, "lock_share_drift_rev14.csv"));

% ---- Y2: chance of a long run (periods up to 4) on a grid of shares, at
% each cell's median trial length; then the 200 trial shares are resampled
% (bootstrap), with the grid's own Monte Carlo error drawn too.
G = {};
Bt = {};
M = ov.n_grid_sims;
cell_id = 0;
for m = modes'
    tv = t_vol(m);
    for c = 1:size(ov.cells, 1)
        a = ov.cells(c, 1);
        th = ov.cells(c, 2);
        cell_id = cell_id + 1;
        L = T_lock(T_lock.mode == m & T_lock.el_deg == a & T_lock.threshold_dBi == th, :);
        pv = L.trial_share;
        n_v = round(median(L.trial_n_vols));
        pg = (floor(min(pv)/ov.p_step):ceil(max(pv)/ov.p_step))'*ov.p_step;
        pg = min(max(pg, 0), 1);
        Ns = ceil([4 8]*3600/tv);
        q = zeros(numel(pg), numel(Ns));
        for g = 1:numel(pg)
            r4 = zeros(M, 1);
            base = 1e9*cell_id + 1e6*g;
            pgg = pg(g);
            if use_par
                parfor i = 1:M
                    r4(i) = null_periodic_trial_rev1(pgg, n_v, i, base, 20261114);
                end
            else
                for i = 1:M
                    r4(i) = null_periodic_trial_rev1(pgg, n_v, i, base, 20261114);
                end
            end
            q(g, :) = mean(r4 >= Ns, 1);
            G(end + 1, :) = {m, a, th, n_v, pg(g), M, 100*q(g, 1), 100*q(g, 2)}; %#ok<SAGROW>
        end
        % Bootstrap: resample the trials, redraw the grid chances.
        prev = rng(20261115 + cell_id, "twister");
        nb = ov.n_boot;
        nullb = zeros(nb, numel(Ns));
        nt = numel(pv);
        for b = 1:nb
            idx = randi(nt, nt, 1);
            qb = binornd(M, q)/M;
            for h = 1:numel(Ns)
                nullb(b, h) = mean(interp1(pg, qb(:, h), pv(idx), "linear"));
            end
        end
        rng(prev);
        for h = 1:numel(Ns)
            H = [4 8];
            H = H(h);
            qj = interp1(pg, q(:, h), pv, "linear");
            lg = L.per_vols >= Ns(h);
            fr = T_fit(T_fit.mode == m & T_fit.el_deg == a & T_fit.threshold_dBi == th & T_fit.hours == H & ...
                T_fit.max_period == 4, :);
            lo = prctile(nullb(:, h), 2.5);
            hi = prctile(nullb(:, h), 97.5);
            [~, iq] = sort(qj, "descend");
            top10 = sum(qj(iq(1:min(10, nt))))/max(sum(qj), eps);
            Bt(end + 1, :) = {m, a, th, H, Ns(h), fr.trials_periodic_ge_pct, fr.periodic_lo_pct, fr.periodic_hi_pct, ...
                fr.null_periodic_ge_pct, 100*mean(qj), 100*lo, 100*hi, 100*mean(qj(~lg)), nnz(lg), 100*top10, ...
                fr.periodic_lo_pct/(100*hi), fr.periodic_hi_pct/(100*lo)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(G, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "n_vols", "share", "n_sims", ...
    "null_ge4h_pct", "null_ge8h_pct"]), fullfile(out_dir, "lock_null_grid_rev14.csv"));
writetable(cell2table(Bt, 'VariableNames', ["mode", "el_deg", "threshold_dBi", "hours", "n_vols", ...
    "model_pct", "model_lo_pct", "model_hi_pct", "null_mc_rev13_pct", "null_grid_pct", "null_boot_lo_pct", ...
    "null_boot_hi_pct", "null_without_long_pct", "n_long", "top10_share_of_null_pct", "ratio_lo", "ratio_hi"]), ...
    fullfile(out_dir, "lock_null_boot_rev14.csv"));

% ---- Y7: hits and recurrence intervals of the long-period runs.
K = readtable(fullfile(ov.in_dir, "lock_long_k_rev13.csv"), "TextType", "string");
tvk = arrayfun(t_vol, K.mode);
K.interval_hours = K.per16_k.*tvk/3600;
K.hits_per_h = K.per16_hits./K.per16_hours;
writetable(K, fullfile(out_dir, "lock_long_k_rev14.csv"));

% ---- Y10: random-walk table with the per-trial chance in percent.
W = readtable(fullfile(ov.in_dir, "lock_walk_null_rev13.csv"), "TextType", "string");
W.null_trial_ge_pct = 100*W.null_trial_ge_pct;              % revision 13 wrote a fraction under this name
W.chance_any_of_trials_pct = 100*W.chance_any_of_trials;
W.chance_any_of_trials = [];
writetable(W, fullfile(out_dir, "lock_walk_null_rev14.csv"));

report_figures_rev14(out_dir, ov.in_dir);
