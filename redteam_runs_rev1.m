function rows = redteam_runs_rev1(ids, opts)
%REDTEAM_RUNS_REV1 Red-team attacks on both MAG models (plan step 8).
%   rows = redteam_runs_rev1(ids) runs the attacks of the plan's section 7.2
%   (docs/mag_optimization_redteam_plan_rev1.md) named in ids and returns
%   one register row per attack. Each attack names a claim, a test and a
%   pass criterion fixed before the run; it writes its evidence to
%   <out_dir>/<id>_rev1.csv and its register row to <out_dir>/<id>_row.csv.
%   ids = "register" assembles docs/redteam_register_rev1.csv from the row
%   files plus the attacks closed in earlier steps (D1, D3, D4, E1, E2, F2),
%   the blocked ones (F3) and the scope limits (G).
%
%   ids   string array: "A1" "A2" "A3" "A4" "B1" "B2" "B4" "B5" "C1" "C2"
%         "C3" "C4" "C5" "D2" "F1" "NOAA", "all" (every attack in that
%         order), or "register". B2 also answers B3.
%   opts  struct, all fields optional:
%     .repo_dir  repository folder (default: this file's folder)
%     .out_dir   evidence folder (default <repo_dir>/docs/figures/redteam)
%     .quick     true: tiny runs that check the code paths only (default
%                false; the tests use it). Quick rows are not evidence.
%     .inputs    "rev5" (default: NSSL measured pattern and the 74 modes,
%                so the register reproduces) | "rev6" (ROC specification
%                pattern and all 216 modes: the revision 6 defaults, red-team
%                finding P7). With "rev6" the default out_dir is
%                <repo_dir>/docs/figures/redteam_rev6/attacks and "register"
%                is not available.
%
%   Run (long; about an hour in total on a laptop, so in batch):
%     j = batch(parcluster("Processes"), @redteam_runs_rev1, 1, ...
%         {["A3" "A4" "B4" "C2" "C3" "C4" "C5"]}, "AutoAddClientPath", false);
%     redteam_runs_rev1("register")
%
%   Sources: public and cleared material only. Attacks that need restricted
%   or unavailable data name it and stop (C3, C5, F3).
%
%   See also redteam_fault_seed_rev1, mag_pair_cdf_rev1, mag_exact_cdf_rev1,
%   nexrad_asr11_mag_rev1.

arguments
    ids (1, :) string = "all"
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
if ~isfolder(opts.out_dir)
    mkdir(opts.out_dir);
end
all_ids = ["A1" "A2" "A3" "A4" "B1" "B2" "B4" "B5" "C1" "C2" "C3" "C4" "C5" "D2" "F1" "NOAA"];
if isequal(ids, "register")
    if opts.inputs == "rev6"
        error("redteam:registerRev5", "The register holds the revision 5 runs; run ""register"" with inputs ""rev5"".");
    end
    rows = make_register(opts);
    return
end
if isequal(ids, "all")
    ids = all_ids;
end
bad = setdiff(ids, all_ids);
if ~isempty(bad)
    error("redteam:unknownAttack", "Unknown attack(s): %s.", strjoin(bad, ", "));
end
ctx = base_context(opts);
rows = cell(numel(ids), 1);
for i = 1:numel(ids)
    t0 = tic;
    r = feval("run_" + lower(ids(i)), ctx, opts);
    r.runtime_s = toc(t0);
    r.quick = opts.quick;
    writetable(struct2table(r, "AsArray", true), fullfile(opts.out_dir, ids(i) + "_row.csv"));
    rows{i} = r;
    fprintf("%s: %s (%.0f s)\n", ids(i), r.result, r.runtime_s);
end
rows = struct2table([rows{:}]', "AsArray", true);
end

% =========================================================================
function ctx = base_context(opts)
% Systems, scan modes and gain tables shared by the attacks.
% The step 8 evidence used the measured pattern (the default then) and the
% 74 modes; opts.inputs = "rev5" keeps them so the register reproduces, and
% "rev6" uses the revision 6 defaults. C1 and the NOAA case set the pattern.
if opts.inputs == "rev6"
    [ctx.pattern, mode_set] = deal("roc_spec", "all");
else
    [ctx.pattern, mode_set] = deal("nssl_measured", "lowest");
end
ob = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
    "cdf_method", "exact", "nexrad_pattern", ctx.pattern));
ctx.sys_n = ob.sys_nexrad;
ctx.sys_a = ob.sys_asr11;
[ctx.scans, ctx.names] = nexrad_mode_scans_rev1(mode_set);
c = vcp_mode_cases_rev1(mode_set);
ctx.cases = c(2:end);
ctx.k_worst = find(ctx.names == "VCP 212 SAILSx3 AVSET 6.4", 1);
ctx.pct = [50 90 99 99.9 99.99 99.999];
if opts.quick
    grid = -1:0.05:1;
else
    grid = -5:0.05:10;
end
ctx.T_n = gain_dist_table_rev1(ctx.sys_n.tx, ctx.scans, struct("name", "WSR-88D", ...
    "mode_names", ctx.names, "el_grid_deg", grid));
ctx.T_a = gain_dist_table_rev1(ctx.sys_a.tx, ctx.sys_a.scan, struct("name", "ASR-11", "el_grid_deg", grid));
end

% -------------------------------------------------------------------------
function o = mode_opts(ctx, k, extra)
% nexrad_asr11_mag_rev1 options for scan mode k plus the fields of extra.
o = rmfield(ctx.cases{k}, "label");
o.nexrad_pattern = ctx.pattern;                  % see base_context
o.cdf_method = "monte_carlo";                    % the step 8 method; attacks set "exact" where they use it
f = fieldnames(extra);
for i = 1:numel(f)
    o.(f{i}) = extra.(f{i});
end
end

% -------------------------------------------------------------------------
function r = row(id, claim, test, criterion, result, key, evidence, severity, next)
r = struct("id", string(id), "claim", string(claim), "test", string(test), ...
    "criterion", string(criterion), "result", string(result), "key_result", string(key), ...
    "evidence", string(evidence), "severity", string(severity), "next_step", string(next), ...
    "runtime_s", NaN, "quick", false);
end

% -------------------------------------------------------------------------
function f = evid(opts, id)
f = fullfile(opts.out_dir, id + "_rev1.csv");
end

function s = rel(opts, f)
s = string(erase(f, opts.repo_dir + filesep));
end

% =========================================================================
% A. Is the convolution valid?
% =========================================================================
function r = run_a1(ctx, opts)
% Exact (independent-pointing) percentiles against time-series percentiles
% of the shared scan clock, every mode. Five independent trials give a
% t-interval for the time-series value; the exact value is conservative
% by at most its bias, so the truth lies in [exact - bias, exact].
pct = [50 90 99 99.9];
if opts.quick
    modes = [ctx.k_worst 1];
    n_tr = 2;
    dur = 120;
else
    modes = 1:numel(ctx.names);
    n_tr = 5;
    dur = 1800;
end
n_m = numel(modes);
n_p = numel(pct);
ex = NaN(n_m, n_p);
ts = NaN(n_m, n_p, n_tr);
bias = NaN(n_m, 1);
for ii = 1:n_m
    for j = 1:n_tr
        o = mode_opts(ctx, modes(ii), struct("seed", 1, "first_trial", j, "n_trials", 1, ...
            "duration_s", dur, "cdf_method", "exact", "n_samples", 1e3, "pct", pct, "ts_method", "fixed"));
        out = nexrad_asr11_mag_rev1(o);
        ts(ii, :, j) = out.ts_pct.gain_dBi';
        if j == 1
            ex(ii, :) = out.cdf.pct21.gain_dBi';
            bias(ii) = out.cdf.bias_dB;
        end
    end
end
m = mean(ts, 3);
s = std(ts, 0, 3);
half = tinv(0.975, n_tr - 1)*s/sqrt(n_tr);
tol = 0.011;                                     % time-series histogram step (0.01 dB) plus round-off
lo = ex - bias - tol;
hi = ex + tol;
pass = (m + half >= lo) & (m - half <= hi);
[mm, pp] = ndgrid(ctx.names(modes)', pct);
T = table(mm(:), pp(:), ex(:), m(:), s(:), half(:), pass(:), m(:) - ex(:), 'VariableNames', ...
    ["mode", "pct", "exact_dBi", "ts_mean_dBi", "ts_sd_dB", "ts_ci95_half_dB", "pass", "ts_minus_exact_dB"]);
writetable(T, evid(opts, "A1"));
n_fail = sum(~pass, "all");
bad_modes = unique(T.mode(~T.pass));
by_p = compose("%g%%: %d", pct', sum(~pass, 1)');
key = sprintf("%d of %d mode-percentile checks outside the 95%% interval (%s). Largest |ts - exact| %.2f dB.", ...
    n_fail, numel(pass), strjoin(by_p, ", "), max(abs(m - ex), [], "all"));
if n_fail == 0
    res = "Pass";
    sev = "None";
    nxt = "None";
elseif n_fail <= ceil(0.05*numel(pass))
    res = "Pass (within the 5% expected by chance)";
    sev = "Low";
    nxt = "Failing modes: " + strjoin(bad_modes, "; ");
else
    res = "Fail";
    sev = "High";
    nxt = "Failing modes: " + strjoin(bad_modes, "; ") + ". Use the time series for these modes.";
end
r = row("A1", "Pointing states are independent (convolution valid)", ...
    sprintf("Exact CDF vs time-percentile CDF of the scan clock, %d modes, %d trials of %g s, 0/0 deg", n_m, n_tr, dur), ...
    "Exact value inside the time-series 95% interval (allowing the exact bias) for every mode", ...
    res, key, rel(opts, evid(opts, "A1")), sev, nxt);
end

% -------------------------------------------------------------------------
function r = run_a2(ctx, opts)
% ATC rotation rate searched across the tolerance band at a fixed rate (no
% jitter): a rate whose revolution time repeats with the WSR-88D cut times
% would make the time series depart from the exact (independent) CDF.
if opts.quick
    rates = [74.5 75 75.5];
    dur = 120;
    n_tr = 1;
    k_top = ctx.k_worst;
else
    rates = 74.5:0.025:75.5;
    dur = 3600;
    n_tr = 2;
    k_top = top_modes(ctx, 3);
end
pct = [99 99.9];
rows_ = {};
for k = k_top
    sn = ctx.sys_n;
    sn.scan = ctx.scans{k};
    ex = mag_exact_cdf_rev1(sn, ctx.sys_a, struct("pct", pct));
    for rr = rates
        sa = ctx.sys_a;
        if isfield(sa.scan, "rate_rpm")
            sa.scan = rmfield(sa.scan, "rate_rpm");
        end
        sa.scan.rate_deg_s = rr;
        sa.scan.rate_jitter_deg_s = 0;
        for j = 1:n_tr
            w = warning("off", "all");
            out = radar_mutual_gain_cdf_rev1(sn, sa, struct("seed", j, "mode", "time_series", ...
                "n_trials", 1, "duration_s", dur, "pct", pct));
            warning(w);
            d = out.pct21.gain_dBi' - ex.pct21.gain_dBi';
            rows_(end + 1, :) = {ctx.names(k), rr, j, ex.pct21.gain_dBi(1), out.pct21.gain_dBi(1), ...
                ex.pct21.gain_dBi(2), out.pct21.gain_dBi(2), d(1), d(2)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows_, 'VariableNames', ["mode", "atc_rate_deg_s", "trial", "exact_p99_dBi", ...
    "ts_p99_dBi", "exact_p99_9_dBi", "ts_p99_9_dBi", "dev_p99_dB", "dev_p99_9_dB"]);
writetable(T, evid(opts, "A2"));
flag = abs(T.dev_p99_9_dB) > 1;
key = sprintf("%d rates x %d modes x %d trials of %g s. Largest |deviation| at 99.9%%: %.2f dB (%s, %.3f deg/s); %d runs over 1 dB.", ...
    numel(rates), numel(k_top), n_tr, dur, max(abs(T.dev_p99_9_dB)), ...
    T.mode(find(abs(T.dev_p99_9_dB) == max(abs(T.dev_p99_9_dB)), 1)), ...
    T.atc_rate_deg_s(find(abs(T.dev_p99_9_dB) == max(abs(T.dev_p99_9_dB)), 1)), sum(flag));
if ~any(flag)
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    res = "Fail";
    sev = "Medium";
    nxt = "Rates over 1 dB: " + strjoin(compose("%.3f", unique(T.atc_rate_deg_s(flag))'), ", ") + ...
        " deg/s. Run the time series for pairs whose ATC rate is known to sit there.";
end
r = row("A2", "The ATC rate tolerance does not matter", ...
    "Fixed ATC rates across 74.5 to 75.5 deg/s (0.025 deg/s step), time series vs exact CDF, worst three modes", ...
    "No rate where clock and exact CDFs differ by more than 1 dB at 99.9%", res, key, rel(opts, evid(opts, "A2")), sev, nxt);
end

% -------------------------------------------------------------------------
function k = top_modes(ctx, n)
% Modes with the largest exact 99.9% MAG at 0/0.
R = mag_pair_cdf_rev1(ctx.T_n, ctx.T_a, 0, 0, struct("pct", 99.9, "per_mode", true));
[~, o] = sort(squeeze(R.mode_pct(1, :, 1)), "descend");
k = o(1:n);
end

% -------------------------------------------------------------------------
function [C, k0, dg] = tilt_mag(ctx, el_n, el_a)
% MAG distribution (path ATC -> WSR-88D) for each WSR-88D tilt alone at
% the table node nearest el_n, el_a: column t is conv(p_tilt, p_atc).
Tn = ctx.T_n;
Ta = ctx.T_a;
[~, in] = min(abs(Tn.el_grid_deg - el_n));
[~, ia] = min(abs(Ta.el_grid_deg - el_a));
Pn = Tn.P_rx;
if isempty(Pn)
    Pn = Tn.P_tx;
end
pn = double(reshape(Pn(:, in, :), size(Pn, 1), []));
pa = double(Ta.P_tx(:, ia, 1));
C = zeros(size(pn, 1) + numel(pa) - 1, size(pn, 2));
for t = 1:size(pn, 2)
    C(:, t) = conv(pn(:, t), pa);
end
k0 = Tn.k0 + Ta.k0;
dg = Tn.bin_dB;
end

function q = pct_of(P, k0, dg, pct)
% Percentiles of the columns of P (pmf by columns): smallest gain with
% exceedance <= 1 - p, as mag_pair_cdf_rev1.
E = cumsum(P, 1, "reverse") - P;
gain = (k0 + (0:size(P, 1) - 1)')*dg;
q = zeros(numel(pct), size(P, 2));
for j = 1:numel(pct)
    [~, b] = max(E <= 1 - pct(j)/100 + 1e-13, [], 1);
    q(j, :) = gain(b);
end
end

function r = run_a3(ctx, opts)
% (1) Random mixtures of modes never exceed the envelope. (2) A mode change
% in the middle of a volume: cuts 1..s of mode X, then a full volume of Y.
[C, k0, dg] = tilt_mag(ctx, 0, 0);
W = ctx.T_n.W;
PM = C*W';
E_env = max(cumsum(PM, 1, "reverse") - PM, [], 2);
rs = RandStream("mrg32k3a", "Seed", 11);
n_mix = 5000;
if opts.quick
    n_mix = 50;
end
g = -log(rand(rs, size(W, 1), n_mix));           % Dirichlet(1) weights over modes
wm = (g./sum(g, 1))'*W;
Pmix = C*wm';
Emix = cumsum(Pmix, 1, "reverse") - Pmix;
excess_mix = max(Emix - E_env, [], "all");
% Transition windows
pct = [99 99.9 99.99];
q_env = pct_of_env(E_env, k0, dg, pct);
dur = cell(numel(ctx.scans), 1);
it = cell(numel(ctx.scans), 1);
for m = 1:numel(ctx.scans)
    s = ctx.scans{m};
    nr = ones(numel(s.tilt_deg), 1);
    if isfield(s, "n_rev") && ~isempty(s.n_rev)
        nr = s.n_rev(:).*ones(numel(s.tilt_deg), 1);
    end
    rate = s.rate_deg_s(:).*ones(numel(s.tilt_deg), 1);
    dur{m} = nr*360./abs(rate);
    [~, it{m}] = ismembertol(s.tilt_deg(:), ctx.T_n.tilt_deg, 1e-6, "DataScale", 1);
end
nT = numel(ctx.T_n.tilt_deg);
nm = numel(ctx.scans);
wy = zeros(nT, nm);
for y = 1:nm
    wy(:, y) = accumarray(it{y}, dur{y}, [nT 1]);
end
best = -Inf(1, numel(pct));
arg = strings(1, numel(pct));
for x = 1:nm
    for s = 1:numel(dur{x}) - 1
        wx = accumarray(it{x}(1:s), dur{x}(1:s), [nT 1]);
        Wwin = (wx + wy)./sum(wx + wy, 1);
        q = pct_of(C*Wwin, k0, dg, pct);
        d = q - q_env(:);
        for j = 1:numel(pct)
            [v, y] = max(d(j, :));
            if v > best(j)
                best(j) = v;
                arg(j) = sprintf("%s, %d cuts, then %s", ctx.names(x), s, ctx.names(y));
            end
        end
    end
end
T = table(["random mixtures (max exceedance excess)"; "transition window " + string(pct') + "%"], ...
    [excess_mix; best'], ["probability"; repmat("dB above envelope", numel(pct), 1)], ["-"; arg'], ...
    'VariableNames', ["check", "value", "unit", "worst_case"]);
writetable(T, evid(opts, "A3"));
key = sprintf("Mixtures: largest exceedance above the envelope %.2g (round-off). Mode change mid-volume, window of the partial volume plus the next: up to %s dB above the envelope at 99/99.9/99.99%%.", ...
    excess_mix, strjoin(compose("%.2f", best), "/"));
if excess_mix <= 1e-12 && all(best <= 1e-9)
    [res, sev, nxt] = deal("Pass", "None", "None");
elseif excess_mix <= 1e-12
    res = "Finding";
    sev = "Low";
    nxt = "The envelope bounds every mixture of whole volumes (proof: the mixture exceedance is a weighted mean of mode exceedances). A volume cut short by a mode change repeats the low cuts, so that window can exceed the envelope; it lasts one volume, so the long-run CDF moves little. State it in the report; per-volume harm (D4) is the metric that sees it.";
else
    res = "Fail";
    sev = "High";
    nxt = "Mixture bound violated: check the envelope code.";
end
r = row("A3", "The envelope bounds any mix of modes", ...
    "Random mode mixtures (Dirichlet, 5,000) and every mid-volume mode change (cuts 1..s of X, then Y), 0/0 deg", ...
    "Mixture bound holds; transitions add no exceedance above the envelope", res, key, rel(opts, evid(opts, "A3")), sev, nxt);
end

function q = pct_of_env(E, k0, dg, pct)
gain = (k0 + (0:numel(E) - 1)')*dg;
q = zeros(numel(pct), 1);
for j = 1:numel(pct)
    q(j) = gain(find(E <= 1 - pct(j)/100 + 1e-13, 1));
end
end

% -------------------------------------------------------------------------
function r = run_a4(ctx, opts)
% The mode table uses 0.5 deg as the lowest cut. Public NWS material lists
% sites with a supplemental cut below 0.5 deg (WDTD: 22 sites as of Build
% 19, -0.2 to +0.4 deg; SCN 24-09 and 26-20 add sites at 0.2 to 0.4 deg).
% Effect: add one supplemental cut to the worst mode.
low = [0.4 0.2 0.0 -0.2];
el_n = [-0.2 0 0.5];
pct = [99 99.9 99.99 99.999];
k = ctx.k_worst;
base = ctx.scans{k};
rows_ = {};
for e = el_n
    sn = ctx.sys_n;
    sn.scan = base;
    b = mag_exact_cdf_rev1(sn, ctx.sys_a, struct("el21_deg", e, "el12_deg", 0, "pct", pct));
    for a = low
        sn.scan = add_low_cut(base, a);
        x = mag_exact_cdf_rev1(sn, ctx.sys_a, struct("el21_deg", e, "el12_deg", 0, "pct", pct));
        d = x.pct21.gain_dBi' - b.pct21.gain_dBi';
        rows_(end + 1, :) = [{e, a}, num2cell(b.pct21.gain_dBi'), num2cell(d)]; %#ok<AGROW>
    end
end
vn = ["el_nexrad_deg", "supplemental_cut_deg", "base_" + pnames(pct), "delta_" + pnames(pct)];
T = cell2table(rows_, 'VariableNames', vn);
writetable(T, evid(opts, "A4"));
dmax = max(T{:, startsWith(vn, "delta_")}, [], "all");
key = sprintf("The %d modes use a 0.5 deg lowest cut; public NWS lists sites with a supplemental cut at -0.2 to +0.4 deg (4 sites at -0.2). Adding one such cut to the worst mode raises the MAG by up to %.2f dB (99 to 99.999%%, takeoff -0.2 to 0.5 deg).", numel(ctx.names), dmax);
r = row("A4", "The 74 modes are the complete set", ...
    "Compare with public NWS VCP material; add one supplemental low cut (-0.2 to 0.4 deg) to the worst mode", ...
    "Every operational VCP and lowest angle is in the table, or the gap is listed with its effect", ...
    "Finding", key, rel(opts, evid(opts, "A4")), "Medium", ...
    "Option added (low_el_deg and low_el_mode in nexrad_vcp_rev1, nexrad_asr11_mag_rev1 and nexrad_mode_scans_rev1). Run the affected sites with their lowest angle; site list: training.weather.gov/wdtd/courses/rac/principles/supplemental.php and NWS SCN 24-09, 26-20.");
end

function s = add_low_cut(s, a)
% Prepend one cut at a deg, copying the first (0.5 deg) cut's other fields.
n = numel(s.tilt_deg);
f = fieldnames(s);
for i = 1:numel(f)
    v = s.(f{i});
    if (isnumeric(v) || islogical(v)) && isvector(v) && numel(v) == n && n > 1
        s.(f{i}) = [v(1); v(:)];
        if isrow(v)
            s.(f{i}) = s.(f{i})';
        end
    end
end
s.tilt_deg(1) = a;
end

function s = pnames(pct)
s = "p" + replace(compose("%g", pct), ".", "_");
end

% =========================================================================
% B. Numerical accuracy
% =========================================================================
function r = run_b1(ctx, opts)
% Tables rebuilt with 0.005 dB bins and a 0.0002 deg azimuth step, at grid
% nodes (no interpolation): reported percentiles lie between the refined
% value (less its own bias) and the refined value plus the reported bias.
if opts.quick
    grid = [-0.5 0 0.5];
else
    grid = [-1 0 0.5 1 2 5];
end
Tn0 = gain_dist_table_rev1(ctx.sys_n.tx, ctx.scans, struct("name", "WSR-88D", "mode_names", ctx.names, "el_grid_deg", grid));
Tn1 = gain_dist_table_rev1(ctx.sys_n.tx, ctx.scans, struct("name", "WSR-88D", "mode_names", ctx.names, ...
    "el_grid_deg", grid, "bin_dB", 0.005, "az_step_deg", 0.0002));
Ta0 = gain_dist_table_rev1(ctx.sys_a.tx, ctx.sys_a.scan, struct("name", "ASR-11", "el_grid_deg", grid));
Ta1 = gain_dist_table_rev1(ctx.sys_a.tx, ctx.sys_a.scan, struct("name", "ASR-11", "el_grid_deg", grid, ...
    "bin_dB", 0.005, "az_step_deg", 0.0002));
[et, er] = ndgrid(grid, grid);
R0 = mag_pair_cdf_rev1(Tn0, Ta0, et(:), er(:), struct("pct", ctx.pct));
R1 = mag_pair_cdf_rev1(Tn1, Ta1, et(:), er(:), struct("pct", ctx.pct));
v0 = R0.pct{:, 3:2 + numel(ctx.pct)};
v1 = R1.pct{:, 3:2 + numel(ctx.pct)};
d = v0 - v1;
ok = (v0 >= v1 - 2*0.005 - 1e-9) & (v0 <= v1 + R0.bias_dB + 1e-9);
[ee, pp] = ndgrid(1:numel(et), ctx.pct);
T = table(et(ee(:)), er(ee(:)), pp(:), v0(:), v1(:), d(:), ok(:), 'VariableNames', ...
    ["el_tx_deg", "el_rx_deg", "pct", "reported_dBi", "refined_dBi", "diff_dB", "pass"]);
writetable(T, evid(opts, "B1"));
key = sprintf("%d node pairs x %d percentiles: reported minus refined %.3f to %.3f dB (stated bias %.2f dB); %d outside.", ...
    numel(et), numel(ctx.pct), min(d, [], "all"), max(d, [], "all"), R0.bias_dB, sum(~ok, "all"));
if all(ok, "all")
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    [res, sev, nxt] = deal("Fail", "High", "Find the cases in the evidence file and fix the rounding.");
end
r = row("B1", "Percentiles are conservative by at most bias_dB", ...
    "Tables rebuilt with 0.005 dB bins and 0.0002 deg azimuth step; WSR-88D tx (74 modes) to ASR-11 at grid nodes", ...
    "Reported value between the refined value and the refined value plus the bias, everywhere", res, key, rel(opts, evid(opts, "B1")), sev, nxt);
end

% -------------------------------------------------------------------------
function r = run_b2(ctx, opts)
% Elevation interpolation (quantile, table grid 0.05 deg) against exact
% slices, every 0.005 deg of the grid, for each side in turn; B3 from the
% grid nodes (single-precision table against double exact slices).
pct = [ctx.pct 99.9999];
g = ctx.T_n.el_grid_deg;
if opts.quick
    e = linspace(g(1), g(end), 21)';
else
    e = (g(1):0.005:g(end))';
end
e = round(e, 6);
Tq = cell(2, 1);
for side = 1:2
    if side == 1
        [et, er] = deal(e, zeros(size(e)));
    else
        [et, er] = deal(zeros(size(e)), e);
    end
    Rq = mag_pair_cdf_rev1(ctx.T_n, ctx.T_a, et, er, struct("pct", pct));
    Rx = mag_pair_cdf_rev1(ctx.T_n, ctx.T_a, et, er, struct("pct", pct, "interp", "exact"));
    vq = Rq.pct{:, 3:2 + numel(pct)};
    vx = Rx.pct{:, 3:2 + numel(pct)};
    node = abs(e/0.05 - round(e/0.05)) < 1e-6;
    [ee, pp] = ndgrid(1:numel(e), pct);
    Tq{side} = table(repmat(side, numel(ee), 1), e(ee(:)), node(ee(:)), pp(:), vq(:), vx(:), vq(:) - vx(:), ...
        'VariableNames', ["varied_side", "el_deg", "grid_node", "pct", "interp_dBi", "exact_dBi", "diff_dB"]);
end
T = vertcat(Tq{:});
T.varied_side = categorical(T.varied_side, [1 2], ["WSR-88D", "ASR-11"]);
writetable(T, evid(opts, "B2"));
off = ~T.grid_node & T.pct <= 99.999;
mx = max(abs(T.diff_dB(off)));
neg = mean(T.diff_dB(off) < -1e-9);
nd = T.grid_node & T.pct >= 99.999;
n_b3 = sum(abs(T.diff_dB(nd)) > 1e-9);
key = sprintf("Off-node, 50 to 99.999%%: max |interp - exact| %.2f dB; interp below exact in %.1f%% of checks. B3, at nodes, 99.999 and 99.9999%%: %d of %d differ.", ...
    mx, 100*neg, n_b3, sum(nd));
if mx <= 0.1 + 1e-9
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    res = "Fail";
    n_over = sum(abs(T.diff_dB(off)) > 0.1 + 1e-9);
    if mx <= 0.3
        sev = "Low";
    else
        sev = "Medium";
    end
    nxt = sprintf("%d of %d off-node checks exceed 0.1 dB (largest %.2f dB, inside the pair function's stated bias of %.2f dB); interpolation is one bin low in %.2f%% of checks. Use interp = ""exact"" (33 ms per pair) where 0.1 dB matters.", ...
        n_over, sum(off), mx, Rq.bias_dB, 100*neg);
end
r = row("B2", "Interpolation error is at most 0.1 dB", ...
    "Every 0.005 deg of the -5 to 10 deg grid, each side varied with the other at 0 deg, quantile interpolation vs exact slices, envelope of 74 modes", ...
    "Maximum error stated and inside 0.1 dB", res, key, rel(opts, evid(opts, "B2")), sev, nxt);
if n_b3 == 0
    [res3, sev3, nxt3] = deal("Pass", "None", "None");
else
    [res3, sev3, nxt3] = deal("Fail", "Low", "Store the table in double precision.");
end
r3 = row("B3", "Single-precision storage does not hurt the tail", ...
    "Grid nodes from B2: single-precision table vs double-precision exact slices at 99.999% and 99.9999%", ...
    "Same bin", res3, sprintf("%d of %d node checks differ.", n_b3, sum(nd)), rel(opts, evid(opts, "B2")), sev3, nxt3);
r3.runtime_s = 0;
r3.quick = opts.quick;
writetable(struct2table(r3, "AsArray", true), fullfile(opts.out_dir, "B3_row.csv"));
end

% -------------------------------------------------------------------------
function r = run_b4(ctx, opts)
% Closed-form slice against a fine azimuth grid at the edge cases.
S_cf = gain_slice_prep_rev1(ctx.sys_n.tx);
S_gr = gain_slice_prep_rev1(ctx.sys_n.tx, struct("method", "grid", "az_step_deg", 0.0005));
dg = S_cf.bin_dB;
cases = [0.5 0.5; 19.5 19.5; 0.5 -3; 0.5 -0.5; 19.5 0; 0.5 10; 19.5 -5; 0 0; 4.0 4.0; 0.5 0.4999];
lev = [0 50 90 99 99.9 99.99 99.999 100]/100;
rows_ = zeros(size(cases, 1), 5);
for i = 1:size(cases, 1)
    h1 = gain_slice_rev1(S_cf, cases(i, 1), cases(i, 2));
    h2 = gain_slice_rev1(S_gr, cases(i, 1), cases(i, 2));
    q1 = slice_q(h1, lev, dg);
    q2 = slice_q(h2, lev, dg);
    rows_(i, :) = [cases(i, :), min(q1 - q2), max(q1 - q2), h1.gmax - h2.gmax];
end
T = array2table(rows_, 'VariableNames', ["tilt_deg", "el_deg", "min_cf_minus_grid_dB", "max_cf_minus_grid_dB", "gmax_cf_minus_grid_dB"]);
writetable(T, evid(opts, "B4"));
ok = T.min_cf_minus_grid_dB >= -1e-9 & T.max_cf_minus_grid_dB <= dg + 1e-9;
key = sprintf("%d edge cases (el = tilt, negative el, tilt 19.5 deg): closed form minus grid %.2f to %.2f dB over 0 to 100%%; %d outside [0, one bin].", ...
    height(T), min(T.min_cf_minus_grid_dB), max(T.max_cf_minus_grid_dB), sum(~ok));
if all(ok)
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    [res, sev, nxt] = deal("Fail", "High", "Closed form low or off by more than one bin at the cases listed.");
end
r = row("B4", "The closed-form slice (T1) is exact", ...
    "WSR-88D closed form vs a 0.0005 deg azimuth grid at el = tilt, negative el and tilt 19.5 deg", ...
    "Within one bin, never low", res, key, rel(opts, evid(opts, "B4")), sev, nxt);
end

function q = slice_q(h, lev, dg)
F = cumsum(h.p);
g = (h.k_lo + (0:numel(h.p) - 1)')*dg;
q = zeros(size(lev));
for j = 1:numel(lev)
    q(j) = g(find(F >= lev(j) - 1e-12, 1));
end
end

% -------------------------------------------------------------------------
function r = run_b5(ctx, opts)
% 20 seeds of the original model's baseline (random mode, 5 million
% samples): does each stated 95% interval cover the exact value?
n_seed = 20;
ns = 5e6;
if opts.quick
    n_seed = 3;
    ns = 1e5;
end
k = ctx.k_worst;
ex = nexrad_asr11_mag_rev1(mode_opts(ctx, k, struct("seed", 1, "n_samples", 1e3, "duration_s", 60, ...
    "n_trials", 1, "cdf_method", "exact")));
xv = ex.cdf.pct21.gain_dBi;
b = ex.cdf.bias_dB;
cov = false(numel(xv), n_seed);
lo = NaN(numel(xv), n_seed);
hi = lo;
for s = 1:n_seed
    o = nexrad_asr11_mag_rev1(mode_opts(ctx, k, struct("seed", s, "n_samples", ns, "duration_s", 60, "n_trials", 1)));
    lo(:, s) = o.cdf.pct21.ci95_lo_dBi;
    hi(:, s) = o.cdf.pct21.ci95_hi_dBi;
    cov(:, s) = hi(:, s) >= xv - b - 1e-9 & lo(:, s) <= xv + 1e-9;
end
T = table(ex.cdf.pct21.pct, xv, sum(cov, 2), repmat(n_seed, numel(xv), 1), median(hi - lo, 2), ...
    'VariableNames', ["pct", "exact_dBi", "seeds_covered", "seeds", "median_ci_width_dB"]);
writetable(T, evid(opts, "B5"));
lo_ok = binoinv(0.025, n_seed, 0.95);           % fewest covered seeds consistent with 95%
ok = T.seeds_covered >= lo_ok;
key = sprintf("Seeds covering the exact value (of %d), by percentile 50 to 99.999%%: %s (at least %d expected).", ...
    n_seed, strjoin(string(T.seeds_covered'), ", "), lo_ok);
if all(ok)
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    [res, sev, nxt] = deal("Fail", "Medium", "Intervals too narrow at the listed percentiles: widen them or report the exact CDF.");
end
r = row("B5", "Monte Carlo intervals in the original model are honest", ...
    sprintf("Worst mode at 0/0, %d seeds of %g samples; each 95%% interval vs the exact value", n_seed, ns), ...
    "The stated 95% interval covers about 95% of seeds", res, key, rel(opts, evid(opts, "B5")), sev, nxt);
end

% =========================================================================
% C. Inputs
% =========================================================================
function r = run_c1(ctx, opts)
% WSR-88D pattern: NSSL KOUN measured (baseline) vs ROC specification and
% F.1245, on a random pair set.
pats = ["nssl_measured", "roc_spec", "f1245"];
n_pair = 1000;
if opts.quick
    n_pair = 10;
end
rs = RandStream("mrg32k3a", "Seed", 8);
g = ctx.T_n.el_grid_deg;
et = g(1) + (g(end) - g(1))*rand(rs, n_pair, 1);
er = g(1) + (g(end) - g(1))*rand(rs, n_pair, 1);
V = cell(numel(pats), 1);
for i = 1:numel(pats)
    o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
        "cdf_method", "exact", "nexrad_pattern", pats(i)));
    T = gain_dist_table_rev1(o.sys_nexrad.tx, ctx.scans, struct("name", "WSR-88D " + pats(i), ...
        "mode_names", ctx.names, "el_grid_deg", g));
    R = mag_pair_cdf_rev1(T, ctx.T_a, et, er, struct("pct", ctx.pct));
    V{i} = R.pct{:, 3:2 + numel(ctx.pct)};
end
rows_ = {};
for i = 2:numel(pats)
    d = V{i} - V{1};
    for j = 1:numel(ctx.pct)
        rows_(end + 1, :) = {pats(i), ctx.pct(j), median(d(:, j)), min(d(:, j)), max(d(:, j)), ...
            median(V{1}(:, j))}; %#ok<AGROW>
    end
end
T = cell2table(rows_, 'VariableNames', ["pattern", "pct", "median_delta_vs_nssl_dB", ...
    "min_delta_dB", "max_delta_dB", "median_nssl_dBi"]);
writetable(T, evid(opts, "C1"));
s = T(T.pattern == "roc_spec", :);
key = sprintf("%d pairs, 74-mode envelope. ROC specification minus NSSL measured, median by percentile 50 to 99.999%%: %s dB (range %.1f to %.1f). F.1245 median: %s dB.", ...
    n_pair, strjoin(compose("%.1f", s.median_delta_vs_nssl_dB'), "/"), min(s.min_delta_dB), max(s.max_delta_dB), ...
    strjoin(compose("%.1f", T.median_delta_vs_nssl_dB(T.pattern == "f1245")'), "/"));
r = row("C1", "The NSSL KOUN measured pattern represents the fleet", ...
    "ROC specification and F.1245 WSR-88D tables, same 1,000 random pairs, 74-mode envelope", ...
    "Difference reported for every percentile; the NOAA case uses the specification", "Reported", key, ...
    rel(opts, evid(opts, "C1")), "Medium", "Closed 2026-11-10: roc_spec is now the default pattern (decision 7); the measured pattern stays as a sensitivity case. A fleet sample of measured patterns is NOAA data, not public.");
end

% -------------------------------------------------------------------------
function r = run_c2(ctx, opts)
% ASR-9 remote sidelobe: +0.5, +3.5 (baseline, borrowed from M.1464-2
% Radar C), +6.5 dBi.
sl = [0.5 3.5 6.5];
pct = [50 90 99 99.9 99.99 99.999];
V = zeros(numel(sl), numel(pct));
for i = 1:numel(sl)
    o = mode_opts(ctx, ctx.k_worst, struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
        "cdf_method", "exact", "atc_radar", "asr9", "asr_remote_sl_dBi", sl(i), "pct", pct));
    out = nexrad_asr11_mag_rev1(o);
    V(i, :) = out.cdf.pct21.gain_dBi';
end
slope = (V(3, :) - V(1, :))/(sl(3) - sl(1));
T = array2table([sl(:), V], 'VariableNames', ["remote_sl_dBi", pnames(pct)]);
T = [T; array2table([NaN, slope], 'VariableNames', T.Properties.VariableNames)];
T.row = ["sl " + string(sl(:)); "slope dB per dB"];
writetable(T, evid(opts, "C2"));
key = sprintf("Slope of the MAG per dB of remote sidelobe, 50 to 99.999%%: %s. Baseline (+3.5 dBi) 99.99%%: %.1f dBi.", ...
    strjoin(compose("%.2f", slope), "/"), V(2, 5));
one = abs(slope - 1) <= 0.1;
if any(one)
    res = "Reported (one-for-one at " + strjoin(compose("%g", pct(one)), ", ") + "%)";
else
    res = "Reported (not one-for-one at any percentile)";
end
r = row("C2", "The ASR-9 remote sidelobe (+3.5 dBi, borrowed) is safe", ...
    "ASR-9, worst mode, 0/0 deg, remote sidelobe +0.5, +3.5, +6.5 dBi (exact CDF)", ...
    "One-for-one sensitivity confirmed; value stays TODO until a public source is found", res, key, ...
    rel(opts, evid(opts, "C2")), "Medium", "Closed 2026-11-10: public bound found (WEBER_1986, MIT LL ATC-138: maximum azimuth sidelobes -28 dB), so the ASR-9 default is now +5.5 dBi (33.5 dBi gain - 28 dB); the NOAA case keeps +6.5 dBi.");
end

% -------------------------------------------------------------------------
function r = run_c3(ctx, opts)
% Angle convention: the C# tool (CouplingCalculations.cs,
% CalculateWorstCaseCoupling_dBi) uses offset = takeoff - tilt for an
% axially symmetric pattern, with takeoff angles from the propagation model
% (SpreadsheetIO.cs "Takeoff Angle(deg)" columns). On the path (beam
% azimuth on the other radar), the MATLAB off-axis angle reduces to
% |el - tilt|. Check that identity on three hand cases.
pn = antenna_prepare_rev1(ctx.sys_n.tx);
cases = [0.3 0.5; -0.4 0.5; 2.0 1.5];          % [takeoff, tilt]
rows_ = zeros(size(cases, 1), 4);
a = ctx.sys_n.tx;
for i = 1:size(cases, 1)
    st = struct("mech_az", 0, "tilt", cases(i, 2), "steer_az", 0, "steer_el", 0);
    g_m = antenna_gain_rev1(pn, 0, cases(i, 1), st);
    off = abs(cases(i, 1) - cases(i, 2));
    g_c = a.g_max_dBi + interp1(a.offaxis_deg, a.gain_dB, off);   % C#: G(takeoff - tilt), symmetric
    rows_(i, :) = [cases(i, :), g_m, g_c];
end
T = array2table(rows_, 'VariableNames', ["takeoff_deg", "tilt_deg", "matlab_gain_dBi", "csharp_rule_gain_dBi"]);
T.diff_dB = T.matlab_gain_dBi - T.csharp_rule_gain_dBi;
writetable(T, evid(opts, "C3"));
key = sprintf("Convention matches: both take the takeoff angle from the local horizontal, positive up, and use takeoff - tilt; on three hand cases MATLAB minus the C# rule is at most %.1e dB. Earth curvature sits in the propagation model's takeoff angle; MATLAB takes the angle as given. New: CalculateWorstCaseCoupling_dBi passes otherAntModel1 for antenna 2 (should be otherAntModel2).", ...
    max(abs(T.diff_dB)));
r = row("C3", "Takeoff angles from the RFS tool mean what the model assumes", ...
    "Code reading of CouplingCalculations.cs and SpreadsheetIO.cs; MATLAB on-path gain vs the C# offset rule on three hand cases", ...
    "Same angle within 0.01 deg", "Closed (needs NOAA/FAA data)", key, rel(opts, evid(opts, "C3")), "Low", ...
    "Closed 2026-11-10 as out of scope here. Method for whoever has the RFS workbook: take three pairs, recompute each takeoff angle from the profile (4/3 earth, antenna heights) and check it against the tool's ""Takeoff Angle(deg)"" columns within 0.01 deg; then feed those angles to mag_pair_cdf_rev1. Report the otherAntModel1 bug (X4) to the tool owner.");
end

% -------------------------------------------------------------------------
function r = run_c4(ctx, opts)
tilt = -1:0.5:2.5;                               % wrapper limit: M.1464-2 Radar C range, -2.5 to 2.5 deg
pct = [50 90 99 99.9 99.99 99.999];
V = zeros(numel(tilt), numel(pct));
for i = 1:numel(tilt)
    o = mode_opts(ctx, ctx.k_worst, struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
        "cdf_method", "exact", "asr_tilt_deg", tilt(i), "pct", pct));
    out = nexrad_asr11_mag_rev1(o);
    V(i, :) = out.cdf.pct21.gain_dBi';
end
T = array2table([tilt(:), V], 'VariableNames', ["asr_tilt_deg", pnames(pct)]);
writetable(T, evid(opts, "C4"));
rng_ = max(V) - min(V);
[~, iw] = max(V(:, 5));
key = sprintf("ASR-11 tilt -1 to +2.5 deg, worst mode, 0/0: MAG range by percentile 50 to 99.999%%: %s dB. Largest 99.99%% at tilt %g deg (%.1f dBi; 0 deg gives %.1f).", ...
    strjoin(compose("%.1f", rng_), "/"), tilt(iw), V(iw, 5), V(tilt == 0, 5));
r = row("C4", "Site tilt does not matter", "ASR-11 beam tilt -1 to +2.5 deg in 0.5 deg steps (the M.1464-2 Radar C range caps the plan's +3 deg), worst mode, 0/0 deg (exact CDF)", ...
    "Effect stated; tilt becomes a per-site input", "Finding", key, rel(opts, evid(opts, "C4")), "Medium", ...
    "Tilt is already a per-pair input (asr_tilt_deg; gain_dist_table_rev1 per tilt). Site tilts are FAA data, not public: the NOAA case uses the worst tilt.");
end

% -------------------------------------------------------------------------
function r = run_c5(ctx, opts) %#ok<INUSD>
% Far-field distance 2*D^2/lambda. WSR-88D: 8.53 m (28 ft) reflector
% [ROC_TECH]. ASR-11 and ASR-9: horizontal aperture estimated from the
% 1.4 deg azimuth beamwidth, D = 70*lambda/theta (estimate, no public
% dimension found).
f = [2.7 2.8 2.9 3.0]*1e9;
lam = physconst("LightSpeed")./f;
nm = ["WSR-88D (8.53 m reflector)"; "ASR-11/ASR-9 (aperture from 1.4 deg beamwidth)"];
rows_ = {};
for i = 1:numel(f)
    D = [8.53; 70*lam(i)/1.4];
    for a = 1:2
        rows_(end + 1, :) = {nm(a), f(i)/1e9, D(a), 2*D(a)^2/lam(i)/1e3}; %#ok<AGROW>
    end
    rows_(end + 1, :) = {"pair, 2(D1+D2)^2/lambda", f(i)/1e9, sum(D), 2*sum(D)^2/lam(i)/1e3}; %#ok<AGROW>
end
T = cell2table(rows_, 'VariableNames', ["antenna", "freq_GHz", "aperture_m", "far_field_km"]);
writetable(T, evid(opts, "C5"));
key = sprintf("Far field: WSR-88D %.1f to %.1f km (2.7 to 3.0 GHz); ASR %.2f km; strict pair criterion up to %.1f km. Pairs closer than %.1f km need near-field coupling.", ...
    min(T.far_field_km(startsWith(T.antenna, "WSR"))), max(T.far_field_km(startsWith(T.antenna, "WSR"))), ...
    max(T.far_field_km(startsWith(T.antenna, "ASR"))), max(T.far_field_km(startsWith(T.antenna, "pair"))), ...
    max(T.far_field_km(startsWith(T.antenna, "WSR"))));
r = row("C5", "Far-field patterns apply", "2D^2/lambda for each antenna at 2.7 to 3.0 GHz, against the shortest pair distance in the study", ...
    "Every pair beyond the far-field distance, or flagged", "Closed (needs NOAA/FAA data)", key, rel(opts, evid(opts, "C5")), "Low", ...
    "Closed 2026-11-10 as out of scope here. Rule for the study: flag every pair closer than 1.5 km (WSR-88D far field at 3.0 GHz; 3.7 km under the strict two-antenna criterion) from the RFS workbook's ""Ground Distance"" column and treat it with a near-field method. The ASR aperture is an estimate.");
end

% =========================================================================
% D. Time series
% =========================================================================
function r = run_d2(ctx, opts)
% Bootstrap over hours (trials) of the 120 h persistence statistics.
n_tr = 120;
dur = 3600;
n_b = 2000;
if opts.quick
    [n_tr, dur, n_b] = deal(3, 300, 50);
end
o = mode_opts(ctx, ctx.k_worst, struct("seed", 1, "n_samples", 1e3, "n_trials", n_tr, "duration_s", dur, ...
    "cdf_method", "exact", "ts_method", "events"));
out = nexrad_asr11_mag_rev1(o);
thr = out.opts.persist_thr_dBi;
A = out.scan_max.asr_rev;
V = out.scan_max.nexrad_volume;
A = A(A.complete, :);
V = V(V.complete, :);
nt = n_tr;
[ra, na, rv, nv, lr] = deal(zeros(nt, numel(thr)));
for j = 1:nt
    a = A.max_g21_dBi(A.trial == j);
    v = V.max_g21_dBi(V.trial == j);
    for k = 1:numel(thr)
        hit = a >= thr(k);                      % NaN (censored, below floor) counts as no hit
        ra(j, k) = sum(hit);
        na(j, k) = numel(a);
        rv(j, k) = sum(v >= thr(k));
        nv(j, k) = numel(v);
        lr(j, k) = longest_run(hit);
    end
end
rs = RandStream("mrg32k3a", "Seed", 21);
B = randi(rs, nt, nt, n_b);
stat = @(idx) [sum(ra(idx, :), 1)./sum(na(idx, :), 1); sum(rv(idx, :), 1)./sum(nv(idx, :), 1); max(lr(idx, :), [], 1)];
S = zeros(3, numel(thr), n_b);
for b = 1:n_b
    S(:, :, b) = stat(B(:, b));
end
est = stat((1:nt)');
lo = prctile(S, 2.5, 3);
hi = prctile(S, 97.5, 3);
nmz = ["rev_reach_share", "vol_reach_share", "longest_run_revs"];
rows_ = {};
for s = 1:3
    for k = 1:numel(thr)
        rows_(end + 1, :) = {nmz(s), thr(k), est(s, k), lo(s, k), hi(s, k), hi(s, k)/lo(s, k)}; %#ok<AGROW>
    end
end
T = cell2table(rows_, 'VariableNames', ["statistic", "thr_dBi", "estimate", "ci95_lo", "ci95_hi", "hi_over_lo"]);
writetable(T, evid(opts, "D2"));
vol = T(T.statistic == "vol_reach_share", :);
wide = vol.thr_dBi(~(vol.hi_over_lo <= 2));
if isempty(wide)
    wtxt = "no threshold";
else
    wtxt = strjoin(compose("%g", wide'), ", ") + " dBi";
end
key = sprintf("%d h, worst mode, 0/0: volume share reaching 70 dBi %.1f%% (95%% %.1f to %.1f%%). Lock-probability interval wider than a factor of two at %s.", ...
    nt*dur/3600, 100*vol.estimate(vol.thr_dBi == 70), 100*vol.ci95_lo(vol.thr_dBi == 70), ...
    100*vol.ci95_hi(vol.thr_dBi == 70), wtxt);
if isempty(wide)
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    res = "Finding";
    sev = "Low";
    nxt = "Report the interval at the thresholds listed, or extend the run (the events method costs about 18 ms per hour).";
end
r = row("D2", "120 hours is enough", "Bootstrap over 120 one-hour trials (2,000 resamples) of revolution and volume reach shares and the longest run", ...
    "Intervals stated; extend the run if the lock probability interval is wider than a factor of two", res, key, ...
    rel(opts, evid(opts, "D2")), sev, nxt);
end

function n = longest_run(hit)
d = diff([0; hit(:); 0]);
s = find(d == 1);
e = find(d == -1);
n = max([0; e - s]);
end

% =========================================================================
% F. Software
% =========================================================================
function r = run_f1(ctx, opts)
% Property test: random patterns, tilts, scan rates, sectors and takeoff
% angles; the exact CDF must fall inside the Monte Carlo order-statistic
% band. Bonferroni over all checks keeps the chance of any false failure
% near 1%.
n_case = 500;
ns = 1e6;
if opts.quick
    n_case = 3;
    ns = 2e4;
end
pct = [50 90 99 99.9];
alpha = 0.01/(n_case*numel(pct)*2);
rs = RandStream("mrg32k3a", "Seed", 31);
res_ = zeros(n_case, 6);
desc = strings(n_case, 1);
w = warning("off", "all");
cleanup = onCleanup(@() warning(w));
for c = 1:n_case
    [s1, d1] = random_sys(rs, ctx);
    [s2, d2] = random_sys(rs, ctx);
    e12 = -3 + 9*rand(rs);
    e21 = -3 + 9*rand(rs);
    desc(c) = d1 + " | " + d2;
    try
        x = mag_exact_cdf_rev1(s1, s2, struct("el12_deg", e12, "el21_deg", e21, "pct", pct));
        m = radar_mutual_gain_cdf_rev1(s1, s2, struct("seed", c, "mode", "random", "n_samples", ns, ...
            "el12_deg", e12, "el21_deg", e21, "pct", pct));
        [ok12, z12] = in_band(sort(m.g12_dBi), x.pct12.gain_dBi, pct, alpha, x.bias_dB);
        [ok21, z21] = in_band(sort(m.g21_dBi), x.pct21.gain_dBi, pct, alpha, x.bias_dB);
        res_(c, :) = [c, e12, e21, ok12 && ok21, max([z12 z21]), x.bias_dB];
    catch err
        % A crash is a failure; keep the systems so it can be reproduced.
        res_(c, :) = [c, e12, e21, 0, NaN, NaN];
        desc(c) = desc(c) + " | ERROR " + err.identifier + ": " + err.message;
        crash.(sprintf("case%d", c)) = struct("s1", s1, "s2", s2, "e12", e12, "e21", e21, "err", err.message);
    end
end
if exist("crash", "var")
    save(fullfile(opts.out_dir, "F1_crash_cases.mat"), "crash");
end
T = array2table(res_, 'VariableNames', ["case", "el12_deg", "el21_deg", "pass", "worst_excess_dB", "bias_dB"]);
T.systems = desc;
writetable(T, evid(opts, "F1"));
n_fail = sum(~T.pass);
n_err = sum(contains(T.systems, "| ERROR"));
key = sprintf("%d crashed. ", n_err) + sprintf("%d random cases (statgain, F.1245, random off-axis tables, NSSL WSR-88D, ASR-11 CSC^2; 1 to 4 tilts; 20%% sectors), %g samples each: %d outside the band (alpha %.1e per check).", ...
    n_case, ns, n_fail, alpha);
if n_fail == 0
    [res, sev, nxt] = deal("Pass", "None", "None");
else
    [res, sev, nxt] = deal("Fail", "High", "Failing cases in the evidence file.");
end
r = row("F1", "The two models agree", "Property test: random patterns, tilts, rates, sectors and takeoff angles; exact vs Monte Carlo order-statistic band", ...
    "500 random cases, no failure", res, key, rel(opts, evid(opts, "F1")), sev, nxt);
end

function [ok, worst] = in_band(s, q, pct, alpha, bias)
% Exact percentile q (conservative by at most bias) inside the two-sided
% order-statistic band of the sorted samples s.
n = numel(s);
ok = true;
worst = -Inf;
for j = 1:numel(pct)
    p = pct(j)/100;
    i_lo = max(1, binoinv(alpha/2, n, p));
    i_hi = min(n, binoinv(1 - alpha/2, n, p) + 1);
    lo = s(i_lo);
    hi = s(i_hi);
    ex = max(lo - q(j), q(j) - bias - hi);       % > 0: outside
    worst = max(worst, ex);
    ok = ok && ex <= 1e-9;
end
end

function [s, d] = random_sys(rs, ctx)
% A random mechanically rotating radar.
u = rand(rs);
if u < 0.2
    ant = struct("type", "statgain", "g_max_dBi", 15 + 31*rand(rs));
    d = sprintf("statgain %.1f", ant.g_max_dBi);
elseif u < 0.4
    ant = struct("type", "f1245", "g_max_dBi", 25 + 21*rand(rs));
    d = sprintf("f1245 %.1f", ant.g_max_dBi);
elseif u < 0.7
    bw = 0.8 + 4*rand(rs);
    x = sort([0, bw/2, bw*(1 + 3*rand(rs, 1, 4)), 90, 180]);
    g = [0, -3, -3 - cumsum(5 + 15*rand(rs, 1, 4)), -50 - 10*rand(rs), -60 - 5*rand(rs)];
    g = min(g, [0 -3 cummin(g(3:end))]);
    ant = struct("type", "offaxis_table", "g_max_dBi", 20 + 26*rand(rs), "offaxis_deg", x, "gain_dB", g);
    d = sprintf("table %.1f bw %.2f", ant.g_max_dBi, bw);
elseif u < 0.85
    ant = ctx.sys_n.tx;
    d = "NSSL WSR-88D";
else
    ant = ctx.sys_a.tx;
    d = "ASR-11 CSC2";
end
n = randi(rs, 4);
sc = struct("type", "mech_rotating", "rate_deg_s", 5 + 85*rand(rs, n, 1), ...
    "tilt_deg", -1 + 9*rand(rs, n, 1), "n_rev", ones(n, 1));
if string(ant.type) == "rect"
    sc.tilt_deg = 3*rand(rs, n, 1);
end
if rand(rs) < 0.2
    sc.sector_deg = 30 + 270*rand(rs);
    sc.sector_center_deg = 360*rand(rs) - 180;
    if rand(rs) < 0.5
        sc.sweep_mode = "bidirectional";
        sc.turnaround_s = rand(rs);
    end
    d = d + sprintf(", sector %.0f", sc.sector_deg);
end
s = struct("name", d, "tx", ant, "rx", ant, "scan", sc);
d = string(d) + sprintf(", %d tilts", n);
end

% =========================================================================
% NOAA case
% =========================================================================
function r = run_noaa(ctx, opts)
% Baseline next to a "NOAA case": the inputs a pessimistic NOAA radar
% engineer would pick from the C attacks: ROC specification pattern (C1),
% ASR remote sidelobe +6.5 dBi (C2), the worst ASR tilt in -1 to +2.5 deg
% (C4, -1 to +2.5 deg), per-trial rate offset (the most persistent drift model, D1),
% worst mode, 0/0 deg.
n_tr = 120;
if opts.quick
    n_tr = 2;
end
tilt = -1:0.5:2.5;                               % wrapper limit: M.1464-2 Radar C range, -2.5 to 2.5 deg
p = zeros(size(tilt));
for i = 1:numel(tilt)
    o = mode_opts(ctx, ctx.k_worst, struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
        "cdf_method", "exact", "asr_tilt_deg", tilt(i), "nexrad_pattern", "roc_spec", "asr_remote_sl_dBi", 6.5));
    out = nexrad_asr11_mag_rev1(o);
    p(i) = out.cdf.pct21.gain_dBi(out.cdf.pct21.pct == 99.99);
end
[~, iw] = max(p);
common = struct("seed", 1, "n_samples", 1e3, "n_trials", n_tr, "duration_s", 3600, ...
    "cdf_method", "exact", "ts_method", "events");
if opts.quick
    common.duration_s = 300;
end
cfg = {struct(), struct("nexrad_pattern", "roc_spec", "asr_remote_sl_dBi", 6.5, "asr_tilt_deg", tilt(iw))};
lab = ["baseline", "noaa_case"];
M = cell(2, 1);
for c = 1:2
    o = mode_opts(ctx, ctx.k_worst, common);
    f = fieldnames(cfg{c});
    for i = 1:numel(f)
        o.(f{i}) = cfg{c}.(f{i});
    end
    out = nexrad_asr11_mag_rev1(o);
    P = out.scan_max.persistence;
    PV = out.scan_max.persistence_volume;
    met = ["MAG " + string(out.cdf.pct21.pct) + "% (dBi)"; ...
        "revolutions reaching " + string(P.threshold_dBi) + " dBi (%)"; ...
        "volumes reaching " + string(PV.threshold_dBi) + " dBi (%)"];
    val = [out.cdf.pct21.gain_dBi; P{:, 2}; PV{:, 2}];
    M{c} = table(met, val, 'VariableNames', ["metric", lab(c)]);
end
T = join(M{1}, M{2});
T.delta = T.noaa_case - T.baseline;
writetable(T, evid(opts, "NOAA"));
k = find(T.metric == "MAG 99.99% (dBi)");
kv = find(T.metric == "volumes reaching 70 dBi (%)");
key = sprintf("NOAA case (ROC spec pattern, ASR remote sidelobe +6.5 dBi, ASR tilt %g deg, worst mode, 0/0): 99.99%% MAG %.1f dBi vs %.1f baseline; volumes reaching 70 dBi %.1f%% vs %.1f%%.", ...
    tilt(iw), T.noaa_case(k), T.baseline(k), T.noaa_case(kv), T.baseline(kv));
r = row("NOAA", "NOAA case run beside the baseline", ...
    sprintf("Exact CDF and %d h event time series, baseline vs NOAA-case inputs", n_tr), ...
    "Both reported side by side in the Rev 6 report", "Reported", key, rel(opts, evid(opts, "NOAA")), "Info", ...
    "Carry both columns into the Rev 6 report; NOAA and FAA reviewers set which inputs stand.");
end

% =========================================================================
% Register
% =========================================================================
function T = make_register(opts)
f = dir(fullfile(opts.out_dir, "*_row.csv"));
R = cell(numel(f), 1);
for i = 1:numel(f)
    R{i} = readtable(fullfile(f(i).folder, f(i).name), "TextType", "string", "Delimiter", ",");
    R{i} = R{i}(:, ["id", "claim", "test", "criterion", "result", "key_result", "evidence", "severity", "next_step", "quick"]);
end
S = static_rows();
T = [vertcat(R{:}); S];
if any(T.quick == 1)
    warning("redteam:quickRows", "Register holds quick (code-path) rows: %s.", strjoin(T.id(T.quick == 1), ", "));
end
order = ["A1" "A2" "A3" "A4" "B1" "B2" "B3" "B4" "B5" "C1" "C2" "C3" "C4" "C5" "D1" "D2" "D3" "D4" ...
    "E1" "E2" "F1" "F2" "F3" "NOAA" "X1" "X2" "X3" "X4" "L2" "R1" "R2" "R3" "R4" "R5" "R6" "G1" "G2" "G3" "G4"];
[~, o] = ismember(T.id, order);
o(o == 0) = numel(order) + 1;
[~, s] = sort(o);
T = T(s, :);
T.quick = [];
writetable(T, fullfile(opts.repo_dir, "docs", "redteam_register_rev1.csv"));
end

function T = static_rows()
% Attacks closed in earlier steps, blocked attacks and scope limits.
c = {
"D1", "Independent rate jitter is realistic", "Persistence tables rerun with correlated drift (random walk) and per-revolution draws (step 6)", "Lock probability and longest run reported for both; the larger one is used", "Finding", "Volumes reaching 70 dBi: 14.2% per-trial offset (baseline), 13.9% random walk, 10.6% per-revolution draws; the per-trial baseline is the most persistent and stays the reference. The spread between drift models is about the size of D2's sampling interval (10.6 to 14.2%), so 120 h separates them only roughly.", "docs/figures/drift_persistence_rev1.csv", "Low", "None; the drift step size has no public source (TODO)."
"D3", "The event method misses nothing", "Section 5.4 tests (test_timeseries_rev1): event-window maxima vs the fixed grid; floor check", "No maximum above the floor outside a window", "Pass", "Maxima above the 15.4 dBi floor equal the fixed grid; reproduces the Rev 5 persistence tables exactly.", "test_timeseries_rev1.m", "None", "None"
"D4", "Time percentiles describe harm", "Corrupted 0.5 deg radials per cut (step 6), checked against a brute-force count", "Harm metric reported beside the percentiles", "Pass", "Harm metric implemented and tested; to be reported per pair in Rev 6.", "test_timeseries_rev1.m", "None", "Report it beside the percentiles in Rev 6."
"E1", "Centered on the path is the worst case", "Bearing swept across and beyond the sector (step 5)", "Maximum at the center, or the true worst bearing reported", "Finding", "Holds for a restart sweep (99.99%: 69.95 dBi centered, 64.40 at the edge); fails for a bidirectional sweep with turnaround: the edge is 7.4 dB worse than the center at 99.99%.", "docs/figures/sector_bearing_sweep_rev1.csv", "Medium", "Sector runs sweep the bearing; never assume the center."
"E2", "Restarting sweep equals back-and-forth sweep", "Both sweep laws, with turnaround, center and edge paths (step 5)", "Difference stated; a path at the sector edge is tested", "Finding", "They differ most at the sector edge (see E1); the sweep law is an input with no public operational value (TODO).", "docs/figures/sector_bearing_sweep_rev1.csv", "Medium", "Ask NOAA for the sector sweep law if sector scans are proposed."
"F2", "The tests would catch a bug", "Ten planted faults across both models (redteam_fault_seed_rev1)", "Every fault fails at least one test", "Pass", "10 of 10 caught (rerun after the review fixes). Fault 6 (swapped interpolation weights), once caught only by the pinned regression reference, now also fails testInterpolationWeightsPointTheRightWay and testUpperAndRoundingNeverBelowQuantile.", "docs/figures/fault_seed_rev1.csv", "None", "None"
"F3", "The model agrees with the C# tool", "Extend crosscheck_csharp_rev1 to the tool's gain-distribution and MACD outputs for one pair with nonzero takeoff angles", "Same percentiles within the tool's bin width", "Closed (needs NOAA/FAA data)", "Closed 2026-11-10 as out of scope here: needs an RFS workbook with takeoff-angle and coupling-distribution columns (the project workbook is over the 20 MB in-session limit).", "-", "Medium", "Method for whoever has the workbook: pick one pair with nonzero takeoff angles, run the tool's ScanRateAnalysis and mag_exact_cdf_rev1 with the same patterns, rates and angles, and compare percentiles within the tool's bin width (crosscheck_csharp_rev1 already matches the antenna patterns)."
"X1", "The wrapper handles every ATC geometry", "Found by the NOAA run: ASR tilt 2.5 deg with a +6.5 dBi floor", "Runs without error", "Fixed", "The ATC gain toward the path was flat (floor at every azimuth), so the lobe-region edges were empty and nexrad_asr11_mag_rev1 crashed. Fixed: no main-beam region in that case.", "test_nexrad_asr11_rev1.m (testFlatAtcCutRuns)", "Medium", "None"
"X2", "Closed-form slices handle any sector", "Found by F1: 10 of 500 random cases crashed (sector inside one segment of a coarse table)", "Runs and matches the grid", "Fixed", "refine_cells built a row instead of a column when a slice spanned one pattern cell (repelem of a scalar). Fixed; the F1 rerun has no crash.", "test_gain_slice_rev1.m (testSectorInsideOnePatternSegment)", "High", "None"
"X3", "The stated grid bias bounds the exact CDF error", "Found by F1: 3 cases above the old 2-bin bound by 0.004 to 0.011 dB", "Error inside the stated bias", "Fixed", "Closed-form slices round up twice (each cell's larger end gain, then the bin edge), so the bound is 2 bins per antenna: mag_exact_cdf_rev1 bias_dB 0.2 dB (was 0.1), mag_pair_cdf_rev1 0.3 dB (was 0.2). Values unchanged; only the stated bound. F1 rerun: 0 outside.", "docs/figures/redteam/F1_rev1.csv", "Low", "None"
"X4", "The C# tool applies each antenna's own fallback model", "Code reading for C3 (CouplingCalculations.cs)", "Antenna 2 uses otherAntModel2", "Finding", "CalculateWorstCaseCoupling_dBi passes otherAntModel1 for antenna 2, in the live and the commented-out versions. It matters only when the two antennas use different fallback models below the CSC^2 start angle.", "CouplingCalculations.cs (project)", "Low", "Report to the tool owner."
"L2", "Layer 2: independent review", "A separate agent that had not seen the plan, the register or the code history tried to break the results (2026-11-10)", "Findings fixed or reported", "Done", "Four of seven claims broken (R1 to R6 below); exact slices, bias at nodes, exact vs Monte Carlo, FFT vs direct and random mode vs time series held. R1, R3 and R4 fixed; R2 is C1.", "-", "Info", "Human review by NOAA and FAA subject-matter experts (layer 3)."
"R1", "Bidirectional sector scans are modeled correctly in the time series", "Review: sample the scan clock and the engine time series for a 30 deg bidirectional sector with 0.5 s turnaround", "Time series visits the sector as random mode assumes", "Fixed", "Each cut restarted the back-and-forth sweep at the same point, so part of the sector was never visited: per-volume maxima 13 to 39 dB low; engine time series 10 dB low at 99% with the path 10 deg off center. Now the sweep phase runs on across cuts and volumes. CDFs (random, exact) were not affected. Rerun of the review's check (30 deg sector, 10 deg off center, 4 trials of 300 s): per-trial maximum of the volume maxima 68.3, 45.5, 71.4, 45.5 dBi (was 29.3, 45.5, 32.6, 40.6; restart 66.3, 61.7, 71.5, 45.9). No committed output used a bidirectional time series (sector_bearing_sweep_rev1.csv holds CDF percentiles).", "test_sectors_rev1.m (testBidirectionalClockSweepsWholeSector, testEngineBidirectionalOddPassesMatchExact)", "High", "None"
"R2", "Inputs are conservative enough", "Review: ROC specification vs NSSL measured WSR-88D pattern", "Protection default is conservative", "Fixed", "Same as C1: the specification gives 10 to 15 dB more at 50 to 99%. Decision 7 (2026-11-10): nexrad_pattern now defaults to roc_spec; Rev 5 reproductions and pinned references set nssl_measured explicitly.", "docs/figures/redteam/C1_rev1.csv", "High", "None"
"R3", "The worst-mode envelope bounds any mix of scan modes", "Review: AVSET terminations above 6.4 deg against the 74-mode envelope, takeoff angles 6 to 10 deg", "Envelope at or above every operational mode", "Fixed", "A termination above 6.4 deg gives the cuts just above 6.4 deg more time than any of the 74 modes: envelope up to 4.35 dB low at 99.9% (2,000 random pairs; every pair with a WSR-88D takeoff angle above 6 deg changed, none below 5.5 deg by more than 0.55 dB). The default mode set now holds every termination: 216 modes; pair time about 2.0 ms against 1.75 ms. The 74-mode set stays available as ""lowest"" for the Rev 5 reproductions.", "test_mag_lookup_rev1.m (testAvsetTerminationsCloseEnvelopeGap)", "Medium", "Use the default (216-mode) set for Rev 6."
"R4", "Nearest-node and rounded angles are conservative", "Review: interp ""nearest"" and el_round_deg 0.01 and 0.1 against exact slices", "Never below the stated bias", "Fixed", "Nearest was up to 0.50 dB low, rounding to 0.1 deg up to 0.85 dB low. ""nearest"" is replaced by ""upper"" (larger of the two neighbors at each probability level); el_round_deg now gives each cell's upper envelope. Both are never below the unrounded quantile value.", "test_mag_lookup_rev1.m (testUpperAndRoundingNeverBelowQuantile)", "Medium", "None"
"R5", "Separable sector slices report the exact peak gain", "Review: brute-force peak of separable sector slices", "Doc matches behavior", "Fixed", "h.gmax was the largest azimuth sample, up to 0.007 dB below the true peak. Separable and grid slices now also evaluate az_off = 0 (when inside) and the interval ends: exact when the peak is the main beam or an end; the help says so.", "test_gain_slice_rev1.m (testSeparableSectorPeakIsExact)", "Low", "None"
"R6", "Input validation", "Review: duplicate percentiles, empty pair list, angle 1e-13 past the grid", "Clear errors", "Fixed", "Were cryptic errors (no silent wrong answer). Now mag_pair:duplicatePct and mag_pair:noPairs, and angles within the 1e-12 tolerance past the grid edge are clamped and run.", "test_mag_lookup_rev1.m (testInputErrors)", "Low", "None"
"G1", "Scope: aggregate interference", "None possible", "Stated in the report", "Scope limit", "One ATC radar per WSR-88D; no aggregate from several ATC radars.", "-", "Info", "State in Rev 6."
"G2", "Scope: polarization and dual-pol products", "None possible", "Stated in the report", "Scope limit", "No polarization mismatch and no dual-polarization product effects.", "-", "Info", "State in Rev 6."
"G3", "Scope: pulse-level coupling and FDR", "None possible", "Stated in the report", "Scope limit", "Antenna coupling only; no pulse-level coupling or frequency-dependent rejection.", "-", "Info", "State in Rev 6."
"G4", "Scope: validation against observed interference", "None possible with public data", "Stated in the report", "Scope limit", "No validation against observed interference; ROC interference records are restricted and only named here.", "-", "Info", "State in Rev 6."
};
T = cell2table(c, 'VariableNames', ["id", "claim", "test", "criterion", "result", "key_result", "evidence", "severity", "next_step"]);
T = convertvars(T, T.Properties.VariableNames, "string");
T.quick = zeros(height(T), 1);
end

% =========================================================================
function opts = fill_opts(opts)
here = string(fileparts(mfilename("fullpath")));
d = struct("repo_dir", here, "out_dir", "", "quick", false, "inputs", "rev5");
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || (isempty(opts.(f{i})) && ~islogical(opts.(f{i})))
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("redteam:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
opts.repo_dir = string(opts.repo_dir);
opts.inputs = string(opts.inputs);
if ~ismember(opts.inputs, ["rev5", "rev6"])
    error("redteam:badInputs", "opts.inputs must be ""rev5"" or ""rev6"".");
end
if opts.out_dir == "" && opts.inputs == "rev6"
    opts.out_dir = fullfile(opts.repo_dir, "docs", "figures", "redteam_rev6", "attacks");
elseif opts.out_dir == ""
    opts.out_dir = fullfile(opts.repo_dir, "docs", "figures", "redteam");
end
opts.out_dir = string(opts.out_dir);
end
