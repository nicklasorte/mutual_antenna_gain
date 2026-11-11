function tests = test_regression_rev5
%TEST_REGRESSION_REV5 Pinned reference outputs for both MAG models.
%   Guards the optimization work (docs/mag_optimization_redteam_plan_rev1.md,
%   package K3):
%     - Original model: rows of docs/figures/vcp_sweep_rev5.csv, rerun with
%       the sweep's base settings (example_nexrad_asr11_rev1), must match
%       to 1e-10 relative. Any change to sampling or the clock fails here.
%     - Two-part model: test_data/mag_lookup_reference_rev1.csv (pinned by
%       make_mag_lookup_reference_rev1 at commit 8c980b7). New gain-slice
%       methods may only move a percentile up, by at most one bin per
%       table (they bound the pattern from above in each cell).
%
%   Run: results = runtests("test_regression_rev5");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
here = fileparts(which("test_regression_rev5"));
tc.TestData.sweep = readtable(fullfile(here, "docs", "figures", "vcp_sweep_rev5.csv"), "TextType", "string", ...
    "VariableNamingRule", "preserve");                % "case" is a keyword
tc.TestData.ref = readtable(fullfile(here, "test_data", "mag_lookup_reference_rev1.csv"), "TextType", "string");
end

% =========================================================================
function testOriginalModelSweepRows(tc)
% Paper baseline, a short clear-air mode, and the Rev 5 worst mode.
c = vcp_mode_cases_rev1("lowest");      % the 74 modes of the Rev 5 study
labels = string(cellfun(@(x) x.label, c, "UniformOutput", false));
want = ["Paper: F.1245, 3 rpm, horizon", "VCP 35", "VCP 212 SAILSx3 AVSET 6.4"];
base = struct("seed", 20261105, "n_samples", 5e6, "duration_s", 3600, "n_trials", 4, "nexrad_pattern", "nssl_measured", "cdf_method", "monte_carlo");   % Rev 5 inputs and method
thr = [40 50 60];
S = tc.TestData.sweep;
for lab = want
    i = find(labels == lab, 1);
    o = base;
    f = setdiff(fieldnames(c{i}), "label");
    for k = 1:numel(f)
        o.(f{k}) = c{i}.(f{k});
    end
    r = nexrad_asr11_mag_rev1(o);
    q = r.cdf.pct21;
    R = r.scan_max.asr_rev(r.scan_max.asr_rev.complete, :);
    got = [q.gain_dBi(q.pct == 99), q.gain_dBi(q.pct == 99.9), q.gain_dBi(q.pct == 99.99), ...
        100*mean(R.max_g21_dBi >= thr(1)), 100*mean(R.max_g21_dBi >= thr(2)), ...
        100*mean(R.max_g21_dBi >= thr(3)), max(r.scan_max.nexrad_volume.max_g21_dBi)];
    row = S(S.("case") == lab, :);
    ref = [row.p99, row.p99_9, row.p99_99, row.rev_pct_ge_40, row.rev_pct_ge_50, ...
        row.rev_pct_ge_60, row.vol_max];
    verifyEqual(tc, got, ref, "RelTol", 1e-10, "Rev 5 sweep row changed: " + lab);
end
end

function testTwoPartModelReference(tc)
ref = tc.TestData.ref;
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "nexrad_pattern", "nssl_measured"));
[scans, names] = nexrad_mode_scans_rev1("lowest");   % the 74 modes of the Rev 5 study
grid = -1:0.05:3;
Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("name", "WSR-88D", ...
    "mode_names", names, "el_grid_deg", grid));
Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("name", "ASR-11", "el_grid_deg", grid));
pairs = unique(ref(:, ["el_tx_deg", "el_rx_deg"]), "stable");
pct = unique(ref.pct)';
R = mag_pair_cdf_rev1(Ta, Tn, pairs.el_tx_deg, pairs.el_rx_deg, struct("pct", pct, "per_mode", true));
lab = [R.mode_pairs; "envelope"];
bin = Tn.bin_dB;
for ip = 1:height(pairs)
    v = [reshape(R.mode_pct(ip, :, :), numel(R.mode_pairs), numel(pct)); R.pct{ip, 3:2 + numel(pct)}];
    sel = ref.el_tx_deg == pairs.el_tx_deg(ip) & ref.el_rx_deg == pairs.el_rx_deg(ip);
    r = ref(sel, :);
    [~, im] = ismember(r.mode_pair, lab);
    [~, jp] = ismember(r.pct, pct);
    got = v(sub2ind(size(v), im, jp));
    d = got - r.mag_dBi;
    msg = sprintf("pair %.2f/%.2f deg", pairs.el_tx_deg(ip), pairs.el_rx_deg(ip));
    verifyGreaterThanOrEqual(tc, min(d), -1e-9, "Two-part percentile moved down: " + msg);
    verifyLessThanOrEqual(tc, max(d), 2*bin + 1e-9, "Two-part percentile moved up more than one bin per table: " + msg);
    % The pinned peak is the azimuth-grid maximum; exact slices may only
    % raise it, and only by the grid's peak loss (well under 0.01 dB).
    verifyGreaterThanOrEqual(tc, R.pct.max_on_path_dBi(ip), r.max_on_path_dBi(1) - 1e-9, "max on path: " + msg);
    verifyLessThanOrEqual(tc, R.pct.max_on_path_dBi(ip), r.max_on_path_dBi(1) + 0.01, "max on path: " + msg);
end
end
