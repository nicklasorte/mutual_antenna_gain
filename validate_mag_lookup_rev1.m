%VALIDATE_MAG_LOOKUP_REV1 Exact two-step MAG against random-mode Monte Carlo.
%   Compares mag_pair_cdf_rev1 (exact convolution of two gain tables) with
%   radar_mutual_gain_cdf_rev1 random mode (n_mc samples) for the ASR-11 ->
%   WSR-88D pair in the Rev 5 worst mode (VCP 212 SAILSx3 AVSET 6.4), at
%   three takeoff-angle pairs: 0/0 (the Rev 5 baseline), 2.0/1.7 (on grid
%   nodes) and 1.73/2.02 (between nodes). Writes
%   docs/figures/mag_lookup_validation_rev1.csv with the exact value, the
%   Monte Carlo percentile and its order-statistic band (+/-3.5 sigma).
%
%   Settings: n_mc (default 2e7), out_csv.
%
%   See also test_mag_lookup_rev1, gain_dist_table_rev1, mag_pair_cdf_rev1.

if ~exist("n_mc", "var"), n_mc = 2e7; end
if ~exist("out_csv", "var")
    out_csv = fullfile(fileparts(which("validate_mag_lookup_rev1")), "docs", "figures", ...
        "mag_lookup_validation_rev1.csv");
end
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1));
scan_n = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
sys_a = o.sys_asr11;
sys_n = struct("name", "WSR-88D", "tx", o.sys_nexrad.tx, "rx", o.sys_nexrad.tx, "scan", scan_n);
grid = -5:0.05:10;
tic
Ta = gain_dist_table_rev1(sys_a.tx, sys_a.scan, struct("name", "ASR-11", "el_grid_deg", grid));
Tn = gain_dist_table_rev1(sys_n.tx, scan_n, struct("name", "WSR-88D", ...
    "mode_names", "VCP 212 SAILSx3 AVSET 6.4", "el_grid_deg", grid));
t_tab = toc;

pct = [50 90 99 99.9 99.99 99.999];
cases = [0 0; 2.0 1.7; 1.73 2.02];
rows = cell(size(cases, 1), 1);
t_ex = zeros(size(cases, 1), 1);
t_mc = zeros(size(cases, 1), 1);
for k = 1:size(cases, 1)
    tic
    R = mag_pair_cdf_rev1(Ta, Tn, cases(k, 1), cases(k, 2), struct("pct", pct));
    t_ex(k) = toc;
    tic
    mc = radar_mutual_gain_cdf_rev1(sys_a, sys_n, struct("seed", 11, "mode", "random", ...
        "n_samples", n_mc, "el12_deg", cases(k, 1), "el21_deg", cases(k, 2)));
    t_mc(k) = toc;
    g = sort(mc.g12_dBi);
    n = numel(g);
    p = pct(:)/100;
    h = 3.5*sqrt(n*p.*(1 - p));
    mc_q = prctile(g, pct(:));
    mc_lo = g(max(1, floor(n*p - h)));
    mc_hi = g(min(n, ceil(n*p + h)));
    exact = R.pct{1, 3:2 + numel(pct)}';
    rows{k} = table(repmat(cases(k, 1), numel(pct), 1), repmat(cases(k, 2), numel(pct), 1), pct(:), ...
        exact, mc_q, mc_lo(:), mc_hi(:), exact - mc_q, (exact >= mc_lo(:)) & (exact <= mc_hi(:) + R.bias_dB), ...
        'VariableNames', ["el_tx_deg", "el_rx_deg", "pct", "exact_dBi", "mc_dBi", "mc_lo_dBi", ...
        "mc_hi_dBi", "exact_minus_mc_dB", "in_band"]);
end
V = vertcat(rows{:});
writetable(V, out_csv);
disp(V)
fprintf("Tables %.1f s (once). Exact pair %.3f s avg; Monte Carlo (%g samples) %.1f s avg.\n", ...
    t_tab, mean(t_ex), n_mc, mean(t_mc));
fprintf("Grid bias bound %.2f dB; all in band: %d\n", R.bias_dB, all(V.in_band));
