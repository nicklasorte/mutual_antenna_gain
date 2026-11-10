%MAKE_MAG_LOOKUP_REFERENCE_REV1 Pin the two-part model's outputs (run once).
%   Writes test_data/mag_lookup_reference_rev1.csv: per-mode-pair and
%   envelope percentiles of mag_pair_cdf_rev1 for the ASR-11 -> WSR-88D
%   pair (all 74 VCP modes, Rev 5 study patterns) at five takeoff-angle
%   pairs, from tables built by gain_dist_table_rev1 with the azimuth-grid
%   method of commit 8c980b7. test_regression_rev5 checks later code
%   against this file. Rerun only to re-pin on purpose, and say why in the
%   commit message.
%
%   See also test_regression_rev5, gain_dist_table_rev1, mag_pair_cdf_rev1.

out_csv = fullfile(fileparts(which("make_mag_lookup_reference_rev1")), "test_data", ...
    "mag_lookup_reference_rev1.csv");
if ~isfolder(fileparts(out_csv)), mkdir(fileparts(out_csv)); end
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1));
[scans, names] = nexrad_mode_scans_rev1();
grid = -1:0.05:3;
Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("name", "WSR-88D", ...
    "mode_names", names, "el_grid_deg", grid));
Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("name", "ASR-11", "el_grid_deg", grid));
pairs = [0 0; 2.0 1.7; 1.73 2.02; -0.6 2.9; 0.5 0.5];
pct = [50 90 99 99.9 99.99 99.999];
R = mag_pair_cdf_rev1(Ta, Tn, pairs(:, 1), pairs(:, 2), struct("pct", pct, "per_mode", true));
n_p = size(pairs, 1);
n_m = numel(R.mode_pairs);
rows = cell(n_p, 1);
for ip = 1:n_p
    v = reshape(R.mode_pct(ip, :, :), n_m, numel(pct));
    env = R.pct{ip, 3:2 + numel(pct)};
    lab = [R.mode_pairs; "envelope"];
    vals = [v; env];
    [mi, pj] = ndgrid(1:n_m + 1, 1:numel(pct));
    rows{ip} = table(pairs(ip, 1) + zeros(numel(mi), 1), pairs(ip, 2) + zeros(numel(mi), 1), ...
        lab(mi(:)), pct(pj(:))', vals(:), 'VariableNames', ["el_tx_deg", "el_rx_deg", "mode_pair", "pct", "mag_dBi"]);
end
ref = vertcat(rows{:});
ref.max_on_path_dBi = repelem(R.pct.max_on_path_dBi, (n_m + 1)*numel(pct));
writetable(ref, out_csv);
fprintf("Wrote %d rows to %s\n", height(ref), out_csv);
