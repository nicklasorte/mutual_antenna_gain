function [ref_pair, ref_ts] = rev6_reference_rev1()
%REV6_REFERENCE_REV1 Outputs pinned by make_rev6_reference_rev1.
%   [ref_pair, ref_ts] = rev6_reference_rev1() recomputes the two pinned
%   tables (see make_rev6_reference_rev1) with the current code.
%
%   See also make_rev6_reference_rev1, test_regression_rev6.
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
    "nexrad_pattern", "roc_spec"));
[scans, names] = nexrad_mode_scans_rev1("all");
grid = -1:0.05:3;
Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("mode_names", names, "el_grid_deg", grid));
Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("el_grid_deg", grid));
pairs = [0 0; 0.5 0.5; 2.0 1.7; 1.73 2.02; -0.6 2.9; 0.27 0.83];
pct = [50 90 99 99.9 99.99 99.999];
R = mag_pair_cdf_rev1(Ta, Tn, pairs(:, 1), pairs(:, 2), struct("pct", pct));
n_p = size(pairs, 1);
[ip, pj] = ndgrid(1:n_p, 1:numel(pct));
env = R.pct{:, 3:2 + numel(pct)};
ref_pair = table(pairs(ip(:), 1), pairs(ip(:), 2), pct(pj(:))', env(:), R.worst_mode(sub2ind(size(R.worst_mode), ip(:), pj(:))), ...
    R.pct.max_on_path_dBi(ip(:)), 'VariableNames', ["el_tx_deg", "el_rx_deg", "pct", "mag_dBi", "worst_mode", "max_on_path_dBi"]);

c = vcp_mode_cases_rev1();
lab = cellfun(@(x) string(x.label), c);
k = find(lab == "VCP 212 SAILSx3 AVSET 6.4", 1);
ow = rmfield(c{k}, "label");
ow.seed = 1;
ow.n_samples = 1e3;
ow.duration_s = 60;
ow.n_trials = 1;
ow.cdf_method = "exact";
w = nexrad_asr11_mag_rev1(ow);
lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 20261105, "n_trials", 4, "duration_s", 3600));
S = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0.5, 0.3, ...
    struct("radial_thr_dBi", [40 50 60 70]));
mk = @(T, kind, scan) table(repmat(string(kind), height(T), 1), T.trial, scan, T.max_g12_dBi, T.max_g21_dBi, T.complete, ...
    'VariableNames', ["kind", "trial", "scan", "max_g12_dBi", "max_g21_dBi", "complete"]);
C = S.nexrad_cut;
n_cut = max(C.cut);
H = S.radial_cuts;
hv = H{:, startsWith(H.Properties.VariableNames, "n_radials_ge_")};
harm = table(repmat("harm", height(H), 1), H.trial, (H.volume - 1)*n_cut + H.cut, hv(:, 2), hv(:, 3), H.complete, ...
    'VariableNames', ["kind", "trial", "scan", "max_g12_dBi", "max_g21_dBi", "complete"]);   % harm rows: n at 50 and 60 dBi
ref_ts = [mk(S.asr_rev, "rev", S.asr_rev.rev); mk(C, "cut", (C.volume - 1)*n_cut + C.cut); ...
    mk(S.nexrad_volume, "volume", S.nexrad_volume.volume); harm];
end
