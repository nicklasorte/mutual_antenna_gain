%EXAMPLE_MAG_LOOKUP_REV1 Two-step MAG: build gain tables once, then any pair.
%   Step 1 builds the exact gain distribution of the WSR-88D antenna (all
%   74 scan modes) and of the ATC radar antenna over the elevation grid
%   -5 to 10 deg in 0.05 deg steps, and saves both tables. Step 2 gets the
%   mutual antenna gain (MAG) CDF of any pair from the two takeoff angles
%   alone, as the worst-mode envelope over every WSR-88D scan mode.
%
%   Settings (set any of these in the workspace before running to change
%   them):
%     out_dir        folder for the saved tables (default tempdir/mag_lookup;
%                    the tables are about 30 MB, so keep them out of git)
%     atc_radar      "asr11" (default) or "asr9"
%     nexrad_pattern "nssl_measured" (default), "roc_spec" or "f1245"
%     el_tx, el_rx   takeoff angles of the worked pair (default 2.0, 1.7 deg)
%     n_batch        random pairs for the timing run (default 1000)
%     parallel_flag  1 builds the tables and runs the batch in parfor
%
%   The antenna patterns are the Rev 5 study patterns, taken from
%   nexrad_asr11_mag_rev1 (a short run returns its systems).
%
%   See also gain_dist_table_rev1, mag_pair_cdf_rev1, nexrad_mode_scans_rev1.

if ~exist("out_dir", "var"), out_dir = fullfile(tempdir, "mag_lookup"); end
if ~exist("atc_radar", "var"), atc_radar = "asr11"; end
if ~exist("nexrad_pattern", "var"), nexrad_pattern = "nssl_measured"; end
if ~exist("el_tx", "var"), el_tx = 2.0; end
if ~exist("el_rx", "var"), el_rx = 1.7; end
if ~exist("n_batch", "var"), n_batch = 1000; end
if ~exist("parallel_flag", "var"), parallel_flag = 0; end
if ~isfolder(out_dir), mkdir(out_dir); end

%% Step 1: one table per radar system (run once, then load)
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, ...
    "n_trials", 1, "atc_radar", atc_radar, "nexrad_pattern", nexrad_pattern));
[scans, names] = nexrad_mode_scans_rev1();
el_grid = -5:0.05:10;

tic
T_nexrad = gain_dist_table_rev1(o.sys_nexrad.tx, scans, struct("name", "WSR-88D " + nexrad_pattern, ...
    "mode_names", names, "el_grid_deg", el_grid, "parallel_flag", parallel_flag));
t_nexrad = toc;
tic
T_atc = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("name", o.sys_asr11.name, ...
    "el_grid_deg", el_grid, "parallel_flag", parallel_flag));
t_atc = toc;
f_nexrad = fullfile(out_dir, "gain_table_wsr88d_" + nexrad_pattern + "_rev1.mat");
f_atc = fullfile(out_dir, "gain_table_" + atc_radar + "_rev1.mat");
save(f_nexrad, "T_nexrad");
save(f_atc, "T_atc");
fprintf("Step 1: WSR-88D table %.1f s (%d modes, %d tilts), %s table %.1f s\n", ...
    t_nexrad, numel(names), numel(T_nexrad.tilt_deg), o.sys_asr11.name, t_atc);

%% Step 2: one pair, ATC radar (tx) -> WSR-88D (rx), worst-mode envelope
R = mag_pair_cdf_rev1(T_atc, T_nexrad, el_tx, el_rx, struct("return_curve", true));
fprintf("\nStep 2: %s -> WSR-88D, takeoff angles %.2f / %.2f deg (envelope of %d modes)\n", ...
    o.sys_asr11.name, el_tx, el_rx, numel(R.mode_pairs));
disp(R.pct)
disp(table(R.pct_list(:), R.worst_mode(1, :)', 'VariableNames', ["pct", "worst_mode"]))

%% Step 2 at scale: many pairs with random takeoff angles
rs = RandStream("mrg32k3a", "Seed", 1);
et = -1 + 6*rand(rs, n_batch, 1);
er = -1 + 6*rand(rs, n_batch, 1);
tic
Rb = mag_pair_cdf_rev1(T_atc, T_nexrad, et, er, struct("parallel_flag", parallel_flag));
t_batch = toc;
fprintf("\n%d pairs in %.1f s (%.1f ms per pair)\n", n_batch, t_batch, 1e3*t_batch/n_batch);
