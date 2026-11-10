function lib = scan_event_library_rev1(spec, opts)
%SCAN_EVENT_LIBRARY_REV1 Event windows of a WSR-88D / ATC scan pair (run once).
%   lib = scan_event_library_rev1(spec, opts) builds, once per pair of scan
%   modes, the time-series event library: for each trial, the scan clock
%   (scan_clock_rev1) and the fixed-grid steps at which either beam is
%   within a few degrees of the path (scan_event_windows_rev1). The library
%   does not depend on the antenna patterns or on the takeoff angles, so
%   mag_pair_timeseries_rev1 reuses it for any number of radar pairs that
%   share the two scan modes, evaluating gains only inside the windows.
%
%   spec  clock spec (scan_clock_rev1); nexrad_asr11_mag_rev1 returns the
%         one it used as out.scan_max.clock_spec
%   opts  struct:
%     .seed         required RNG seed
%     .n_trials     trials (default 2)
%     .first_trial  first trial number (default 1)
%     .duration_s   seconds per trial (default 3600)
%     .dt_s         fixed time step (default 0.002, the wrapper's step)
%     .delta_n_deg, .delta_a_deg, .dt_coarse_s  window options
%                   (scan_event_windows_rev1; defaults 4, 4, 0.02)
%
%   lib struct: spec, seed, trials, duration_s, dt_s, runs (cell, one
%   n-by-2 int32 array of [first last] grid indices per trial),
%   delta_eff_n_deg, delta_eff_a_deg (smallest over trials), share (mean
%   share of the grid kept), window_opts.
%
%   Example:
%     o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, ...
%         "duration_s", 60, "n_trials", 1));
%     lib = scan_event_library_rev1(o.scan_max.clock_spec, ...
%         struct("seed", 1, "n_trials", 10, "duration_s", 3600));
%
%   See also mag_pair_timeseries_rev1, scan_event_windows_rev1, scan_clock_rev1.

arguments
    spec (1, 1) struct
    opts (1, 1) struct
end
if ~isfield(opts, "seed")
    error("scan_events:missingSeed", "opts.seed is required.");
end
d = struct("n_trials", 2, "first_trial", 1, "duration_s", 3600, "dt_s", 0.002, ...
    "delta_n_deg", 4, "delta_a_deg", 4, "dt_coarse_s", 0.02);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
wo = struct("delta_n_deg", opts.delta_n_deg, "delta_a_deg", opts.delta_a_deg, "dt_coarse_s", opts.dt_coarse_s);
trials = opts.first_trial + (0:opts.n_trials - 1);
runs = cell(opts.n_trials, 1);
[de_n, de_a, share] = deal(zeros(opts.n_trials, 1));
for jj = 1:opts.n_trials
    C = scan_clock_rev1(spec, opts.seed, trials(jj), opts.duration_s);
    [idx, info] = scan_event_windows_rev1(C, opts.duration_s, opts.dt_s, wo);
    b = find(diff(idx) > 1);
    runs{jj} = int32([idx([1; b + 1]), idx([b; numel(idx)])]);
    de_n(jj) = info.delta_eff_n_deg;
    de_a(jj) = info.delta_eff_a_deg;
    share(jj) = info.share;
end
lib = struct("spec", spec, "seed", opts.seed, "trials", trials, "duration_s", opts.duration_s, ...
    "dt_s", opts.dt_s, "runs", {runs}, "delta_eff_n_deg", min(de_n), "delta_eff_a_deg", min(de_a), ...
    "share", mean(share), "window_opts", wo);
end
