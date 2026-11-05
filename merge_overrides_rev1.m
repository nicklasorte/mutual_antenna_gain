function o = merge_overrides_rev1(o, ov)
%MERGE_OVERRIDES_REV1 Apply run-size overrides to an options struct.
%   o = merge_overrides_rev1(o, ov) copies the run-size fields of ov
%   (n_samples, duration_s, n_trials, dt_s) into the options struct o. The
%   example scripts and report_figures_rev5 call it when the workspace
%   holds a struct named example_overrides, so the tests can run each
%   script end to end with small runs. Other fields of ov are ignored here;
%   report_figures_rev5 reads out_dir, n_h and n_list from it directly.
%
%   Example:
%     o = merge_overrides_rev1(struct("seed", 1, "n_samples", 5e6), struct("n_samples", 1e4));
%
%   See also example_nexrad_asr11_rev1, report_figures_rev5.

arguments
    o (1,1) struct
    ov (1,1) struct
end
f = ["n_samples", "duration_s", "n_trials", "dt_s"];
for k = 1:numel(f)
    if isfield(ov, f(k))
        o.(f(k)) = ov.(f(k));
    end
end
end
