function T = mag_sensitivity_rev1(base, cases)
%MAG_SENSITIVITY_REV1 One-at-a-time sensitivity of the WSR-88D/ASR-11 MAG.
%   T = mag_sensitivity_rev1(base, cases) runs nexrad_asr11_mag_rev1 for
%   the options struct base, then once for each case, changing only the
%   fields that case sets. It answers the request at the 2.7 GHz Retuning
%   LOE meeting of November 4, 2026 for a sensitivity section: how far each
%   antenna or scan assumption moves the MAG, in dB, at each percentile and
%   against main beam to main beam. A MAG change applies one for one to the
%   required basic transmission loss.
%
%   Inputs:
%     base   options struct for nexrad_asr11_mag_rev1 (seed required). Every
%            case reuses its seed, so differences come from the assumption,
%            not from sampling.
%     cases  struct array with fields group (string), label (string), and
%            opts (struct of option overrides)
%
%   Output table T, one row per run (baseline first):
%     group, label
%     p99_dBi, p99_9_dBi, p99_99_dBi, p99_999_dBi   time-weighted MAG
%     d99_dB ... d99_999_dB    change from the baseline at each percentile
%     mbmb_peak_dBi            peak-gain main beam to main beam
%     mbmb_path_dBi            main beam to main beam on the path
%     rev_pct_ge_45_dBi        % of ASR-11 revolutions whose maximum
%                              reaches 45 dBi
%     max_run_45_revs          longest run of consecutive such revolutions
%     max_run_45_vols          longest run of consecutive NEXRAD volumes
%                              whose maximum reaches 45 dBi
%
%   Example:
%     b = struct("seed", 1, "vcp", 212, "sails", 3, "n_samples", 1e5, "duration_s", 600, "n_trials", 1);
%     c = struct("group", "ASR-11", "label", "tilt +2.5", "opts", struct("asr_tilt_deg", 2.5));
%     T = mag_sensitivity_rev1(b, c);
%
%   See also nexrad_asr11_mag_rev1, example_sensitivity_rev1.

arguments
    base (1,1) struct
    cases struct
end
for k = 1:numel(cases)
    if ~all(isfield(cases(k), ["group", "label", "opts"]))
        error("radar_coupling:badCase", "Each case needs fields group, label and opts.");
    end
end
if isfield(base, "persist_thr_dBi") && ~ismember(45, base.persist_thr_dBi)
    base.persist_thr_dBi = unique([45, base.persist_thr_dBi(:)']);   % the default grid has 45
end

all_cases = [struct("group", "Baseline", "label", "Baseline", "opts", struct()), cases(:)'];
n = numel(all_cases);
rows = cell(n, 1);
for k = 1:n
    o = base;
    f = fieldnames(all_cases(k).opts);
    for i = 1:numel(f)
        o.(f{i}) = all_cases(k).opts.(f{i});
    end
    r = nexrad_asr11_mag_rev1(o);
    q = r.cdf.pct21;
    p = arrayfun(@(x) q.gain_dBi(abs(q.pct - x) < 1e-9), [99 99.9 99.99 99.999]);
    P = r.scan_max.persistence(r.scan_max.persistence.threshold_dBi == 45, :);
    rows{k} = {string(all_cases(k).group), string(all_cases(k).label), p(1), p(2), p(3), p(4), ...
        r.sensitivity.mag_dBi(1), r.ref_levels.mag_dBi(1), P.rev_pct, P.max_run_revs, ...
        r.scan_max.persistence_volume.max_run_vols(r.scan_max.persistence_volume.threshold_dBi == 45)};
end
T = cell2table(vertcat(rows{:}), "VariableNames", ["group", "label", "p99_dBi", "p99_9_dBi", ...
    "p99_99_dBi", "p99_999_dBi", "mbmb_peak_dBi", "mbmb_path_dBi", "rev_pct_ge_45_dBi", "max_run_45_revs", ...
    "max_run_45_vols"]);
d = T{:, ["p99_dBi", "p99_9_dBi", "p99_99_dBi", "p99_999_dBi"]} - T{1, ["p99_dBi", "p99_9_dBi", "p99_99_dBi", "p99_999_dBi"]};
T = addvars(T, d(:, 1), d(:, 2), d(:, 3), d(:, 4), 'After', "p99_999_dBi", ...
    'NewVariableNames', ["d99_dB", "d99_9_dB", "d99_99_dB", "d99_999_dB"]);
end

