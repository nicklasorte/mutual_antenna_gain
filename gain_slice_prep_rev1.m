function S = gain_slice_prep_rev1(ant, opts)
%GAIN_SLICE_PREP_REV1 Prepare an antenna for exact gain-distribution slices.
%   S = gain_slice_prep_rev1(ant, opts) caches what gain_slice_rev1 needs to
%   return the time-weighted gain distribution of a mechanically rotating
%   antenna at one beam tilt and one elevation of the other station. It
%   picks the fastest exact method the pattern allows and checks it:
%
%   "closed_form"  Patterns that depend only on the off-axis angle gamma
%                  (statgain, f1245, parabolic, offaxis_table). For a beam
%                  sweeping azimuth uniformly, the share of time with
%                  off-axis angle <= gamma is known exactly:
%                    phi(gamma) = acos((cos(gamma) - sin(el)*sin(tilt)) /
%                                      (cos(el)*cos(tilt)))
%                    share     = phi(gamma)/180 deg (full rotation)
%                  (spherical law of cosines, as off_axis_deg in
%                  antenna_gain_rev1). The pattern is cut at breakpoints in
%                  gamma; each cell gets its share of time and the larger of
%                  its two end gains. For a table (linear between nodes)
%                  that is the exact maximum over the cell, so the result
%                  bounds the true distribution from above. No azimuth grid.
%   "separable"    Aperture patterns built from an azimuth cut and an
%                  elevation cut (rect, measured): g = g_max + max(g_az(az)
%                  + g_el(el, tilt), floor), M.1851-2 Eq. (15). The azimuth
%                  cut is sampled once; every slice is that sample shifted
%                  by g_el and clamped at the floor. One sort, then each
%                  slice costs a search over the gain bins only.
%   "grid"         Anything else (planar arrays): evaluate antenna_gain_rev1
%                  on the azimuth grid, as gain_dist_table_rev1 rev 1 did.
%
%   The separable method is used only if the pattern passes a numerical
%   check of the separable form at random angles; otherwise it falls back
%   to "grid" and says so in S.notes.
%
%   ant   antenna model struct (antenna_prepare_rev1 input)
%   opts  struct, all fields optional:
%     .bin_dB       gain bin width, dB (default 0.05)
%     .az_step_deg  azimuth step for "separable" and "grid" (default 0.001)
%     .gam_step_deg off-axis breakpoint step for analytic closed-form
%                   patterns (statgain, f1245, parabolic; default 0.001)
%     .method       "auto" (default), or force "closed_form", "separable"
%                   or "grid" (for tests)
%
%   S struct: .pa (prepared antenna), .method, .bin_dB, and the method's
%   cache (.gam_b and .g_b; or .az, .g_az, .g_az_sorted; or .az). Pass S to
%   gain_slice_rev1.
%
%   See also gain_slice_rev1, gain_dist_table_rev1, antenna_gain_rev1.

arguments
    ant (1, 1) struct
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
pa = antenna_prepare_rev1(ant);
S = struct("pa", pa, "method", "", "bin_dB", opts.bin_dB, "az_step_deg", opts.az_step_deg, ...
    "notes", strings(0, 1));
n_az = round(360/opts.az_step_deg);
az = ((1:n_az)' - 0.5)*360/n_az - 180;          % midpoint grid of az_off, (-180, 180]
S.az = az;

method = opts.method;
if method == "auto"
    switch pa.type
        case {"statgain", "f1245", "parabolic", "offaxis_table"}
            method = "closed_form";
        case {"rect", "measured"}
            method = "separable";
        otherwise
            method = "grid";
    end
end

switch method
    case "closed_form"
        if ~ismember(pa.type, ["statgain", "f1245", "parabolic", "offaxis_table"])
            error("gain_slice:badMethod", "closed_form needs a pattern of the off-axis angle only (got %s).", pa.type);
        end
        if pa.type == "offaxis_table"
            % griddedInterpolant "linear": the table nodes are the exact
            % breakpoints, and each cell is monotone between them.
            gam_b = pa.off_interp.GridVectors{1}(:);
        else
            % Analytic pattern: fine breakpoints. Between them a smooth
            % pattern is monotone except at a lobe peak, where the error is
            % second order in the step (below 1e-4 dB at 0.001 deg).
            gam_b = (0:opts.gam_step_deg:180)';
            if gam_b(end) < 180, gam_b(end + 1) = 180; end
        end
        z = zeros(size(gam_b));
        % Off-axis angle = az_off when el = tilt = 0 (off_axis_deg).
        S.gam_b = gam_b;
        S.g_b = antenna_gain_rev1(pa, 0, 0, struct("mech_az", -gam_b, "tilt", z, ...
            "steer_az", z, "steer_el", z));
    case "separable"
        z = zeros(n_az, 1);
        g0 = antenna_gain_rev1(pa, 0, 0, struct("mech_az", 0, "tilt", 0, "steer_az", 0, "steer_el", 0));
        S.g_max = pa.g_max_dBi;
        S.floor_dB = pa.floor_dB;
        S.g_az = antenna_gain_rev1(pa, 0, 0, struct("mech_az", -az, "tilt", z, "steer_az", z, "steer_el", z)) - S.g_max;
        ok = abs(g0 - S.g_max) < 1e-9 && separable_check(pa, S);
        if ok
            S.g_az_sorted = sort(S.g_az);
            % Unique values and how many samples are <= each, for the
            % binary-search bin counts in gain_slice_rev1.
            [S.u_az, ia] = unique(S.g_az_sorted, "last");
            S.cnt_az = ia;
        else
            method = "grid";
            S.notes(end + 1) = "Separable check failed; using the azimuth grid.";
            S = rmfield(S, ["g_az", "g_max", "floor_dB"]);
        end
    case "grid"
    otherwise
        error("gain_slice:badMethod", "Unknown method ""%s"".", method);
end
S.method = string(method);
end

% =========================================================================
function ok = separable_check(pa, S)
% g(az, el, tilt) must equal g_max + max(g_az(az) + g_el(el, tilt), floor),
% with g_az from el = tilt = 0 and g_el from az = 0, both <= 0 dB.
rs = RandStream("mrg32k3a", "Seed", 5);
n = 400;
az = 360*rand(rs, n, 1) - 180;
el = 20*rand(rs, n, 1) - 5;
tilt = 6*rand(rs, n, 1) - 1;
z = zeros(n, 1);
g = antenna_gain_rev1(pa, 0, el, struct("mech_az", -az, "tilt", tilt, "steer_az", z, "steer_el", z));
g_az = antenna_gain_rev1(pa, 0, 0, struct("mech_az", -az, "tilt", z, "steer_az", z, "steer_el", z)) - S.g_max;
g_el = zeros(n, 1);
for i = 1:n
    g_el(i) = antenna_gain_rev1(pa, 0, el(i), struct("mech_az", 0, "tilt", tilt(i), "steer_az", 0, "steer_el", 0)) - S.g_max;
end
pred = S.g_max + max(g_az + g_el, S.floor_dB);
ok = all(abs(pred - g) < 1e-9) && max(S.g_az) <= 1e-9 && max(g_el) <= 1e-9;
end

% -------------------------------------------------------------------------
function opts = fill_opts(opts)
d = struct("bin_dB", 0.05, "az_step_deg", 0.001, "gam_step_deg", 0.001, "method", "auto");
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i})
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("gain_slice:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
opts.method = string(opts.method);
if ~isscalar(opts.bin_dB) || opts.bin_dB <= 0 || opts.bin_dB > 1
    error("gain_slice:badBin", "opts.bin_dB must be a scalar in (0, 1] dB.");
end
if ~isscalar(opts.az_step_deg) || opts.az_step_deg <= 0 || opts.az_step_deg > 0.1
    error("gain_slice:badAzStep", "opts.az_step_deg must be a scalar in (0, 0.1] deg.");
end
if ~isscalar(opts.gam_step_deg) || opts.gam_step_deg <= 0 || opts.gam_step_deg > 0.1
    error("gain_slice:badGamStep", "opts.gam_step_deg must be a scalar in (0, 0.1] deg.");
end
end
