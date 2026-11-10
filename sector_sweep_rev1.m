function u = sector_sweep_rev1(tau_s, rate_deg_s, sector_deg, u0_deg, mode, turnaround_s)
%SECTOR_SWEEP_REV1 Beam position within a sector scan against time.
%   u = sector_sweep_rev1(tau_s, rate_deg_s, sector_deg, u0_deg, mode,
%   turnaround_s) returns the beam's distance u (deg, 0 to sector_deg) from
%   the sector's leading edge, tau_s seconds after the start of an
%   elevation step, for a beam that sweeps at |rate_deg_s|. It is the one
%   sweep law that radar_mutual_gain_cdf_rev1 and nexrad_asr11_mag_rev1 use
%   for time series.
%
%   mode
%     "restart"        (default) sweep the sector at constant rate and jump
%                      back to the leading edge: u = mod(u0 + rate*tau,
%                      sector). With sector_deg = 360 this is ordinary
%                      rotation.
%     "bidirectional"  sweep forward, dwell turnaround_s at the far edge,
%                      sweep back, dwell turnaround_s at the near edge, and
%                      repeat. One pass (one way) takes sector/rate +
%                      turnaround_s. The beam radiates during the dwell
%                      (bounding case: it is pointed at the edge).
%   No sweep law of an operational sector-scanning radar in scope is known
%   from a public source (TODO); both laws are bounding models.
%
%   tau_s         time since the start of the step, s (vector)
%   rate_deg_s    sweep rate, deg/s (scalar, or one per tau_s; sign ignored)
%   sector_deg    sector width, deg (scalar, 0 < sector_deg <= 360)
%   u0_deg        position at tau_s = 0, deg from the leading edge (scalar)
%   turnaround_s  dwell at each edge, s (bidirectional only; default 0)
%
%   Time share of the edges in random mode, per step: each edge holds
%   turnaround_s/(2*(sector/rate + turnaround_s)) of the time; the rest is
%   uniform across the sector (radar_mutual_gain_cdf_rev1,
%   mag_exact_cdf_rev1, gain_slice_rev1 edge_frac).
%
%   Example:
%     t = (0:0.1:20)';
%     u = sector_sweep_rev1(t, 30, 90, 0, "bidirectional", 1);
%     plot(t, u)
%
%   See also radar_mutual_gain_cdf_rev1, nexrad_asr11_mag_rev1, gain_slice_rev1.

arguments
    tau_s double {mustBeReal}
    rate_deg_s double {mustBeReal, mustBeNonzero}
    sector_deg (1, 1) double {mustBePositive}
    u0_deg (1, 1) double {mustBeReal} = 0
    mode (1, 1) string {mustBeMember(mode, ["restart", "bidirectional"])} = "restart"
    turnaround_s (1, 1) double {mustBeNonnegative} = 0
end
if sector_deg > 360
    error("sector_sweep:badSector", "sector_deg must be in (0, 360].");
end
r = abs(rate_deg_s);
if mode == "restart"
    u = mod(u0_deg + r.*tau_s, sector_deg);
    return
end
if sector_deg >= 360
    error("sector_sweep:bidirectional360", "A bidirectional sweep needs sector_deg below 360.");
end
t_s = sector_deg./r;                             % one-way sweep time
P = t_s + turnaround_s;                          % one pass including its end dwell
t = mod(tau_s + min(max(u0_deg, 0), sector_deg)./r, 2*P);
u = zeros(size(t));                              % near-edge dwell stays at 0
fwd = t < t_s;
u(fwd) = at(r, fwd).*t(fwd);
u(t >= t_s & t < P) = sector_deg;                % far-edge dwell
bwd = t >= P & t < P + t_s;
u(bwd) = sector_deg - at(r, bwd).*(t(bwd) - at(P, bwd));
u = min(max(u, 0), sector_deg);
end

% -------------------------------------------------------------------------
function x = at(v, m)
% v for the masked samples (v scalar or one value per sample).
if isscalar(v)
    x = v;
else
    x = v(m);
end
end
