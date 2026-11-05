function D = alignment_dwell_rev1(rate1_deg_s, rate2_deg_s, bw1_deg, bw2_deg, opts)
%ALIGNMENT_DWELL_REV1 How long two rotating beams stay aligned main beam to main beam.
%   D = alignment_dwell_rev1(rate1_deg_s, rate2_deg_s, bw1_deg, bw2_deg)
%   answers the synchronization concern raised at the 2.7 GHz Retuning LOE
%   meeting of November 4, 2026: radars with equal scan rates, or rates that
%   are integer multiples (for example 5 and 10 rpm), can drift into
%   relative synchronization and couple main beam to main beam for a long
%   time (up to about a month, a meeting statement recorded as an
%   assumption, source LOE_20261104 in nexrad_asr11_mag_rev1).
%
%   Two antennas rotate at rates r1 and r2 (deg/s). Radar 1 points at radar
%   2 when its azimuth is 0; radar 2 points at radar 1 at 180. Main beam to
%   main beam (MB-MB) means both beams are within their 3 dB beamwidths of
%   the path at the same instant, so the MAG is within 6 dB of the
%   peak-gain MB-MB value. Take the nearest ratio r1/r2 ~ p/q (p, q
%   coprime, both at most opts.max_ratio). The joint pointing repeats every
%   joint period T_j = 360 p / r1, with one MB-MB chance per period. The
%   phase psi = q*a1 - p*a2 drifts at
%       delta = q*r1 - p*r2                         (deg/s)
%   and MB-MB occurs in each period while psi stays inside a window of
%       W = q*bw1 + p*bw2                           (deg)
%   So once aligned, MB-MB repeats for an episode of
%       T_ep = W / |delta|  seconds, about T_ep / T_j consecutive periods.
%   Equal nominal rates give delta = 0 and an unbounded episode; only the
%   rate difference between the two drives ends it.
%
%   Rate tolerance: each drive runs at its nominal rate plus a constant
%   offset, uniform in +/- opts.jitter_deg_s (M.2069-0 Annex 1 drive
%   accuracy +/-0.5 deg/s; the FAA ASR-11 tolerance is TODO, see
%   nexrad_asr11_mag_rev1). A constant offset is the persistent case: a
%   wandering rate shortens episodes. Then delta = delta0 + q*d1 - p*d2, and
%   P(T_ep > T) = P(|delta| < W/T), from the distribution of the sum of the
%   two uniform offsets (trapezoidal).
%
%   Inputs:
%     rate1_deg_s, rate2_deg_s  nominal rotation rates (deg/s, > 0)
%     bw1_deg, bw2_deg          3 dB azimuth beamwidths (deg)
%   Name-value options:
%     jitter_deg_s  [j1 j2] half-widths of the uniform rate offsets
%                   (default [0.5 0.5] deg/s)
%     t_query_s     episode lengths to report (default 1 min, 10 min, 1 h,
%                   1 day, 30 days)
%     max_ratio     largest p or q considered (default 10)
%     max_episode_s upper limit on any episode (default Inf). Use the cut
%                   length for a radar that changes rate every cut (a
%                   WSR-88D cut lasts at most about 29 s), since alignment
%                   cannot outlast the rate.
%
%   Output struct D:
%     .ratio            [p q]
%     .joint_period_s   T_j
%     .window_deg       W
%     .delta0_deg_s     nominal phase drift rate
%     .episode0_s       nominal episode (Inf for exactly commensurate rates)
%     .events0          nominal MB-MB periods per episode
%     .p_exceed         table: t_s, p_episode_longer (probability, over
%                       the drive offsets, that an aligned episode lasts
%                       longer than t_s), periods (t_s / T_j)
%
%   Example (two ASR-11 radars, 12.5 rpm, 1.45 deg):
%     D = alignment_dwell_rev1(75, 75, 1.45, 1.45); disp(D.p_exceed)
%
%   See also nexrad_asr11_mag_rev1, radar_mutual_gain_cdf_rev1.

arguments
    rate1_deg_s (1,1) double {mustBePositive}
    rate2_deg_s (1,1) double {mustBePositive}
    bw1_deg (1,1) double {mustBePositive}
    bw2_deg (1,1) double {mustBePositive}
    opts.jitter_deg_s (1,2) double {mustBeNonnegative} = [0.5 0.5]
    opts.t_query_s (:,1) double {mustBePositive} = [60; 600; 3600; 86400; 30*86400]
    opts.max_ratio (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.max_episode_s (1,1) double {mustBePositive} = Inf
end
if any(opts.jitter_deg_s >= [rate1_deg_s rate2_deg_s])
    error("radar_coupling:badJitter", "jitter_deg_s must be below each rotation rate.");
end

% ---- Nearest ratio p:q with the smallest phase drift per joint period
best = Inf;
for q = 1:opts.max_ratio
    p = max(1, round(q*rate1_deg_s/rate2_deg_s));
    if p > opts.max_ratio || gcd(p, q) ~= 1
        continue
    end
    d = abs(q*rate1_deg_s - p*rate2_deg_s);
    if d < best - 1e-12
        best = d;
        pq = [p q];
    end
end
if ~isfinite(best)
    error("radar_coupling:noRatio", "No ratio p:q with p, q <= %d near %.4g.", opts.max_ratio, rate1_deg_s/rate2_deg_s);
end
p = pq(1);
q = pq(2);
D.ratio = pq;
D.joint_period_s = 360*p/rate1_deg_s;
D.window_deg = q*bw1_deg + p*bw2_deg;
D.delta0_deg_s = q*rate1_deg_s - p*rate2_deg_s;
D.episode0_s = min(D.window_deg/abs(D.delta0_deg_s), opts.max_episode_s);
D.events0 = D.episode0_s/D.joint_period_s;

% ---- Probability an episode outlasts t, over the drive offsets
A = q*opts.jitter_deg_s(1);
B = p*opts.jitter_deg_s(2);
a = D.window_deg./opts.t_query_s;               % |delta| must stay below this
P = sum_uniform_cdf(a - D.delta0_deg_s, A, B) - sum_uniform_cdf(-a - D.delta0_deg_s, A, B);
P(opts.t_query_s >= opts.max_episode_s) = 0;
D.p_exceed = table(opts.t_query_s, P, opts.t_query_s/D.joint_period_s, ...
    'VariableNames', ["t_s", "p_episode_longer", "periods"]);
end

% -------------------------------------------------------------------------
function F = sum_uniform_cdf(s, A, B)
% CDF of X + Y, X ~ U(-A, A), Y ~ U(-B, B) independent (trapezoidal
% density). A = B = 0 is a point mass at 0.
F = zeros(size(s));
if A == 0 && B == 0
    F = double(s > 0);
    return
end
if A < B
    [A, B] = deal(B, A);
end
if B == 0
    F = min(max((s + A)/(2*A), 0), 1);
    return
end
for i = 1:numel(s)
    % F(s) = E_Y[ F_X(s - Y) ], F_X uniform CDF, integrated in closed form
    % via the antiderivative of the clamped ramp.
    u = s(i) - [B, -B];                         % s - y at y = B and y = -B
    F(i) = (ramp_int(u(2), A) - ramp_int(u(1), A))/(2*B);
end
F = min(max(F, 0), 1);
end

function R = ramp_int(u, A)
% Integral from -Inf to u of the uniform U(-A, A) CDF.
if u <= -A
    R = 0;
elseif u >= A
    R = A + (u - A);                            % area of the ramp (A) plus the flat part
else
    R = (u + A)^2/(4*A);
end
end
