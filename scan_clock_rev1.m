function C = scan_clock_rev1(spec, seed, trial, duration_s)
%SCAN_CLOCK_REV1 Random start phases and rates of one time-series trial.
%   C = scan_clock_rev1(spec, seed, trial, duration_s) draws the random
%   parts of the WSR-88D and ATC radar scan clocks for one trial: the
%   NEXRAD start azimuth and start time in its volume, the ATC start
%   azimuth, and the ATC rotation rate (one offset per trial, one per
%   revolution, or a correlated random walk). scan_clock_state_rev1(C, t)
%   then returns both radars' pointing at any times t. This is the one scan
%   clock behind the time series of nexrad_asr11_mag_rev1, the event
%   library (scan_event_library_rev1) and the pair time series
%   (mag_pair_timeseries_rev1), so all three see identical scans.
%
%   Trial j uses RandStream("mrg32k3a", Seed = seed) with substream
%   1000 + j, drawn in the order of nexrad_asr11_mag_rev1 rev 1 (Rev 5
%   outputs reproduce bit for bit).
%
%   spec struct (scan_clock_spec in nexrad_asr11_mag_rev1 builds it):
%     .rate_n, .el_n      NEXRAD rate (deg/s) and tilt (deg) per cut
%     .sec_n, .center_n   NEXRAD sector width and center (deg; 360 = full)
%     .sweep_n, .turn_n   NEXRAD sweep law and turnaround (sector_sweep_rev1)
%     .rate_a, .tilt_a    ATC nominal rate (deg/s) and tilt (deg)
%     .jit_a              ATC rate tolerance, deg/s (uniform +/- jit_a)
%     .jitter_mode        "per_trial" | "per_rev" | "random_walk"
%     .drift_step_a       random_walk only: standard deviation of the rate
%                         change from one revolution to the next, deg/s
%                         (the walk is held inside +/- jit_a)
%     .sec_a, .center_a   ATC sector width and absolute center (deg)
%
%   C struct: spec fields plus t_vol, cut_start, az0_n, t0_n, az0_a, and
%   r (per_trial) or r_rev and tb (per_rev, random_walk).
%
%   See also scan_clock_state_rev1, scan_event_library_rev1,
%   mag_pair_timeseries_rev1, nexrad_asr11_mag_rev1.

arguments
    spec (1, 1) struct
    seed (1, 1) double
    trial (1, 1) double {mustBeInteger, mustBePositive}
    duration_s (1, 1) double {mustBePositive}
end
C = spec;
dur_n = spec.sec_n./spec.rate_n + spec.turn_n;  % one revolution or one sector pass per cut
C.t_vol = sum(dur_n);
C.cut_start = [0; cumsum(dur_n(1:end-1))];
rs = RandStream("mrg32k3a", "Seed", seed);
rs.Substream = 1000 + trial;                    % apart from the CDF substreams
C.az0_n = 360*rand(rs);
C.t0_n = C.t_vol*rand(rs);
C.az0_a = 360*rand(rs);
rate_a = spec.rate_a;
jit = spec.jit_a;
switch spec.jitter_mode
    case "per_trial"
        C.r = rate_a + jit*(2*rand(rs) - 1);
    case {"per_rev", "random_walk"}
        n_rev = ceil((C.az0_a + (rate_a + jit)*duration_s)/360) + 1;
        if spec.jitter_mode == "per_rev"
            C.r_rev = rate_a + jit*(2*rand(rs, n_rev, 1) - 1);
        else
            % Correlated drift: a random walk in rate, one step per
            % revolution, held inside the tolerance band. No public
            % source gives the drive's drift statistics (TODO).
            r = zeros(n_rev, 1);
            r(1) = rate_a + jit*(2*rand(rs) - 1);
            step = spec.drift_step_a*randn(rs, n_rev, 1);
            for k = 2:n_rev
                r(k) = min(max(r(k - 1) + step(k), rate_a - jit), rate_a + jit);
            end
            C.r_rev = r;
        end
        C.tb = [0; cumsum(360./C.r_rev)] - C.az0_a/C.r_rev(1);   % revolution start times
    otherwise
        error("scan_clock:badJitterMode", "spec.jitter_mode must be ""per_trial"", ""per_rev"" or ""random_walk"".");
end
end
