function [st_n, st_a, ck] = scan_clock_state_rev1(C, t)
%SCAN_CLOCK_STATE_REV1 WSR-88D and ATC radar pointing at given times.
%   [st_n, st_a, ck] = scan_clock_state_rev1(C, t) evaluates the scan clock
%   C (scan_clock_rev1) at times t (s, column vector, any subset of the
%   trial). Every output is elementwise in t, so a subset of a time grid
%   gives exactly the states of the full grid at those times.
%
%   st_n, st_a  pointing states for antenna_gain_rev1: mech_az, tilt,
%               steer_az, steer_el (N-by-1 each)
%   ck struct:
%     .vol    NEXRAD volume number (1 = the volume that holds t = 0)
%     .cut    NEXRAD cut within its volume
%     .rev    ATC revolution (or sector sweep) number
%     .r_rev  ATC rate of each revolution, deg/s (indexed by rev)
%
%   Geometry: the ATC radar is at azimuth 0 from the NEXRAD and the NEXRAD
%   at azimuth 180 from the ATC radar (radar_mutual_gain_cdf_rev1).
%
%   See also scan_clock_rev1, sector_sweep_rev1.

arguments
    C (1, 1) struct
    t (:, 1) double
end
n_t = numel(t);
% ---- NEXRAD: volume, cut, azimuth (one revolution or sector pass per cut)
tau = t + C.t0_n;
vol = floor(tau/C.t_vol) + 1;
w = tau - (vol - 1)*C.t_vol;
cut = discretize(w, [C.cut_start; C.t_vol]);
cut(isnan(cut)) = numel(C.rate_n);
sec = C.sec_n;
if sec < 360
    u = sector_sweep_rev1(w - C.cut_start(cut), C.rate_n(cut), sec, C.az0_n/360*sec, C.sweep_n, C.turn_n);
    st_n = state(mod(C.center_n + u - sec/2, 360), C.el_n(cut));
else
    st_n = state(mod(C.az0_n + C.rate_n(cut).*(w - C.cut_start(cut)), 360), C.el_n(cut));
end
% ---- ATC radar: revolution counted from azimuth 0; one pass of the NEXRAD
%      bearing per revolution
switch C.jitter_mode
    case "per_trial"
        r = C.r;
        if C.sec_a < 360
            % Sector: distance swept since the first sweep's leading edge;
            % each sweep is one revolution.
            sw = C.az0_a/360*C.sec_a + r*t;
            rev = floor(sw/C.sec_a) + 1;
            u_a = sector_sweep_rev1(sw - (rev - 1)*C.sec_a, 1, C.sec_a);
            ang = mod(C.center_a + u_a - C.sec_a/2, 360);
        else
            ang = C.az0_a + r*t;
            rev = floor(ang/360) + 1;
        end
        r_rev = r + zeros(max(rev), 1);
    otherwise                                   % per_rev, random_walk
        r_rev = C.r_rev;
        rev = discretize(t, C.tb);
        ang = 360*(rev - 1) + (t - C.tb(rev)).*r_rev(rev);
end
st_a = state(mod(ang, 360), C.tilt_a + zeros(n_t, 1));
ck = struct("vol", vol, "cut", cut, "rev", rev, "r_rev", r_rev);
end

% -------------------------------------------------------------------------
function st = state(az, tilt)
% As nexrad_asr11_mag_rev1 rev 1.
z = zeros(numel(az), 1);
st = struct("mech_az", az(:), "tilt", tilt(:), "steer_az", z, "steer_el", z);
end
