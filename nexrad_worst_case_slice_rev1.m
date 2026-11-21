function [slice, info] = nexrad_worst_case_slice_rev1(opts)
%NEXRAD_WORST_CASE_SLICE_REV1 WSR-88D (NEXRAD) worst-case horizontal gain slice.
%   slice = nexrad_worst_case_slice_rev1() returns the WSR-88D antenna gain
%   toward a station on the horizon (0 deg elevation) as the antenna turns
%   through a full revolution, as an N-by-2 array:
%
%     slice(:, 1)  azimuth of the beam relative to the bearing of the
%                  station, deg, 0:step_deg:360 (0 and 360 = beam on the
%                  station)
%     slice(:, 2)  antenna gain toward the station, dBi
%
%   Worst case at the horizon. The default beam elevation (tilt_deg = 0)
%   puts the beam peak on the horizon, so the slice passes through the
%   45.5 dBi peak. Any other tilt gives a lower gain at every azimuth, so
%   tilt_deg = 0 bounds the gain toward any station on the horizon, as a
%   coordination zone that protects the radar needs. The operational
%   lowest cut is 0.5 deg (tilt_deg = 0.5 gives 42.0 dBi on the horizon).
%
%   Pattern (pencil beam, rotationally symmetric about the beam axis):
%
%     G(theta) = G_max - 12 (theta / bw)^2  until it meets the envelope,
%     G(theta) = G_max + E(theta)           beyond,
%
%   theta is the off-axis angle, cos(theta) = cos(tilt) cos(az). The
%   parabolic main beam is the M.1851-2 and M.1849-3 main-lobe form. E is
%   the sidelobe envelope below the peak:
%
%     "roc_spec" (default): -29 dB to 10 deg, -40 dB beyond (the ROC
%        specification; an upper bound, so the worst case)
%     "nssl_measured": -35 dB to 2 deg, linear to -48 dB at 6 deg and
%        -55 dB at 20 deg, -55 dB to 180 deg (KOUN measurements)
%
%   The envelope replaces the main beam from the angle where they meet, so
%   the first null is filled, which is conservative.
%
%   Name-value options:
%     step_deg  azimuth step (0.1); 360 must be a whole number of steps
%     pattern   "roc_spec" (default) | "nssl_measured"
%     tilt_deg  beam elevation above the horizon (0, the worst case)
%
%   info fields: name, g_max_dBi, bw_deg, pattern, tilt_deg, env_deg,
%   env_dB, x_deg (off-axis angle where the main beam meets the
%   envelope), peak_dBi, floor_dBi, far_field_m, params (table of every
%   value with its source), sources (table of citations with URLs).
%
%   Sources: ROC_TECH (45.5 dBi, 0.925 deg, sidelobe specification, 8.53 m
%   reflector at 2850 MHz), NSSL_1998 (measured sidelobes), M1849
%   (45 dBi, 0.92 deg; main-lobe form), M1851 (main-lobe form).
%
%   Example:
%     s = nexrad_worst_case_slice_rev1();
%     plot(s(:,1), s(:,2)); xlabel("Azimuth from station (deg)"); ylabel("Gain (dBi)")
%
%   DRAFT//PRE-DECISIONAL
%
%   See also asr11_worst_case_slice_rev1, asr9_worst_case_slice_rev1.

arguments
    opts.step_deg (1,1) double {mustBePositive} = 0.1
    opts.pattern (1,1) string {mustBeMember(opts.pattern, ["roc_spec", "nssl_measured"])} = "roc_spec"
    opts.tilt_deg (1,1) double {mustBeInRange(opts.tilt_deg, -1, 20)} = 0
end
n_step = 360/opts.step_deg;
if abs(n_step - round(n_step)) > 1e-9
    error("nexrad_worst_case_slice:badStep", "360 must be a whole number of steps; got step_deg = %g.", opts.step_deg);
end

g_max = 45.5;       % ROC_TECH: 45.5 dB at 2850 MHz including radome (M1849: 45.0)
bw = 0.925;         % ROC_TECH: 3 dB beamwidth 0.925 deg at 2850 MHz (M1849: 0.92)
switch opts.pattern
    case "roc_spec"     % ROC_TECH: sidelobes -29 dB beyond the main beam, -40 dB beyond 10 deg
        env_deg = [0 10 180];
        env_dB = [-29 -29 -40];
    case "nssl_measured"    % NSSL_1998 Section II: KOUN measurements
        env_deg = [0 2 6 20 180];
        env_dB = [-35 -35 -48 -55 -55];
end

az = (0:round(n_step))'*opts.step_deg;
az(end) = 360;
% Off-axis angle from the beam axis (azimuth 0, elevation tilt) to the
% station (azimuth az, elevation 0), spherical law of cosines.
if opts.tilt_deg == 0
    theta = abs(mod(az + 180, 360) - 180);  % exact on the horizon
else
    theta = acosd(min(max(cosd(opts.tilt_deg)*cosd(az), -1), 1));
end
g_rel = pencil_dB(theta, bw, env_deg, env_dB, opts.pattern);
slice = [az, g_max + g_rel];

% Angle where the main beam meets the envelope (first segment is flat).
x_deg = bw*sqrt(-env_dB(1)/12);

info.name = "WSR-88D (NEXRAD)";
info.g_max_dBi = g_max;
info.bw_deg = bw;
info.pattern = opts.pattern;
info.tilt_deg = opts.tilt_deg;
info.env_deg = env_deg;
info.env_dB = env_dB;
info.x_deg = x_deg;
info.peak_dBi = max(slice(:, 2));
info.floor_dBi = g_max + min(env_dB);
d_m = 8.53;         % ROC_TECH: reflector diameter 8.53 m (28 ft)
f_Hz = 2850e6;      % ROC_TECH: gain and beamwidth quoted at 2850 MHz
info.far_field_m = 2*d_m^2/(physconst("LightSpeed")/f_Hz);   % 2 D^2 / lambda, about 1.38 km
info.params = cell2table({ ...
    "Peak gain", g_max, "dBi", "ROC_TECH", "At 2850 MHz, including radome; M1849 Radar 1: 45 dBi"
    "3 dB beamwidth", bw, "deg", "ROC_TECH", "At 2850 MHz; M1849 Radar 1: 0.92 deg"
    "Main beam", -12, "dB per (theta/bw)^2", "M1851; M1849", "Parabolic main lobe, 3 dB at theta = bw/2"
    "Sidelobes to 10 deg (roc_spec)", -29, "dB", "ROC_TECH", "Specification, below peak"
    "Sidelobes beyond 10 deg (roc_spec)", -40, "dB", "ROC_TECH", "Specification, below peak, held to 180 deg"
    "First sidelobes to 2 deg (nssl_measured)", -35, "dB", "NSSL_1998", "KOUN, linear to -48 dB at 6 deg and -55 dB at 20 deg"
    "Far sidelobes (nssl_measured)", -55, "dB", "NSSL_1998", "Beyond 20 deg; no back-lobe data, held to 180 deg"
    "Beam elevation", opts.tilt_deg, "deg", "Worst case", "0 puts the peak on the horizon; lowest operational cut 0.5 deg"
    "Reflector diameter", d_m, "m", "ROC_TECH", "Far field 2D^2/lambda about 1.38 km at 2850 MHz"
    }, 'VariableNames', ["parameter", "value", "unit", "source", "note"]);
info.sources = sources_table(["ROC_TECH", "NSSL_1998", "M1849", "M1851"]);
end

% -------------------------------------------------------------------------
function g = pencil_dB(theta, bw, env_deg, env_dB, pattern)
% Main beam until it meets the envelope, then the envelope. Custom: no
% built-in antenna element has this envelope (phased.CustomAntennaElement
% would need the same table).
if pattern == "roc_spec"
    env = env_dB(1)*ones(size(theta));
    env(theta > env_deg(2) + 1e-9) = env_dB(3);    % step at 10 deg (10 deg itself is -29 dB)
else
    env = interp1(env_deg, env_dB, theta, "linear");
end
main = -12*(theta/bw).^2;
x = bw*sqrt(-env_dB(1)/12);                 % meets the flat first segment
g = main;
g(theta >= x) = env(theta >= x);
end

function T = sources_table(keys)
c = {
"ROC_TECH", "NOAA Radar Operations Center, ""NEXRAD Technical Information"" (WSR-88D antenna, transmitter, pedestal and VCP summary).", ...
    "https://roc.noaa.gov/public-documents/wsr88d/NEXRAD-Technical-Information.pdf"
"NSSL_1998", "NSSL, errata and supplement to the 1998 WSR-88D polarimetric upgrade report (D3a): KOUN antenna beamwidth and sidelobe measurements, Section II.", ...
    "https://www.nssl.noaa.gov/publications/wsr88d_reports/ERRATA1998Rpt_D3a(WC).pdf"
"M1849", "Recommendation ITU-R M.1849-3 (02/2023), Technical and operational aspects of ground-based meteorological radars; Annex 2 Table 6 (Radar 1) and Section 7.3.", ...
    "https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1849-3-202302-I!!PDF-E.pdf"
"M1851", "Recommendation ITU-R M.1851-2, Mathematical models for radiodetermination radar systems antenna patterns for use in interference analyses.", ...
    "https://www.itu.int/rec/R-REC-M.1851/en"
};
T = cell2table(c, 'VariableNames', ["key", "citation", "url"]);
T = T(ismember(T.key, keys), :);
end
