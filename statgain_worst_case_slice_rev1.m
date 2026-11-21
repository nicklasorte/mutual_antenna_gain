function [slice, info] = statgain_worst_case_slice_rev1(g_max_dBi, opts)
%STATGAIN_WORST_CASE_SLICE_REV1 Generic statgain worst-case horizontal gain slice.
%   slice = statgain_worst_case_slice_rev1(g_max_dBi) returns the gain of
%   the NTIA statistical gain (statgain) antenna model, NTIA TM-13-489
%   Table 4-1, toward a station on the horizon (0 deg elevation) as the
%   antenna turns through a full revolution, as an N-by-2 array:
%
%     slice(:, 1)  azimuth of the beam relative to the bearing of the
%                  station, deg, 0:step_deg:360 (0 and 360 = beam on the
%                  station)
%     slice(:, 2)  antenna gain toward the station, dBi
%
%   Statgain is a generic pattern: it needs only the peak gain, and it
%   is rotationally symmetric about the beam axis, so it treats every
%   antenna as a pencil beam. Its implied 3 dB beamwidth comes from the
%   main-beam term, 2 sqrt(3 / (4e-4 10^(G/10))) deg: about 0.92 deg at
%   45.5 dBi (the WSR-88D value) but about 3.4 deg at 34 dBi, more than
%   twice the ASR-11 azimuth beamwidth. Use it as a comparison with the
%   radar-specific patterns, or when no radar-specific pattern exists.
%
%   Pattern (theta = off-axis angle, deg; G = g_max_dBi; TM-13-489
%   Table 4-1):
%
%     main beam     0 <= theta <= theta_m:  G - 4e-4 10^(G/10) theta^2
%     first SL      theta_m < theta <= theta_r:  0.75 G - 7
%     sidelobes     theta_r < theta <= theta_b:
%                     29 - 25 log10(theta)              G >= 48
%                     53 - G/2 - 25 log10(theta)        G < 48
%     back          theta > theta_b:
%                     -13 (G >= 48), 11 - G/2 (22 <= G < 48), 0 (G < 22)
%
%     theta_m = 50 sqrt(0.25 G + 7) / 10^(G/20)
%     theta_r = 27.466 10^(-0.3 G/10) (G >= 48), 250 / 10^(G/20) (G < 48)
%     theta_b = 48 (G >= 22), 131.8257 10^(-G/50) (G < 22)
%
%   Worst case at the horizon: tilt_deg = 0 puts the beam peak on the
%   horizon; cos(theta) = cos(tilt) cos(az).
%
%   g_max_dBi  peak gain, dBi, 10 or more (TM-13-489 Section 4.4.2)
%
%   Name-value options:
%     step_deg  azimuth step (0.1); 360 must be a whole number of steps
%     tilt_deg  beam elevation above the horizon (0, the worst case)
%
%   info fields: name, g_max_dBi, category (1: G >= 48, 2: 22 <= G < 48,
%   3: G < 22), theta_m_deg, theta_r_deg, theta_b_deg, bw_implied_deg,
%   back_dBi, peak_dBi, params, sources.
%
%   Example:
%     s = statgain_worst_case_slice_rev1(45.5);
%     plot(s(:,1), s(:,2)); xlabel("Azimuth from station (deg)"); ylabel("Gain (dBi)")
%
%   DRAFT//PRE-DECISIONAL
%
%   See also nexrad_worst_case_slice_rev1, asr11_worst_case_slice_rev1,
%   asr9_worst_case_slice_rev1.

arguments
    g_max_dBi (1,1) double {mustBeGreaterThanOrEqual(g_max_dBi, 10)}
    opts.step_deg (1,1) double {mustBePositive} = 0.1
    opts.tilt_deg (1,1) double {mustBeInRange(opts.tilt_deg, -90, 90)} = 0
end
n_step = 360/opts.step_deg;
if abs(n_step - round(n_step)) > 1e-9
    error("statgain_worst_case_slice:badStep", "360 must be a whole number of steps; got step_deg = %g.", opts.step_deg);
end
G = g_max_dBi;
% TM-13-489 Table 4-1 breakpoints (deg) and category.
theta_m = 50*sqrt(0.25*G + 7)/10^(G/20);
if G >= 48
    category = 1;
    theta_r = 27.466*10^(-0.3*G/10);
    theta_b = 48;
    back = -13;
elseif G >= 22
    category = 2;
    theta_r = 250/10^(G/20);
    theta_b = 48;
    back = 11 - G/2;
else
    category = 3;
    theta_r = 250/10^(G/20);
    theta_b = 131.8257*10^(-G/50);
    back = 0;
end

az = (0:round(n_step))'*opts.step_deg;
az(end) = 360;
if opts.tilt_deg == 0
    th = abs(mod(az + 180, 360) - 180);     % exact on the horizon
else
    th = acosd(min(max(cosd(opts.tilt_deg)*cosd(az), -1), 1));
end
g = back*ones(size(th));
i_main = th <= theta_m;
i_first = th > theta_m & th <= theta_r;
i_side = th > theta_r & th <= theta_b;
g(i_main) = G - 4e-4*10^(G/10)*th(i_main).^2;
g(i_first) = 0.75*G - 7;
if category == 1
    g(i_side) = 29 - 25*log10(th(i_side));
else
    g(i_side) = 53 - G/2 - 25*log10(th(i_side));
end
slice = [az, g];

info.name = sprintf("Statgain %.4g dBi", G);
info.g_max_dBi = G;
info.category = category;
info.theta_m_deg = theta_m;
info.theta_r_deg = theta_r;
info.theta_b_deg = theta_b;
info.bw_implied_deg = 2*sqrt(3/(4e-4*10^(G/10)));
info.back_dBi = back;
info.peak_dBi = max(g);
info.params = cell2table({ ...
    "Peak gain", G, "dBi", "Input", ""
    "Main beam limit theta_m", theta_m, "deg", "TM_13_489", "Table 4-1"
    "First sidelobe limit theta_r", theta_r, "deg", "TM_13_489", "Table 4-1"
    "Sidelobe limit theta_b", theta_b, "deg", "TM_13_489", "Table 4-1"
    "First sidelobe level", 0.75*G - 7, "dBi", "TM_13_489", "Table 4-1"
    "Back level", back, "dBi", "TM_13_489", "Table 4-1"
    "Implied 3 dB beamwidth", info.bw_implied_deg, "deg", "Derived", "From the main-beam term"
    "Beam elevation", opts.tilt_deg, "deg", "Worst case", "0 puts the peak on the horizon"
    }, 'VariableNames', ["parameter", "value", "unit", "source", "note"]);
info.sources = cell2table({"TM_13_489", ...
    "NTIA Technical Memorandum TM-13-489, ""Antenna Models for Electromagnetic Compatibility Analyses,"" October 2012; Section 4.4.2 and Table 4-1 (statistical gain antenna model).", ...
    "https://www.ntia.gov/report/2012/antenna-models-electromagnetic-compatibility-analyses"}, ...
    'VariableNames', ["key", "citation", "url"]);
end
