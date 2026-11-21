function [slice, info] = asr11_worst_case_slice_rev1(opts)
%ASR11_WORST_CASE_SLICE_REV1 ASR-11 worst-case horizontal gain slice.
%   slice = asr11_worst_case_slice_rev1() returns the ASR-11 antenna gain
%   (low beam) toward a station on the horizon (0 deg elevation) as the
%   antenna turns through a full revolution, as an N-by-2 array:
%
%     slice(:, 1)  azimuth of the beam relative to the bearing of the
%                  station, deg, 0:step_deg:360 (0 and 360 = beam on the
%                  station)
%     slice(:, 2)  antenna gain toward the station, dBi
%
%   Worst case at the horizon. The default beam elevation (tilt_deg = 0)
%   puts the peak of the elevation beam on the horizon, so the elevation
%   pattern adds 0 dB there. Any other tilt moves the horizon off the
%   peak (onto the lower edge of the beam for a positive tilt, or onto
%   the cosecant-squared skirt for a negative one), so tilt_deg = 0
%   bounds the gain toward any station on the horizon, as a coordination
%   zone that protects the radar needs.
%
%   Pattern (Rec. ITU-R M.1851-2 rectangular aperture, Eq. (15)):
%
%     G(az) = G_max + max(g_az(az) + g_el(0), floor)
%
%   g_az is the M.1851-2 peak envelope of the cosine taper (n = 1, chosen
%   by M.1851-2 Table 9 from the -24.5 dB first sidelobe):
%     main lobe: 20 log10 |F(mu)|, mu = pi K sin(az) / bw, K = 68.8,
%                F(mu) = [S(mu + pi/2) + S(mu - pi/2)] / 2 / F(0), S(b) = sin(b)/b
%     beyond the -14.4 dB point: -A ln(B |az| / bw), A = 17.51, B = 2.33
%   g_el is the M.1851-2 cosecant-squared pattern (Eqs (22), (25)-(27))
%   peaking at the tilt angle: the uniform-aperture beam up to tilt +
%   bw_el/2 (3 dB down), then the cosecant-squared skirt to 30 deg.
%   floor is the remote sidelobe level, +3.5 dBi (-30.5 dB below peak).
%
%   The peak envelope runs above the M.1464-2 first sidelobe (+9.5 dBi)
%   near the main beam, which is conservative.
%
%   Name-value options:
%     step_deg  azimuth step (0.1); 360 must be a whole number of steps
%     tilt_deg  beam elevation (0, the worst case); M.1464-2 Radar C
%               vertical scan range -2.5 to +2.5 deg, and the start of
%               the cosecant-squared skirt (tilt + 2.4 deg) must stay
%               above 0 deg
%
%   info fields: name, g_max_dBi, az_bw_deg, el_bw_deg, first_sl_dBi,
%   floor_dBi, n, K, A, B, brk_deg, g_el_dB (elevation gain at the
%   horizon), floor_az_deg (azimuth where the slice reaches the floor),
%   peak_dBi, el_fn (handle: g_el in dB versus elevation, deg), params,
%   sources.
%
%   Sources: M1464 (Table 1, Radar C: 34 dBi, 1.45 deg azimuth, 4.8 deg
%   elevation, +9.5 dBi first and +3.5 dBi remote sidelobes,
%   cosecant-squared to +30 deg), M1851 (pattern equations), FI_ASR11 and
%   WEBER_1999 (corroborating gain and beamwidth).
%
%   Example:
%     s = asr11_worst_case_slice_rev1();
%     plot(s(:,1), s(:,2)); xlabel("Azimuth from station (deg)"); ylabel("Gain (dBi)")
%
%   DRAFT//PRE-DECISIONAL
%
%   See also asr9_worst_case_slice_rev1, nexrad_worst_case_slice_rev1.

arguments
    opts.step_deg (1,1) double {mustBePositive} = 0.1
    opts.tilt_deg (1,1) double {mustBeInRange(opts.tilt_deg, -2.5, 2.5)} = 0
end
p.name = "ASR-11";
p.radar = "M.1464-2 Radar C";
p.g_max_dBi = 34;       % M1464 Radar C (FI_ASR11: > 33.5 dB; WEBER_1999: about 34 dB)
p.az_bw_deg = 1.45;     % M1464 Radar C (FI_ASR11: 1.40 deg)
p.el_bw_deg = 4.8;      % M1464 Radar C (WEBER_1999: 5 deg)
p.first_sl_dBi = 9.5;   % M1464 Radar C first sidelobe
p.remote_sl_dBi = 3.5;  % M1464 Radar C remote sidelobes
p.csc2_end_deg = 30;    % M1464 Radar C: cosecant-squared to +30 deg
p.remote_note = "M1464 Radar C remote sidelobe level";
p.remote_src = "M1464";
p.srcs = ["M1464", "M1851", "FI_ASR11", "WEBER_1999"];
[slice, info] = asr_slice(p, opts.step_deg, opts.tilt_deg);
end

% -------------------------------------------------------------------------
function [slice, info] = asr_slice(p, step_deg, tilt_deg)
n_step = 360/step_deg;
if abs(n_step - round(n_step)) > 1e-9
    error("asr_worst_case_slice:badStep", "360 must be a whole number of steps; got step_deg = %g.", step_deg);
end
if tilt_deg + p.el_bw_deg/2 <= 0
    error("asr_worst_case_slice:badTilt", "The cosecant-squared start (tilt + %.1f deg) must be above 0 deg; got tilt_deg = %g.", ...
        p.el_bw_deg/2, tilt_deg);
end
sll = p.first_sl_dBi - p.g_max_dBi;          % first sidelobe below peak
fl = p.remote_sl_dBi - p.g_max_dBi;          % floor below peak
% M.1851-2 Table 9: taper power from the first sidelobe level.
n = find(abs(sll) >= [13.2 20 30 39 45], 1, "last") - 1;
% M.1851-2 Table 4 beamwidth factor and Table 6 peak-envelope constants.
K = [50.8 68.8 83.2 95 106];
A = [8.584 17.51 26.882 35.84 45.88];
B = [2.876 2.33 1.962 1.756 1.56];
brk_pk = [-5.75 -14.4 -22.3 -31.5 -39.4];
cut_floor = [-30 -50 -60 -70 -80];
i = n + 1;
bw = p.az_bw_deg;
% Breakpoint: the main lobe falls to brk_pk between the peak and the
% first null (mu = (n+2) pi/2). fzero chosen: no closed form for n >= 1.
th_null = asind(min((n + 2)*pi/2*bw/(pi*K(i)), 1));
f = @(th) 20*log10(abs(rect_field(n, pi*K(i)*sind(th)/bw))) - brk_pk(i);
brk_deg = fzero(f, [1e-6*bw, th_null*(1 - 1e-9)]);

az = (0:round(n_step))'*step_deg;
az(end) = 360;
a = abs(mod(az + 180, 360) - 180);           % angle off the station, 0..180
g_az = 20*log10(abs(rect_field(n, pi*K(i)*sind(min(a, 90))/bw)));
beyond = a > brk_deg;
g_az(beyond) = -A(i)*log(B(i)*a(beyond)/bw);
g_az = max(g_az, cut_floor(i));
el_fn = @(el) csc2_dB(el, tilt_deg, p.el_bw_deg, p.csc2_end_deg, fl);
g_el = el_fn(0);                              % station on the horizon
slice = [az, p.g_max_dBi + max(g_az + g_el, fl)];

info.name = p.name;
info.radar = p.radar;
info.g_max_dBi = p.g_max_dBi;
info.az_bw_deg = bw;
info.el_bw_deg = p.el_bw_deg;
info.first_sl_dBi = p.first_sl_dBi;
info.floor_dBi = p.remote_sl_dBi;
info.tilt_deg = tilt_deg;
info.n = n;
info.K = K(i);
info.A = A(i);
info.B = B(i);
info.brk_deg = brk_deg;
info.g_el_dB = g_el;
% Envelope reaches the floor where -A ln(B az / bw) = floor (tilt 0).
info.floor_az_deg = bw/B(i)*exp((g_el - fl)/A(i));
info.peak_dBi = max(slice(:, 2));
info.el_fn = el_fn;
info.params = cell2table({ ...
    "Peak gain", p.g_max_dBi, "dBi", "M1464", p.radar
    "Azimuth 3 dB beamwidth", bw, "deg", "M1464", p.radar
    "Elevation 3 dB beamwidth", p.el_bw_deg, "deg", "M1464", p.radar
    "First sidelobe", p.first_sl_dBi, "dBi", "M1464", sprintf("%.1f dB below peak; sets the M.1851-2 taper", sll)
    "Remote sidelobes (floor)", p.remote_sl_dBi, "dBi", p.remote_src, p.remote_note
    "Taper power n", n, "", "M1851", "Table 9 (cosine taper for a -20 to -30 dB first sidelobe)"
    "K", K(i), "", "M1851", "Table 4"
    "A", A(i), "", "M1851", "Table 6, peak envelope"
    "B", B(i), "", "M1851", "Table 6, peak envelope"
    "Envelope breakpoint", brk_pk(i), "dB", "M1851", sprintf("Table 6; at %.4f deg", brk_deg)
    "Cosecant-squared upper limit", p.csc2_end_deg, "deg", "M1464", "M1851 Eqs (22), (25)-(27)"
    "Beam elevation", tilt_deg, "deg", "Worst case", "0 puts the elevation beam peak on the horizon"
    }, 'VariableNames', ["parameter", "value", "unit", "source", "note"]);
info.sources = sources_table(p.srcs);
end

function F = rect_field(n, mu)
% cos^n taper, M.1851-2 Table 4, normalized to F(0) = 1: cos^n expanded
% into exponentials, F = 2^-n sum C(n,k) S(mu + (n-2k) pi/2). Built-in sinc
% supplies S(b) = sinc(b/pi).
F = zeros(size(mu));
F0 = 0;
for k = 0:n
    F = F + nchoosek(n, k)*sinc((mu + (n - 2*k)*pi/2)/pi);
    F0 = F0 + nchoosek(n, k)*sinc((n - 2*k)/2);
end
F = F/F0;
end

function g = csc2_dB(el, tilt, bw, el_end, fl)
% Cosecant-squared elevation pattern, M.1851-2 Eqs (22), (25)-(27):
% uniform-aperture beam (peak at tilt) up to the start of the skirt,
% tilt + bw/2, then 20 log10(sin(start)/sin(el)) + g_uni(start) to
% el_end, floored.
nose = tilt + bw/2;
g_uni = @(x) 20*log10(abs(sinc(50.8*sind(x - tilt)/bw)));
g = fl*ones(size(el));
i1 = el >= tilt - bw/0.88 & el <= nose;
g(i1) = g_uni(el(i1));
i2 = el > nose & el <= el_end;
g(i2) = 20*log10(sind(nose)./sind(el(i2))) + g_uni(nose);
g = max(g, fl);
end

function T = sources_table(keys)
c = {
"M1464", "Recommendation ITU-R M.1464-2 (02/2015), Characteristics of radiolocation radars, and characteristics and protection criteria for sharing studies for aeronautical radionavigation and meteorological radars in the radiodetermination service operating in the frequency band 2 700-2 900 MHz; Table 1, Radar B (ASR-9) and Radar C (ASR-11).", ...
    "https://www.itu.int/dms_pubrec/itu-r/rec/m/R-REC-M.1464-2-201502-I!!PDF-E.pdf"
"M1851", "Recommendation ITU-R M.1851-2, Mathematical models for radiodetermination radar systems antenna patterns for use in interference analyses.", ...
    "https://www.itu.int/rec/R-REC-M.1851/en"
"FI_ASR11", "Forecast International, ""ASR-11 (DASR),"" Radar Forecast, August 2014 (gain > 33.5 dB, 1.40 deg azimuth beamwidth).", ...
    "https://www.forecastinternational.com/archive/disp_pdf.cfm?DACH_RECNO=1032"
"WEBER_1999", "M. E. Weber, ""FAA Surveillance Radar Data as a Complement to the WSR-88D Network,"" MIT Lincoln Laboratory ATC-266, 1999, Section 3.1 (ASR-9 and ASR-11: 1.4 deg azimuth, 5 deg elevation, about 34 dB).", ...
    "https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1999_ATC-266_WW-15318.pdf"
"WEBER_1986", "M. E. Weber, ""Assessment of ASR-9 Weather Channel Performance: Analysis and Simulation,"" MIT Lincoln Laboratory ATC-138, 1986, Section II.A (maximum azimuth sidelobes -28 dB).", ...
    "https://archive.ll.mit.edu/mission/aviation/publications/publication-files/atc-reports/Weber_1986_ATC-138_WW-15318.pdf"
"WIKI_ASR9", "Wikipedia, ""ASR-9"" (1.4 deg azimuth beamwidth, 5 deg elevation, cosecant-squared high and low beams).", ...
    "https://en.wikipedia.org/wiki/ASR-9"
};
T = cell2table(c, 'VariableNames', ["key", "citation", "url"]);
T = T(ismember(T.key, keys), :);
end
