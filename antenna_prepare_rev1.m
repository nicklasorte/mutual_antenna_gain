function ant_p = antenna_prepare_rev1(ant)
%ANTENNA_PREPARE_REV1 Validate a radar antenna model and precompute constants.
%   ant_p = antenna_prepare_rev1(ant) checks the antenna model struct ant
%   and returns ant_p, which antenna_gain_rev1 evaluates. Call it once per
%   antenna, outside any Monte Carlo loop.
%
%   Supported ant.type values and required fields. When ant.type is absent
%   or empty, the default is "statgain".
%
%   "statgain"  (default) Statistical gain model, NTIA TM-13-489 Section
%           4.4.2 and Table 4-1 (from ECAC-TN-85-023, the DoD Joint
%           Spectrum Center default for Gmax >= 10 dBi). Rotationally
%           symmetric; gain depends only on Gmax and the off-axis angle.
%       .g_max_dBi    peak gain (dBi), >= 10
%
%   "f1245"  Rec. ITU-R F.1245 reference pattern, as coded in the 2.7 GHz
%           study tool (AntennaModel_ITU_1245_1.cs) and used for the ASR and
%           WSR-88D antennas in the 2.7 GHz link budget working paper
%           (Table 4.4-1). Rotationally symmetric.
%       .g_max_dBi       peak gain (dBi)
%       .d_over_lambda   optional aperture diameter in wavelengths; default
%                        10^((g_max_dBi - 7.7)/20), the F.1245 estimate
%           For D/lambda <= 100:
%             G = Gmax - 2.5e-3*(D/lambda*phi)^2      0 <= phi < phi_m
%             G = 39 - 5log10(D/lambda) - 25log10(phi) phi_m <= phi < 48
%             G = -3 - 5log10(D/lambda)               48 <= phi <= 180
%           For D/lambda > 100: main beam to phi_m, G1 to max(phi_m, phi_r),
%           29 - 25log10(phi) to 48 deg, -13 dBi beyond, where
%           G1 = 2 + 15log10(D/lambda), phi_m = 20*sqrt(Gmax - G1)/(D/lambda),
%           phi_r = 12.02*(D/lambda)^-0.6.
%
%   "rect"  Rectangular aperture, Rec. ITU-R M.1851-2 Annex 1 Section 2.1
%           (models without pedestal, Tables 4, 6 and 9). Azimuth and
%           elevation cuts are combined in dB, M.1851-2 Eq. (15).
%       .g_max_dBi    peak gain (dBi)
%       .az_bw_deg    azimuth 3 dB beamwidth (deg)
%       .az_sll_dB    azimuth first sidelobe level below peak (negative dB);
%                     selects the cosine taper with M.1851-2 Table 9
%       .el_bw_deg    elevation 3 dB beamwidth (deg)
%       .el_pattern   "rect" (default) or "csc2"
%       .el_sll_dB    elevation first sidelobe level (when el_pattern="rect")
%       .csc2_end_deg     upper angle of the CSC^2 shape (when "csc2")
%       .csc2_floor_dB    CSC^2 floor G0 (when "csc2"), M.1851-2 Table 7
%       .csc2_start_deg   optional absolute CSC^2 start angle, M.1851-2
%                         Eq. (24); default tilt + el_bw_deg/2
%       .envelope     "peak" (single entry), "average" (aggregate), or
%                     "theoretical", M.1851-2 Section 2.1.3
%       .floor_dB     optional total floor below peak; default is the
%                     azimuth mask floor of M.1851-2 Table 6
%
%   "parabolic"  Circular aperture, M.1851-2 Section 4 (Tables 11, 13, 14),
%           rotationally symmetric about the beam axis.
%       .g_max_dBi, .az_bw_deg (3 dB beamwidth), .az_sll_dB, .envelope,
%       .floor_dB (optional; default M.1851-2 Table 13 floor)
%
%   "measured"  Digitized principal-plane cuts (M.1851-2 recommends 5).
%       .g_max_dBi
%       .az_deg, .az_gain_dB   azimuth cut, -180..180 deg, gain relative to
%                              peak (0 dB at boresight)
%       .el_deg, .el_gain_dB   optional elevation cut relative to boresight
%       .floor_dB              floor below peak
%
%   "array_planar"  Rectangular planar array, M.1851-2 Section 6.2,
%           Eqs (50)-(53), with an optional Taylor amplitude taper. Every
%           field has a default, so struct("type", "array_planar") gives a
%           41 dBi, horizon-pointing array.
%       .g_max_dBi        target peak gain (dBi) used to size a square array
%                         when nx is not given (default 41, M.1464-2
%                         Table 2 Radar F1). Give g_max_dBi or nx, not both.
%       .nx, .ny          elements along the horizontal and vertical axes
%                         (default: smallest square array that reaches
%                         g_max_dBi; ny defaults to nx)
%       .freq_MHz         frequency for d/lambda (default 2900, the centre
%                         of the 2700-3100 MHz study band)
%       .dx_m, .dy_m      element spacing (default lambda/2: no grating
%                         lobes at any scan angle, M.1851-2 Section 6.1)
%       .elem_gain_dBi    element peak gain (default 10log10(4*pi*dx*dy/
%                         lambda^2), the gain of one unit cell; 4.97 dBi
%                         at lambda/2 spacing)
%       .elem_cos_exp     element power pattern exponent m, gain ~ cos^m
%                         (default 1, the projected cell area)
%       .taper            "uniform" or "taylor" (default)
%       .taylor_sll_dB    Taylor sidelobe level, both axes (default -25 dB,
%                         M.1464-2 Radar F1 receive; use -17 dB for its
%                         transmit pattern)
%       .taylor_nbar      Taylor nbar (default 4, the taylorwin default)
%       .face_tilt_deg    face normal elevation (deg, default 0, horizon)
%       .eff_dB           array efficiency (dB, <= 0, default 0)
%       .floor_dB         far-sidelobe floor below peak (negative dB)
%       .front_back_dB    level below peak behind the face (negative dB)
%                         Both default to the NTIA TM-13-489 Table 4-1
%                         back level for the same peak gain, minus the
%                         peak (for example -9.5 - 41 = -50.5 dB at 41 dBi).
%       The returned g_max_dBi is the computed broadside gain.
%
%   See also antenna_gain_rev1, radar_mutual_gain_cdf_rev1.

arguments
    ant (1,1) struct
end

ant_p = ant;
if ~isfield(ant, "type") || isempty(ant.type)
    ant_p.type = "statgain";
else
    ant_p.type = string(ant.type);
end

switch ant_p.type
    case "statgain"
        must_have_field(ant, "g_max_dBi");
        ant_p.sg = prep_statgain(ant.g_max_dBi);

    case "f1245"
        must_have_field(ant, "g_max_dBi");
        ant_p.f = prep_f1245(ant.g_max_dBi, ...
            field_or_default(ant, "d_over_lambda", 10^((ant.g_max_dBi - 7.7)/20)));

    case "rect"
        must_have_field(ant, ["g_max_dBi", "az_bw_deg", "az_sll_dB", "el_bw_deg", "envelope"]);
        envelope = check_envelope(ant.envelope);
        ant_p.az_cut = prep_rect_cut(ant.az_bw_deg, ant.az_sll_dB, envelope);
        el_pattern = string(field_or_default(ant, "el_pattern", "rect"));
        switch el_pattern
            case "rect"
                must_have_field(ant, "el_sll_dB");
                ant_p.el_cut = prep_rect_cut(ant.el_bw_deg, ant.el_sll_dB, envelope);
            case "csc2"
                must_have_field(ant, ["csc2_end_deg", "csc2_floor_dB"]);
                ant_p.el_cut = struct( ...
                    "kind", "csc2", ...
                    "bw_deg", ant.el_bw_deg, ...
                    "end_deg", ant.csc2_end_deg, ...
                    "floor_dB", ant.csc2_floor_dB, ...
                    "start_deg", field_or_default(ant, "csc2_start_deg", NaN));
            otherwise
                error("radar_coupling:badElPattern", ...
                    "el_pattern must be ""rect"" or ""csc2"", got ""%s"".", el_pattern);
        end
        ant_p.floor_dB = field_or_default(ant, "floor_dB", ant_p.az_cut.floor_dB);

    case "parabolic"
        must_have_field(ant, ["g_max_dBi", "az_bw_deg", "az_sll_dB", "envelope"]);
        envelope = check_envelope(ant.envelope);
        ant_p.cut = prep_parabolic(ant.az_bw_deg, ant.az_sll_dB, envelope);
        ant_p.floor_dB = field_or_default(ant, "floor_dB", ant_p.cut.floor_dB);

    case "measured"
        must_have_field(ant, ["g_max_dBi", "az_deg", "az_gain_dB", "floor_dB"]);
        ant_p.az_interp = make_cut_interp(ant.az_deg, ant.az_gain_dB, "az");
        if isfield(ant, "el_deg") && ~isempty(ant.el_deg)
            must_have_field(ant, "el_gain_dB");
            ant_p.el_interp = make_cut_interp(ant.el_deg, ant.el_gain_dB, "el");
        else
            ant_p.el_interp = [];
        end

    case "array_planar"
        ant_p = prep_array(ant, ant_p);

    otherwise
        error("radar_coupling:badAntennaType", ...
            "ant.type must be ""statgain"", ""f1245"", ""rect"", ""parabolic"", ""measured"" or ""array_planar"", got ""%s"".", ant_p.type);
end
end

% -------------------------------------------------------------------------
function ant_p = prep_array(ant, ant_p)
% Planar array with defaults; see the "array_planar" help for sources.
freq_MHz = field_or_default(ant, "freq_MHz", 2900);
lambda_m = physconst("LightSpeed") / (freq_MHz * 1e6);
dx = field_or_default(ant, "dx_m", lambda_m/2);
dy = field_or_default(ant, "dy_m", lambda_m/2);
% Unit-cell element: gain 4*pi*dx*dy/lambda^2, cos(theta) power pattern
% (projected cell area).
elem_gain = field_or_default(ant, "elem_gain_dBi", 10*log10(4*pi*dx*dy/lambda_m^2));
elem_cos_exp = field_or_default(ant, "elem_cos_exp", 1);
eff_dB = field_or_default(ant, "eff_dB", 0);
if eff_dB > 0
    error("radar_coupling:badEfficiency", "eff_dB must be <= 0 (M.1851-2 Section 7).");
end
taper = string(field_or_default(ant, "taper", "taylor"));
sll = field_or_default(ant, "taylor_sll_dB", -25);
nbar = field_or_default(ant, "taylor_nbar", 4);
if ~ismember(taper, ["uniform", "taylor"])
    error("radar_coupling:badTaper", "taper must be ""uniform"" or ""taylor"", got ""%s"".", taper);
end
weights = @(n) taper_weights(taper, n, nbar, sll);

has_n = isfield(ant, "nx") && ~isempty(ant.nx);
if has_n
    if isfield(ant, "g_max_dBi") && ~isempty(ant.g_max_dBi)
        error("radar_coupling:arraySizeConflict", ...
            "Give either nx and ny or g_max_dBi for an array, not both.");
    end
    nx = ant.nx;
    ny = field_or_default(ant, "ny", nx);
else
    % Size a square array to reach the target peak gain.
    g_target = field_or_default(ant, "g_max_dBi", 41);
    n_need = 10^((g_target - elem_gain - eff_dB)/10);
    nx = max(2, ceil(sqrt(n_need)));
    while taper_n_eff(weights(nx), weights(nx)) < n_need
        nx = nx + 1;
    end
    ny = nx;
end
wx = weights(nx);
wy = weights(ny);
% Peak gain at broadside, M.1851-2 Eq. (50) with taper efficiency:
% D = De * |AF|^2 / sum(w.^2); at broadside |AF|^2 = sum(w)^2.
g_max = elem_gain + 10*log10(taper_n_eff(wx, wy)) + eff_dB;

% Far-sidelobe floor and back level: NTIA TM-13-489 Table 4-1 back level
% for the same peak gain (11 - G/2 for 22 <= G < 48 dBi).
if g_max >= 10
    back_dBi = statgain_back_dBi(g_max);
    floor_default = back_dBi - g_max;
else
    floor_default = [];
end
floor_dB = field_or_default(ant, "floor_dB", floor_default);
front_back_dB = field_or_default(ant, "front_back_dB", floor_default);
if isempty(floor_dB) || isempty(front_back_dB)
    error("radar_coupling:missingField", ...
        "Array peak gain %.1f dBi is below 10 dBi; give floor_dB and front_back_dB.", g_max);
end

ant_p.nx = nx;
ant_p.ny = ny;
ant_p.dx_m = dx;
ant_p.dy_m = dy;
ant_p.freq_MHz = freq_MHz;
ant_p.elem_gain_dBi = elem_gain;
ant_p.elem_cos_exp = elem_cos_exp;
ant_p.floor_dB = floor_dB;
ant_p.front_back_dB = front_back_dB;
ant_p.afx = prep_array_factor(wx);
ant_p.afy = prep_array_factor(wy);
ant_p.kx = 2*pi*dx/lambda_m;
ant_p.ky = 2*pi*dy/lambda_m;
ant_p.face_tilt_deg = field_or_default(ant, "face_tilt_deg", 0);
ant_p.eff_dB = eff_dB;
ant_p.g_max_dBi = g_max;
ant_p.lambda_m = lambda_m;
ant_p.taper = taper;
if max(dx, dy)/lambda_m > 0.5
    warning("radar_coupling:gratingLobes", ...
        "Element spacing %.3f lambda exceeds 0.5 lambda; grating lobes can appear at large scan angles (M.1851-2 Section 6.1).", ...
        max(dx, dy)/lambda_m);
end
end

function w = taper_weights(taper, n, nbar, sll)
if taper == "uniform"
    w = ones(n, 1);
else
    w = taylorwin(n, nbar, sll);
end
end

function n_eff = taper_n_eff(wx, wy)
n_eff = (sum(wx)^2/sum(wx.^2)) * (sum(wy)^2/sum(wy.^2));
end

function g = statgain_back_dBi(G)
% NTIA TM-13-489 Table 4-1 level beyond the last breakpoint.
if G >= 48
    g = -13;
elseif G >= 22
    g = 11 - G/2;
else
    g = 0;
end
end

% -------------------------------------------------------------------------
function f = prep_f1245(g, dl)
% Rec. ITU-R F.1245 constants, as in AntennaModel_ITU_1245_1.cs.
if ~(isscalar(dl) && dl > 0)
    error("radar_coupling:f1245Range", "d_over_lambda must be a positive scalar.");
end
f.g = g;
f.dl = dl;
f.g1 = 2 + 15*log10(dl);
if g <= f.g1
    error("radar_coupling:f1245Range", ...
        "F.1245 needs Gmax > G1 = 2 + 15log10(D/lambda) = %.2f dBi, got %.2f dBi.", f.g1, g);
end
f.phi_m = 20*sqrt(g - f.g1)/dl;
f.phi_r = 12.02*dl^-0.6;
f.large = dl > 100;
end

% -------------------------------------------------------------------------
function sg = prep_statgain(g)
% NTIA TM-13-489 Table 4-1 breakpoints (deg) and gain category.
if g < 10
    error("radar_coupling:statgainRange", ...
        "Statgain applies to Gmax >= 10 dBi (NTIA TM-13-489 Section 4.4.2), got %.2f dBi.", g);
end
sg.g = g;
sg.theta_m = 50*sqrt(0.25*g + 7)/10^(g/20);
if g >= 48
    sg.category = 1;
    sg.theta_r = 27.466*10^(-0.3*g/10);
    sg.theta_b = 48;
elseif g >= 22
    sg.category = 2;
    sg.theta_r = 250/10^(g/20);
    sg.theta_b = 48;
else
    sg.category = 3;
    sg.theta_r = 250/10^(g/20);
    sg.theta_b = 131.8257*10^(-g/50);
end
end

% -------------------------------------------------------------------------
function cut = prep_rect_cut(bw_deg, sll_dB, envelope)
% M.1851-2 Table 9: taper from first sidelobe level (models without pedestal).
a = abs(sll_dB);
if a < 13.2
    error("radar_coupling:badSLL", ...
        "Rectangular aperture first sidelobe must be <= -13.2 dB (M.1851-2 Table 9), got %.2f dB.", sll_dB);
end
n = find(a >= [13.2 20 30 39 45], 1, "last") - 1;     % 0..4
% M.1851-2 Table 4 beamwidth factors and Table 6 mask constants.
K       = [50.8 68.8 83.2 95 106];
A       = [8.584 17.51 26.882 35.84 45.88];
B       = [2.876 2.33 1.962 1.756 1.56];
brk_pk  = [-5.75 -14.4 -22.3 -31.5 -39.4];
brk_av  = [-12.16 -20.6 -29.0 -37.6 -42.5];
c_av    = [-3.72 -4.32 -4.6 -4.2 -2.61];
floors  = [-30 -50 -60 -70 -80];
i = n + 1;
cut = struct("kind", "rect", "n", n, "bw_deg", bw_deg, "K", K(i), ...
    "A", A(i), "B", B(i), "c_av", c_av(i), "floor_dB", floors(i), ...
    "envelope", envelope, "brk_deg", NaN);
if envelope ~= "theoretical"
    brk_level = brk_pk(i);
    if envelope == "average"
        brk_level = brk_av(i);
    end
    % First null of the cos^n taper pattern is at mu = (n+2)*pi/2.
    arg = (n + 2)*pi/2 * bw_deg / (pi*K(i));
    th_null = asind(min(arg, 1));
    % fzero chosen over a closed form: the main-lobe crossing has no
    % closed form for n >= 1.
    f = @(th) 20*log10(abs(m1851_rect_field_rev1(n, pi*K(i)*sind(th)/bw_deg))) - brk_level;
    cut.brk_deg = fzero(f, [1e-6*bw_deg, th_null*(1 - 1e-9)]);
end
end

% -------------------------------------------------------------------------
function cut = prep_parabolic(bw_deg, sll_dB, envelope)
% M.1851-2 Table 14: parabolic taper power from peak sidelobe level.
a = abs(sll_dB);
if a < 15
    error("radar_coupling:badSLL", ...
        "Circular aperture first sidelobe must be <= -15 dB (M.1851-2 Table 14), got %.2f dB.", sll_dB);
end
n = find(a >= [15 20 27 33 38], 1, "last") - 1;        % 0..4
% M.1851-2 Table 11 beamwidth factors and Table 13 mask constants.
K       = [58.2125 72.5938 84.0529 96.3142 108.2317];
C       = [28.9 49.0 69.13 89.0 108.8];
D       = [11.9 14.4 15.46 16.12 16.27];
brk_pk  = [0.8537 0.9893 1.13 1.2165 1.2835];
brk_av  = [1.051 1.161 1.273 1.339 1.3906];
floors  = [-35 -50 -60 -70 -80];
i = n + 1;
brk = brk_pk(i);
if envelope == "average"
    brk = brk_av(i);
end
cut = struct("kind", "parabolic", "n", n, "bw_deg", bw_deg, "K", K(i), ...
    "C", C(i), "D", D(i), "brk_deg", brk*bw_deg, "floor_dB", floors(i), ...
    "envelope", envelope);
end

% -------------------------------------------------------------------------
function af = prep_array_factor(w)
% |AF(psi)|^2 / sum(w.^2) on a periodic grid, from a zero-padded FFT of the
% element weights. fft is used instead of summing exponentials per sample;
% the interpolant then serves millions of samples cheaply.
n_fft = 2^nextpow2(max(2^16, 256*numel(w)));
p = abs(fft(w(:), n_fft)).^2 / sum(w.^2);
p = fftshift(p);
psi = (-n_fft/2:n_fft/2-1)' * (2*pi/n_fft);           % [-pi, pi)
psi = [psi; pi];
p = [p; p(1)];
af = griddedInterpolant(psi, p, "linear", "none");
end

% -------------------------------------------------------------------------
function gi = make_cut_interp(ang_deg, gain_dB, which_cut)
ang_deg = ang_deg(:);
gain_dB = gain_dB(:);
if numel(ang_deg) ~= numel(gain_dB) || numel(ang_deg) < 2
    error("radar_coupling:badMeasuredCut", "%s angle and gain vectors must match and have >= 2 points.", which_cut);
end
if ~issorted(ang_deg, "strictascend")
    error("radar_coupling:badMeasuredCut", "%s angles must be strictly increasing.", which_cut);
end
gi = griddedInterpolant(ang_deg, gain_dB, "linear", "nearest");
end

% -------------------------------------------------------------------------
function envelope = check_envelope(envelope)
envelope = string(envelope);
if ~ismember(envelope, ["peak", "average", "theoretical"])
    error("radar_coupling:badEnvelope", ...
        "envelope must be ""peak"", ""average"" or ""theoretical"", got ""%s"".", envelope);
end
end

% -------------------------------------------------------------------------
function v = field_or_default(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

% -------------------------------------------------------------------------
function must_have_field(s, names)
missing = names(~isfield(s, names));
if ~isempty(missing)
    error("radar_coupling:missingField", "Antenna model is missing field(s): %s.", strjoin(missing, ", "));
end
end
