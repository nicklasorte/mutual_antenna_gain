function g_dBi = antenna_gain_rev1(ant_p, az_to_deg, el_to_deg, st)
%ANTENNA_GAIN_REV1 Radar antenna gain toward a direction, per pointing state.
%   g_dBi = antenna_gain_rev1(ant_p, az_to_deg, el_to_deg, st) returns the
%   gain (dBi) of the antenna prepared by antenna_prepare_rev1 toward the
%   direction (az_to_deg, el_to_deg), for each pointing state in st.
%
%   az_to_deg  azimuth of the other station, deg clockwise from north
%   el_to_deg  elevation of the other station, deg above local horizontal
%   st         struct of N-by-1 column vectors (one row per sample):
%      .mech_az   azimuth of the reflector boresight or array face normal
%      .tilt      mechanical beam elevation of a reflector (deg); ignored
%                 for arrays, which use ant_p.face_tilt_deg
%      .steer_az  array beam steering azimuth relative to the face normal
%      .steer_el  array beam steering elevation, absolute (deg)
%
%   Statgain (the default) uses the off-axis angle from the beam axis,
%   NTIA TM-13-489 Table 4-1; "f1245" uses the same angle with the
%   Rec. ITU-R F.1245 pattern. Aperture antennas combine azimuth and
%   elevation cuts in dB,
%   Rec. ITU-R M.1851-2 Eq. (15). Planar arrays use M.1851-2 Eqs (50)-(53).
%
%   See also antenna_prepare_rev1.

az_off = wrapTo180(az_to_deg - st.mech_az);

switch ant_p.type
    case "statgain"
        g_dBi = statgain_dBi(ant_p.sg, off_axis_deg(az_off, el_to_deg, st.tilt));

    case "f1245"
        g_dBi = f1245_dBi(ant_p.f, off_axis_deg(az_off, el_to_deg, st.tilt));

    case "rect"
        g_az = rect_cut_dB(ant_p.az_cut, az_off);
        if ant_p.el_cut.kind == "csc2"
            g_el = csc2_cut_dB(ant_p.el_cut, el_to_deg, st.tilt);
        else
            g_el = rect_cut_dB(ant_p.el_cut, el_to_deg - st.tilt);
        end
        g_dBi = ant_p.g_max_dBi + max(g_az + g_el, ant_p.floor_dB);

    case "parabolic"
        % Off-axis angle from the beam axis (spherical law of cosines);
        % the circular-aperture pattern is rotationally symmetric
        % (M.1851-2 Section 4.1).
        gam = off_axis_deg(az_off, el_to_deg, st.tilt);
        g_dBi = ant_p.g_max_dBi + max(parabolic_cut_dB(ant_p.cut, gam), ant_p.floor_dB);

    case "measured"
        g = ant_p.az_interp(az_off);
        if ~isempty(ant_p.el_interp)
            g = g + ant_p.el_interp(el_to_deg - st.tilt);
        end
        g_dBi = ant_p.g_max_dBi + max(g, ant_p.floor_dB);

    case "array_planar"
        g_dBi = planar_array_dBi(ant_p, az_to_deg, el_to_deg, st);

    otherwise
        error("radar_coupling:badAntennaType", "Unprepared antenna type ""%s"".", ant_p.type);
end
end

% -------------------------------------------------------------------------
function gam = off_axis_deg(az_off, el, tilt)
% Angle between the beam axis (azimuth 0, elevation tilt) and the
% direction (az_off, el), spherical law of cosines.
cos_g = cosd(el).*cosd(tilt).*cosd(az_off) + sind(el).*sind(tilt);
gam = acosd(min(max(cos_g, -1), 1));
end

% -------------------------------------------------------------------------
function g = statgain_dBi(sg, th)
% NTIA TM-13-489 Table 4-1. th is the off-axis angle (deg), 0..180.
G = sg.g;
g = zeros(size(th));
i_main = th <= sg.theta_m;
i_first = th > sg.theta_m & th <= sg.theta_r;
i_side = th > sg.theta_r & th <= sg.theta_b;
i_back = th > sg.theta_b;
g(i_main) = G - 4e-4*10^(G/10)*th(i_main).^2;
g(i_first) = 0.75*G - 7;
switch sg.category
    case 1                                      % Gmax >= 48 dBi
        g(i_side) = 29 - 25*log10(th(i_side));
        g(i_back) = -13;
    case 2                                      % 22 <= Gmax < 48 dBi
        g(i_side) = 53 - G/2 - 25*log10(th(i_side));
        g(i_back) = 11 - G/2;
    otherwise                                   % 10 <= Gmax < 22 dBi
        g(i_side) = 53 - G/2 - 25*log10(th(i_side));
        g(i_back) = 0;
end
end

% -------------------------------------------------------------------------
function g = f1245_dBi(f, phi)
% Rec. ITU-R F.1245 pattern. phi is the off-axis angle (deg), 0..180.
g = zeros(size(phi));
main = phi < f.phi_m;
g(main) = f.g - 2.5e-3*(f.dl*phi(main)).^2;
if f.large
    i = ~main & phi < max(f.phi_m, f.phi_r);
    g(i) = f.g1;
    i = phi >= max(f.phi_m, f.phi_r) & phi < 48;
    g(i) = 29 - 25*log10(phi(i));
    g(phi >= 48) = -13;
else
    i = ~main & phi < 48;
    g(i) = 39 - 5*log10(f.dl) - 25*log10(phi(i));
    g(phi >= 48) = -3 - 5*log10(f.dl);
end
end

% -------------------------------------------------------------------------
function g = rect_cut_dB(cut, th)
% Rectangular aperture cut, M.1851-2 Tables 4 and 6 (no pedestal).
a = abs(th);
mu = pi*cut.K*sind(min(a, 90))/cut.bw_deg;
g = 20*log10(abs(m1851_rect_field_rev1(cut.n, mu)));
switch cut.envelope
    case "peak"
        beyond = a > cut.brk_deg;
        g(beyond) = -cut.A*log(cut.B*a(beyond)/cut.bw_deg);
    case "average"
        beyond = a > cut.brk_deg;
        g(beyond) = -cut.A*log(cut.B*a(beyond)/cut.bw_deg) + cut.c_av;
    otherwise                                   % "theoretical"
        g(a > 90) = cut.floor_dB;
end
g = max(g, cut.floor_dB);
end

% -------------------------------------------------------------------------
function g = parabolic_cut_dB(cut, gam)
% Circular aperture, M.1851-2 Eq. (34) and Table 13 (no pedestal).
n = cut.n;
u = pi*cut.K*sind(min(gam, 90))/cut.bw_deg;
F = 2^(n+1)*factorial(n+1)*besselj(n+1, u)./u.^(n+1);
F(u == 0) = 1;
g = 20*log10(abs(F));
switch cut.envelope
    case "peak"
        beyond = gam > cut.brk_deg;
        g(beyond) = -cut.C*log10(gam(beyond)/cut.bw_deg) - cut.D;
    case "average"
        % M.1851-2 Section 4.2 step 9: average mask = peak mask - 4 dB.
        beyond = gam > cut.brk_deg;
        g(beyond) = -cut.C*log10(gam(beyond)/cut.bw_deg) - cut.D - 4;
    otherwise
        g(gam > 90) = cut.floor_dB;
end
g = max(g, cut.floor_dB);
end

% -------------------------------------------------------------------------
function g = csc2_cut_dB(cut, th, tilt)
% Ground-radar cosecant-squared elevation pattern, M.1851-2 Eqs (22),
% (25)-(27) and Table 7. th is the elevation angle (deg); tilt is the beam
% pointing elevation (deg), scalar or per sample.
bw = cut.bw_deg;
th = th + zeros(size(tilt));
if isnan(cut.start_deg)
    start = tilt + bw/2;
else
    start = cut.start_deg + zeros(size(tilt));
end
if any(start <= 0)
    error("radar_coupling:csc2Start", ...
        "CSC^2 start angle must be above 0 deg (got %.2f); check tilt and csc2_start_deg.", min(start));
end
th_null = tilt - bw/0.88;
g_uni = @(x, t) 20*log10(abs(sinc(50.8*sind(x - t)/bw)));
g = cut.floor_dB*ones(size(th));
i1 = th >= th_null & th <= start;
g(i1) = g_uni(th(i1), tilt(i1));
i2 = th > start & th <= cut.end_deg;
g(i2) = 20*log10(sind(start(i2))./sind(th(i2))) + g_uni(start(i2), tilt(i2));
g = max(g, cut.floor_dB);
end

% -------------------------------------------------------------------------
function g = planar_array_dBi(ant_p, az_to, el_to, st)
% Planar array gain, M.1851-2 Eqs (50)-(53), written with direction
% cosines in the array frame: u, v in the face plane, c along the normal.
% Face frame (ENU): normal n = [cb*sp, cb*cp, sb]; horizontal axis
% h = [cp, -sp, 0]; vertical axis w = h x n = [-sb*sp, -sb*cp, cb].
sp = sind(st.mech_az);
cp = cosd(st.mech_az);
sb = sind(ant_p.face_tilt_deg);
cb = cosd(ant_p.face_tilt_deg);

[e, nn, up] = enu_dir(az_to, el_to);
u = e.*cp - nn.*sp;
v = -e.*sb.*sp - nn.*sb.*cp + up.*cb;
c = e.*cb.*sp + nn.*cb.*cp + up.*sb;

[e0, n0, up0] = enu_dir(st.mech_az + st.steer_az, st.steer_el);
u0 = e0.*cp - n0.*sp;
v0 = -e0.*sb.*sp - n0.*sb.*cp + up0.*cb;

psi_x = wrap_pi(ant_p.kx*(u - u0));
psi_y = wrap_pi(ant_p.ky*(v - v0));
af2 = ant_p.afx(psi_x).*ant_p.afy(psi_y);       % |AF|^2 / sum(w.^2), both axes

front = c > 0;
g = ant_p.g_max_dBi + ant_p.front_back_dB + zeros(size(c));
elem_dB = 10*ant_p.elem_cos_exp*log10(c(front));
g(front) = ant_p.elem_gain_dBi + elem_dB + 10*log10(af2(front)) + ant_p.eff_dB;
g(front) = max(g(front), ant_p.g_max_dBi + ant_p.floor_dB);
end

% -------------------------------------------------------------------------
function [e, n, u] = enu_dir(az, el)
e = cosd(el).*sind(az);
n = cosd(el).*cosd(az);
u = sind(el);
end

% -------------------------------------------------------------------------
function p = wrap_pi(p)
% |AF|^2 is 2*pi periodic in psi; fold into [-pi, pi].
p = mod(p + pi, 2*pi) - pi;
end
