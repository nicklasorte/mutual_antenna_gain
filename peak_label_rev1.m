function txt = peak_label_rev1(r, P, g_asr_max)
%PEAK_LABEL_REV1 Short label saying which beams make a per-revolution peak.
%   txt = peak_label_rev1(r, P, g_asr_max) labels one row r of
%   out.scan_max.asr_rev (nexrad_asr11_mag_rev1) for Figure 5 of the
%   methodology report. "Main beam" means within 3 dB of that antenna's
%   peak gain on the path: for the NEXRAD, its gain on the path at the
%   cut's tilt (from the plateau table P = out.scan_max.plateaus); for the
%   ASR-11, g_asr_max (its peak, 34 dBi).
%
%   Example:
%     txt = peak_label_rev1(out.scan_max.asr_rev(k, :), out.scan_max.plateaus, 34);
%
%   See also report_figures_rev4, nexrad_asr11_mag_rev1.

arguments
    r (1,:) table
    P table
    g_asr_max (1,1) double
end
g_n_mb = P.nexrad_gain_on_path_dBi(P.el_deg == r.nexrad_el_at_max_deg);
if isempty(g_n_mb)
    g_n_mb = NaN;                               % cut not in the plateau table: not main beam
end
n_mb = r.nexrad_gain_at_max_dBi >= g_n_mb - 3;
a_mb = r.asr11_gain_at_max_dBi >= g_asr_max - 3;
if n_mb && a_mb
    txt = sprintf("%.1f dBi: both main beams (NEXRAD %.1f deg cut)", r.max_g21_dBi, r.nexrad_el_at_max_deg);
elseif a_mb
    txt = sprintf("%.1f dBi: ASR-11 main beam + NEXRAD %.1f deg cut, beam %.1f deg off path", ...
        r.max_g21_dBi, r.nexrad_el_at_max_deg, abs(r.nexrad_az_off_at_max_deg));
elseif n_mb
    txt = sprintf("%.1f dBi: NEXRAD main beam (%.1f deg cut) + ASR-11 beam %.1f deg off path", ...
        r.max_g21_dBi, r.nexrad_el_at_max_deg, abs(r.asr11_az_off_at_max_deg));
else
    txt = sprintf("%.1f dBi: NEXRAD %.1f deg cut + ASR-11, both off axis", r.max_g21_dBi, r.nexrad_el_at_max_deg);
end
txt = string(txt);
end
