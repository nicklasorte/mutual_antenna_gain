function scan = nexrad_vcp_rev1(vcp, opts)
%NEXRAD_VCP_REV1 WSR-88D volume coverage pattern as a mech_rotating scan.
%   scan = nexrad_vcp_rev1(vcp) returns the scan struct of a WSR-88D
%   (NEXRAD) volume coverage pattern for radar_mutual_gain_cdf_rev1. Each
%   elevation cut is one step of a mech_rotating volume scan: one full
%   revolution at the cut's azimuth rate and elevation.
%
%   scan = nexrad_vcp_rev1(vcp, Name=Value) adds the volume-scan options.
%
%   vcp  12, 35, 112, 212 or 215: the VCPs in operational use (ROC VCP
%        Info Sheet, February 2026; NEXRAD Technical Information, Build 23).
%
%   References:
%     [ROC_VCP_SHEET] NOAA ROC, "Convective VCP Decision Matrix" (VCP Info
%       Sheet), Feb 2026: VCPs in use, SAILS/MRLE/AVSET/MPDA definitions,
%       volume times.
%       https://www.roc.noaa.gov/public-documents/operations-branch/VCP_Info_Sheet.pdf
%     [ROC_TECH] NOAA ROC, "NEXRAD Technical Information": cuts per VCP,
%       elevation range, SAILS/MRLE/AVSET limits, pedestal limit 36 deg/s.
%       https://roc.noaa.gov/public-documents/wsr88d/NEXRAD-Technical-Information.pdf
%     [ROC_VCP112] NOAA ROC, "Theory and Concept of Operations for VCP 112,"
%       19 Mar 2019: Table 1 cut list, Table 3 volume times.
%       https://roc.noaa.gov/WSR88D/PublicDocs/NewTechnology/Theory_ConOps_VCP112.pdf
%     [ROC_VCP_PARADIGM] NOAA ROC, "New VCP Paradigm," Oct 2015: VCP 35
%       cut table (Figure 1).
%       https://www.weather.gov/media/roc/NewScience/New_VCP_Paradigm_Public_Oct_2015.pdf
%     [LEVEL2_AWS] NOAA NEXRAD Level II archive (s3://unidata-nexrad-level2):
%       https://www.ncei.noaa.gov/products/radar/next-generation-weather-radar
%       https://unidata-nexrad-level2.s3.amazonaws.com/2026/11/04/KDIX/KDIX20261104_154008_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/01/09/KDIX/KDIX20240109_195311_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/05/07/KLOT/KLOT20240507_171514_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/05/28/KFWS/KFWS20240528_000104_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/10/09/KTBW/KTBW20241009_035113_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/04/26/KOAX/KOAX20240426_150526_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2024/07/08/KHGX/KHGX20240708_013308_V06
%       https://unidata-nexrad-level2.s3.amazonaws.com/2025/02/16/KLWX/KLWX20250216_045233_V06
%     [METPY] Unidata MetPy 1.7, metpy.io.Level2File (VCP message decoder).
%       https://unidata.github.io/MetPy/
%
%   Name-value options:
%     sails          SAILS supplemental 0.5 deg scans inserted in the volume
%                    (default 0). Maximum per VCP (ROC): 35, 112, 215: 1;
%                    12, 212: 3.
%     mrle           MRLE: rescan the lowest 2, 3 or 4 elevations midway
%                    through the volume (default 0). VCPs 12, 212, 215.
%                    Cannot be used with SAILS.
%     avset_top_deg  AVSET termination elevation (default Inf, off). AVSET
%                    ends the volume after the cut at this elevation when no
%                    significant echoes lie above it; 6.4 deg is the lowest
%                    termination angle. VCPs 12, 112, 212, 215.
%     cd_rate_deg_s  azimuth rate of the Doppler cuts of the 0.5-1.3 deg
%                    split cuts (default: the decoded value in the table).
%                    It follows the site's Doppler PRF; decoded values are
%                    12.524, 14.458, 15.831 and 17.106 deg/s. VCPs 12, 35,
%                    212, 215.
%
%   Source of the cut tables: VCP messages (Message Type 5) in operational
%   Level II files (NOAA NEXRAD Level II on AWS, unidata-nexrad-level2),
%   decoded with MetPy 1.7 Level2File:
%     VCP 35   KDIX 2026-11-04 15:40Z
%     VCP 215  KDIX 2024-01-09 19:53Z
%     VCP 212  KLOT 2024-05-07 17:15Z (SAILSx1 removed for the base table)
%     VCP 12   KFWS 2024-05-28 00:01Z
%     VCP 112  KTBW 2024-10-09 03:51Z (SAILSx1 removed for the base table)
%   KOAX, KHGX and KLWX files give the same cuts and rates except the
%   split-cut Doppler rate (cd_rate_deg_s). Elevations are the nominal
%   angles; the files store them quantized (0.5 deg is 0.4834 deg).
%
%   Volume time is the sum of 360/rate over the cuts. Elevation
%   transitions are not modeled; the VCP 112 Theory and ConOps (2019)
%   Table 3 gives 322 s for the 20-cut volume against 320 s here.
%
%   SAILS and MRLE inserts go at the elevation-group boundary nearest
%   k/(n+1) of the base volume time, k = 1..n. This reproduces the decoded
%   VCP 212 SAILSx1 sequences (KLOT, KOAX, KHGX). The RDA places the VCP
%   112 SAILS cut one group later (after 1.8 deg, 13.6 s) than this rule.
%
%   scan fields:
%     .type           "mech_rotating"
%     .rate_deg_s     per-cut azimuth rate (deg/s, clockwise)
%     .tilt_deg       per-cut elevation (deg)
%     .n_rev          ones
%     .vcp            VCP number
%     .cuts           table: cut, el_deg, waveform, rate_deg_s, dur_s, kind
%                     ("base" | "sails" | "mrle")
%     .volume_s       volume time (s)
%     .options        options used
%
%   Waveforms: CS contiguous surveillance, CD contiguous Doppler with
%   range-ambiguity resolution (SZ-2), B batch, CDX contiguous Doppler
%   without ambiguity resolution.
%
%   Example:
%     scan = nexrad_vcp_rev1(212, sails=1);
%     sys = struct("tx", struct("type","f1245","g_max_dBi",46), "scan", scan);
%
%   See also radar_mutual_gain_cdf_rev1.

arguments
    vcp (1,1) double
    opts.sails (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.mrle (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.avset_top_deg (1,1) double = Inf
    opts.cd_rate_deg_s double = []
end

[cuts, lim] = vcp_table(vcp);

% ---- Options checks: fail loudly on combinations the RDA does not run.
if opts.sails > lim.sails
    error("radar_coupling:vcpSails", "VCP %d supports at most SAILSx%d (got %d).", vcp, lim.sails, opts.sails);
end
if opts.mrle ~= 0
    if ~lim.mrle
        error("radar_coupling:vcpMrle", "VCP %d does not support MRLE.", vcp);
    end
    if ~ismember(opts.mrle, [2 3 4])
        error("radar_coupling:vcpMrle", "mrle must be 0, 2, 3 or 4 (got %d).", opts.mrle);
    end
    if opts.sails > 0
        error("radar_coupling:vcpSailsMrle", "MRLE cannot be used with SAILS (ROC VCP Info Sheet).");
    end
end
if isfinite(opts.avset_top_deg)
    if ~lim.avset
        error("radar_coupling:vcpAvset", "VCP %d does not support AVSET.", vcp);
    end
    if opts.avset_top_deg < 6.4
        error("radar_coupling:vcpAvset", "AVSET terminates at or above 6.4 deg (got %.2f).", opts.avset_top_deg);
    end
end
if ~isempty(opts.cd_rate_deg_s)
    if vcp == 112
        error("radar_coupling:vcpCdRate", "cd_rate_deg_s is not defined for VCP 112 (two Doppler cuts per split elevation).");
    end
    if ~isscalar(opts.cd_rate_deg_s) || opts.cd_rate_deg_s <= 0 || opts.cd_rate_deg_s > 36
        error("radar_coupling:vcpCdRate", "cd_rate_deg_s must be a scalar in (0, 36] deg/s (pedestal limit 36 deg/s).");
    end
    cuts.rate_deg_s(cuts.waveform == "CD") = opts.cd_rate_deg_s;
end

% ---- AVSET: drop the cuts above the termination elevation.
cuts = cuts(cuts.el_deg <= opts.avset_top_deg + 1e-9, :);
cuts.dur_s = 360./cuts.rate_deg_s;
cuts.kind = repmat("base", height(cuts), 1);

% ---- SAILS (lowest elevation group) or MRLE (lowest mrle groups).
els = unique(cuts.el_deg, "stable");
if opts.sails > 0
    cuts = insert_groups(cuts, cuts(cuts.el_deg == els(1), :), opts.sails, "sails");
elseif opts.mrle > 0
    cuts = insert_groups(cuts, cuts(ismember(cuts.el_deg, els(1:opts.mrle)), :), 1, "mrle");
end
cuts.cut = (1:height(cuts))';
cuts = cuts(:, ["cut", "el_deg", "waveform", "rate_deg_s", "dur_s", "kind"]);

scan = struct();
scan.type = "mech_rotating";
scan.rate_deg_s = cuts.rate_deg_s;
scan.tilt_deg = cuts.el_deg;
scan.n_rev = ones(height(cuts), 1);
scan.vcp = vcp;
scan.cuts = cuts;
scan.volume_s = sum(cuts.dur_s);
scan.options = opts;
end

% =========================================================================
function cuts = insert_groups(cuts, grp, n, kind)
% Insert n copies of grp at the elevation-group boundaries nearest
% k/(n+1) of the base volume time.
grp.kind = repmat(string(kind), height(grp), 1);
t_end = cumsum(cuts.dur_s);
is_bnd = [cuts.el_deg(1:end-1) ~= cuts.el_deg(2:end); false];   % group ends, not volume end
bnd_idx = find(is_bnd);
pos = zeros(n, 1);
for k = 1:n
    [~, i] = min(abs(t_end(bnd_idx) - k*t_end(end)/(n + 1)));
    pos(k) = bnd_idx(i);
end
if numel(unique(pos)) < n
    error("radar_coupling:vcpInsert", "Not enough elevation groups to place %d %s inserts.", n, kind);
end
out = cuts([], :);
prev = 0;
for p = sort(pos)'
    out = [out; cuts(prev+1:p, :); grp]; %#ok<AGROW> at most 3 inserts
    prev = p;
end
cuts = [out; cuts(prev+1:end, :)];
end

% =========================================================================
function [cuts, lim] = vcp_table(vcp)
% Base cut tables (no SAILS or MRLE), decoded from Level II VCP messages;
% see the help text for the files. Rates in deg/s.
lo3 = [0.5 0.5 0.9 0.9 1.3 1.3]';          % split cuts: CS then CD
cs_cd = ["CS" "CD" "CS" "CD" "CS" "CD"]';
batch = [1.8 2.4 3.1 4.0 5.1 6.4]';
hi212 = [8.0 10.0 12.5 15.6 19.5]';
switch vcp
    case 35     % clear air, 0.5-6.4 deg, KDIX 2026-11-04
        el = [lo3; batch];
        wf = [cs_cd; repmat("B", 6, 1)];
        r = [4.966 14.458 4.966 14.458 5.471 14.458 ...
            15.491 17.754 16.930 18.073 18.073 18.073]';
        lim = struct("sails", 1, "mrle", false, "avset", false);
    case 215    % default precipitation, KDIX 2024-01-09
        el = [lo3; batch; 8.0; 10.0; 12.0; 14.0; 16.7; 19.5];
        wf = [cs_cd; repmat("B", 6, 1); repmat("CDX", 6, 1)];
        r = [11.470 15.831 13.370 15.831 15.919 15.831 ...
            16.776 20.654 19.534 20.237 20.237 20.237 ...
            24.862 25.642 25.642 25.642 25.642 25.642]';
        lim = struct("sails", 1, "mrle", true, "avset", true);
    case 212    % rapidly evolving convection, KLOT 2024-05-07
        el = [lo3; batch; hi212];
        wf = [cs_cd; repmat("B", 6, 1); repmat("CDX", 5, 1)];
        r = [21.149 12.524 21.149 12.524 23.027 12.524 ...
            26.389 27.334 28.224 26.400 26.400 26.400 ...
            28.422 28.411 28.740 28.740 28.740]';
        lim = struct("sails", 3, "mrle", true, "avset", true);
    case 12     % rapidly evolving convection, KFWS 2024-05-28
        el = [lo3; batch; hi212];
        wf = [cs_cd; repmat("B", 6, 1); repmat("CDX", 5, 1)];
        r = [21.149 24.994 21.149 24.994 23.027 24.994 ...
            25.719 25.939 26.741 27.598 27.664 27.620 ...
            28.400 28.806 28.488 28.488 28.488]';
        lim = struct("sails", 3, "mrle", true, "avset", true);
    case 112    % large-scale high-velocity systems (MPDA), KTBW 2024-10-09
        el = [repmat(0.5, 3, 1); repmat(0.9, 3, 1); repmat(1.3, 3, 1); batch; hi212];
        wf = [repmat(["CS" "CD" "CD"]', 3, 1); repmat("B", 6, 1); repmat("CDX", 5, 1)];
        r = [18.677 20.028 14.458 19.841 20.028 15.831 21.555 20.028 17.106 ...
            26.389 27.334 28.224 26.400 26.400 26.400 ...
            28.422 28.411 28.740 28.740 28.740]';
        lim = struct("sails", 1, "mrle", false, "avset", true);
    otherwise
        error("radar_coupling:badVcp", ...
            "VCP %g is not tabulated; use 12, 35, 112, 212 or 215.", vcp);
end
cuts = table(el, wf, r, 'VariableNames', ["el_deg", "waveform", "rate_deg_s"]);
end
