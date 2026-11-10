function c = vcp_mode_cases_rev1(avset)
%VCP_MODE_CASES_REV1 Case list for the WSR-88D scan-mode sweep.
%   c = vcp_mode_cases_rev1() returns a cell array of option structs for
%   nexrad_asr11_mag_rev1. Each struct has a label field plus the options
%   that set the case. c{1} is the 2.7 GHz working paper baseline (F.1245
%   patterns, NEXRAD at 3 rpm on the horizon). c{2:end} are every WSR-88D
%   VCP in operational Level II data, 2024-2026, with every option
%   combination the RDA allows.
%
%   c = vcp_mode_cases_rev1(avset) sets the AVSET terminations:
%     "all"     (default) every cut elevation from 6.4 deg up to the cut
%               below the top: 216 modes
%     "lowest"  6.4 deg only: the 74 modes of the Rev 5 study (pinned
%               references)
%
%   The limits come from the VCP table (nexrad_vcp_rev1, scan.limits). The
%   ROC VCP Info Sheet says AVSET is compatible with SAILS and MRLE, and
%   SAILS and MRLE cannot be combined. AVSET at its lowest termination
%   (6.4 deg) removes the most high cuts and so gives the 0.5 deg cuts
%   their largest time share. A termination higher up gives the cuts just
%   above 6.4 deg more time than any 74-mode case, which matters for
%   takeoff angles near those cuts (review R3: up to 2.6 dB at 99.9% near
%   7 deg), so "all" includes every termination. The slowest decoded
%   split-cut Doppler rate (12.524 deg/s) is added where a VCP's default is
%   faster: the longest dwell at low tilt.
%
%   Example:
%     c = vcp_mode_cases_rev1();
%     o = rmfield(c{2}, "label"); o.seed = 1;
%     r = nexrad_asr11_mag_rev1(o);
%
%   See also example_nexrad_asr11_rev1, example_asr9_rev1, nexrad_vcp_rev1.

arguments
    avset (1, 1) string {mustBeMember(avset, ["all", "lowest"])} = "all"
end
c = {};
c{end+1} = struct("label", "Paper: F.1245, 3 rpm, horizon", "nexrad_scan", "constant", ...
    "nexrad_pattern", "f1245", "asr_pattern", "f1245");
vcp_list = [35 31 32 34 215 212 12 112];
cd_slow = 12.524;
for vcp = vcp_list
    s0 = nexrad_vcp_rev1(vcp);
    lim = s0.limits;
    adapt = struct("label", "", "sails", 0, "mrle", 0);          % SAILS or MRLE (or neither)
    for k = 1:lim.sails
        adapt(end+1) = struct("label", sprintf(" SAILSx%d", k), "sails", k, "mrle", 0); %#ok<AGROW>
    end
    if lim.mrle
        for k = [2 3 4]
            adapt(end+1) = struct("label", sprintf(" MRLEx%d", k), "sails", 0, "mrle", k); %#ok<AGROW>
        end
    end
    avs = Inf;
    if lim.avset
        el = unique(s0.cuts.el_deg);
        tops = el(el >= 6.4 - 1e-9 & el < max(el) - 1e-9)';   % a top-cut termination is the full volume
        if avset == "lowest"
            tops = 6.4;
        end
        avs = [Inf tops];
    end
    cds = {[]};
    cd_default = s0.cuts.rate_deg_s(find(s0.cuts.waveform == "CD", 1));
    if lim.cd_rate && abs(cd_default - cd_slow) > 1e-9
        cds = {[], cd_slow};
    end
    for a = adapt
        for t = avs
            for d = cds
                lab = sprintf("VCP %d%s", vcp, a.label);
                if isfinite(t)
                    lab = lab + sprintf(" AVSET %.1f", t);
                end
                if ~isempty(d{1})
                    lab = lab + sprintf(", CD %.3f", d{1});
                end
                c{end+1} = struct("label", lab, "vcp", vcp, "sails", a.sails, "mrle", a.mrle, ...
                    "avset_top_deg", t, "cd_rate_deg_s", d{1}); %#ok<AGROW>
            end
        end
    end
end
end
