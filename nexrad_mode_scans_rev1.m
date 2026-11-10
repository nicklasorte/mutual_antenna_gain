function [scans, names] = nexrad_mode_scans_rev1(avset, opts)
%NEXRAD_MODE_SCANS_REV1 Scan structs for every WSR-88D scan mode in the sweep.
%   [scans, names] = nexrad_mode_scans_rev1() returns one nexrad_vcp_rev1
%   scan struct per WSR-88D scan mode in vcp_mode_cases_rev1 (every VCP in
%   operational Level II data, 2024-2026, with every option combination
%   the RDA allows, with every AVSET termination). Case 1 of vcp_mode_cases_rev1 (the working
%   paper's 3 rpm horizon scan with F.1245 patterns) is not a VCP and is
%   left out.
%
%   [scans, names] = nexrad_mode_scans_rev1("lowest") returns the 74 modes
%   of the Rev 5 study (AVSET at 6.4 deg only). See vcp_mode_cases_rev1.
%
%   nexrad_mode_scans_rev1(avset, low_el_deg=x, low_el_mode=m) gives every
%   mode the site's lowest elevation below 0.5 deg (nexrad_vcp_rev1), for
%   a per-site gain table.
%
%   scans  1-by-N cell array of scan structs, for gain_dist_table_rev1
%   names  1-by-N string array of mode labels ("VCP 212 SAILSx3 AVSET 6.4")
%
%   Example:
%     [scans, names] = nexrad_mode_scans_rev1();
%     T = gain_dist_table_rev1(ant_nexrad, scans, struct("mode_names", names));
%
%   See also vcp_mode_cases_rev1, nexrad_vcp_rev1, gain_dist_table_rev1.

arguments
    avset (1, 1) string = "all"
    opts.low_el_deg double = []
    opts.low_el_mode (1, 1) string = "add"
end
c = vcp_mode_cases_rev1(avset);
c = c(2:end);                                   % c{1} is the paper's constant-rate case
scans = cell(1, numel(c));
names = strings(1, numel(c));
for k = 1:numel(c)
    s = c{k};
    scans{k} = nexrad_vcp_rev1(s.vcp, sails=s.sails, mrle=s.mrle, ...
        avset_top_deg=s.avset_top_deg, cd_rate_deg_s=s.cd_rate_deg_s, ...
        low_el_deg=opts.low_el_deg, low_el_mode=opts.low_el_mode);
    names(k) = string(s.label);
end
end
