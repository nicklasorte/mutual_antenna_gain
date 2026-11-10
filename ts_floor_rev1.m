function [f21, f12] = ts_floor_rev1(pn, pa, tilts_n, el_n, tilt_a, el_a, delta_n, delta_a)
%TS_FLOOR_REV1 Largest MAG while both beams are off the path (event floor).
%   [f21, f12] = ts_floor_rev1(pn, pa, tilts_n, el_n, tilt_a, el_a,
%   delta_n, delta_a) returns the largest mutual antenna gain (dBi) the
%   WSR-88D and ATC radar can reach while the NEXRAD beam azimuth is at
%   least delta_n deg off the bearing to the ATC radar and the ATC beam
%   azimuth at least delta_a deg off the bearing to the NEXRAD: outside
%   the windows of scan_event_windows_rev1 (use its delta_eff values).
%     f21 = max G_tx,ATC + max G_rx,NEXRAD   (path ATC -> NEXRAD)
%     f12 = max G_tx,NEXRAD + max G_rx,ATC
%   Each maximum is over the azimuths beyond delta (0.01 deg grid, plus the
%   edge itself) and, for the NEXRAD, over every tilt in its scan.
%
%   pn, pa      structs with fields tx and rx (antenna_prepare_rev1 output)
%   tilts_n     NEXRAD tilts (deg)
%   el_n, el_a  takeoff angles at the NEXRAD and the ATC radar (deg)
%   tilt_a      ATC beam tilt (deg)
%
%   See also scan_event_windows_rev1, mag_pair_timeseries_rev1.

x = unique([delta_n; (ceil(delta_n*100)/100:0.01:180)'; 180]);
y = unique([delta_a; (ceil(delta_a*100)/100:0.01:180)'; 180]);
x = [x; -x];
y = [y; -y];
tu = unique(tilts_n(:));
[gnt, gnr] = deal(-Inf);
for i = 1:numel(tu)
    st = struct("mech_az", -x, "tilt", tu(i) + 0*x, "steer_az", 0*x, "steer_el", 0*x);
    gnt = max(gnt, max(antenna_gain_rev1(pn.tx, 0, el_n, st)));
    gnr = max(gnr, max(antenna_gain_rev1(pn.rx, 0, el_n, st)));
end
st = struct("mech_az", 180 - y, "tilt", tilt_a + 0*y, "steer_az", 0*y, "steer_el", 0*y);
gat = max(antenna_gain_rev1(pa.tx, 180, el_a, st));
gar = max(antenna_gain_rev1(pa.rx, 180, el_a, st));
f21 = gat + gnr;
f12 = gnt + gar;
end
