%EXAMPLE_ASR11_DRIFT_REV1 How long equal or multiple scan rates stay aligned.
%   At the 2.7 GHz Retuning LOE meeting of November 4, 2026, participants
%   noted that radars with the same scan rate can drift into relative
%   synchronization and couple main beam to main beam for long periods (up
%   to about a month), and that rates that are integer multiples (5 and
%   10 rpm) behave the same way. alignment_dwell_rev1 turns that into
%   numbers. Cases (ASR-11 values from nexrad_asr11_mag_rev1: 75 deg/s,
%   1.45 deg; rate tolerance +/-0.5 deg/s, M.2069-0, FAA value TODO):
%     1. ASR-11 and ASR-11, both 12.5 rpm (equal rates)
%     2. ASR-11 at 12.5 rpm and a radar at 6.25 rpm with the same beam
%        (2:1, the "multiples" case; illustrative, not a specific system)
%     3. ASR-11 and the WSR-88D 0.5 deg Doppler cut of VCP 212 (12.524
%        deg/s, 0.925 deg): no near-integer ratio, and the NEXRAD changes
%        rate every cut (at most 29 s), so alignment cannot persist
%     4. ASR-11 and a whole NEXRAD volume of VCP 212 SAILSx3 AVSET 6.4
%        (353.7 s). Each cut crosses the path once per volume, so the
%        volume acts as a slow rotation at 360/T_volume deg/s. When 74
%        ASR-11 revolutions nearly equal one volume (ASR-11 rate near
%        75.3 deg/s, inside the +/-0.5 deg/s tolerance), the same coupling
%        comes back once per volume. The effective NEXRAD window is the
%        0.925 deg beam crossed at the slowest 0.5 deg rate (12.524 deg/s),
%        in volume-phase degrees. Found in trial 61 of the 120 h run.
%   For each case: the nominal episode and the probability, over the drive
%   offsets, that an aligned episode outlasts 1 min to 30 days.

dw = struct("label", {}, "D", {});
dw(end+1) = struct("label", "ASR-11 / ASR-11, 12.5 rpm", "D", alignment_dwell_rev1(75, 75, 1.45, 1.45));
dw(end+1) = struct("label", "ASR-11 12.5 rpm / radar 6.25 rpm (2:1)", "D", alignment_dwell_rev1(75, 37.5, 1.45, 1.45));
dw(end+1) = struct("label", "ASR-11 / WSR-88D 0.5 deg Doppler cut", ...
    "D", alignment_dwell_rev1(75, 12.524, 1.45, 0.925, jitter_deg_s=[0.5 0], max_episode_s=360/12.524));
s_v = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
r_v = 360/s_v.volume_s;                         % volume phase rate (deg/s)
dw(end+1) = struct("label", sprintf("ASR-11 / WSR-88D volume (VCP 212 SAILSx3 AVSET 6.4, %.1f s)", s_v.volume_s), ...
    "D", alignment_dwell_rev1(75, r_v, 1.45, 0.925*r_v/12.524, jitter_deg_s=[0.5 0], max_ratio=100));
rows = cell(numel(dw), 1);
for k = 1:numel(dw)
    D = dw(k).D;
    pe = D.p_exceed.p_episode_longer';
    rows{k} = [{string(dw(k).label), sprintf("%d:%d", D.ratio), D.joint_period_s, D.window_deg, ...
        abs(D.delta0_deg_s), D.episode0_s}, num2cell(pe)];
end
Tdrift = cell2table(vertcat(rows{:}), "VariableNames", ["pair", "ratio", "joint_period_s", ...
    "window_deg", "nominal_drift_deg_s", "nominal_episode_s", "p_gt_1min", "p_gt_10min", ...
    "p_gt_1h", "p_gt_1d", "p_gt_30d"]);
disp(Tdrift)
