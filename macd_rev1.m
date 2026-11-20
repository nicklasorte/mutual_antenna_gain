function out = macd_rev1(ant1, ant2, opts)
%MACD_REV1 MATLAB port of the C# Mutual Antenna Coupling Distribution (MACD).
%   out = macd_rev1(ant1, ant2) follows the C# MACD calculation step by
%   step, so its figures can walk a reader through what the RFS tool does:
%
%     1. Sample each antenna's pointing over its scan at a step of one
%        quarter of its 3 dB beamwidth, and take its gain toward the other
%        antenna at each pointing (AntennaGainDistributions.cs,
%        DetermineAntennaGainDistribution: bearingIncrement =
%        BeamWidth_deg / 4, numBearings = (int)(scanwidth / increment) + 1).
%     2. Put the gains in a histogram with 1 dB bins, from the backlobe gain
%        truncated to an integer up to the mainbeam gain (binValMin =
%        (int)backLobeGain_dBi, binWidth_dB = 1.0), one count per pointing.
%     3. For an antenna with several modes or elevation cuts, average the
%        mode histograms, each weighted by P(mode, el) = P_mode * P_el
%        (2.7 GHz link budget working paper, Appendix B, Eq. (B-1) and
%        Section B.7). The project's C# files don't hold this step; it's
%        ported from the appendix.
%     4. Convolve the two antennas' histograms to get the MAG histogram
%        (DistributionAnalysis.ConvolveHistograms).
%     5. Read the MAG at a percentage point (ComputeValueAtPercentagePoint).
%
%   DistributionAnalysis.cs, which holds the histogram update, the
%   convolution, and the percentage-point lookup, isn't in the project, so
%   three conventions here are assumptions, stated in out.assumptions: a
%   gain goes in the bin whose lower edge is at or below it (floor), a gain
%   outside the bin range goes in the end bin, and the percentage point is
%   given two ways (the upper edge of the bin where the cumulative
%   probability reaches it, and a linear interpolation inside that bin).
%
%   ant1, ant2  structs:
%     name              label
%     g_max_dBi         mainbeam gain
%     bw_deg            3 dB azimuth beamwidth (sets the sampling step)
%     gain_fn           @(rel_az_deg, tilt_deg) gain toward the other
%                       antenna, dBi, for pointing rel_az_deg from it with
%                       the beam at tilt_deg (vectorized in rel_az_deg)
%     backlobe_dBi      gain at 180 deg (sets the lowest bin)
%     scan_width_deg    azimuth scan width (360 for a full rotation)
%     mainbeam_az_deg   center of the scan (0)
%     bearing_deg       bearing to the other antenna (0; the other
%                       antenna's is this plus 180)
%     modes             table with tilt_deg and prob (one row, tilt 0,
%                       prob 1, if absent); prob is P_mode * P_el
%
%   opts (name-value):
%     step_div   divisor of the beamwidth for the sampling step (4, as C#)
%     bin_dB     histogram bin width (1, as C#)
%     pcts       percentage points (50 90 99 99.9 99.99 99.999)
%
%   out fields:
%     samples    1-by-2 cell: table of rel_az_deg (signed pointing from
%                the other antenna) and gain_dBi, first mode
%     hist       1-by-2 cell: table of bin_lo_dBi and prob (mode-averaged)
%     mode_hist  1-by-2 cell: matrix, one column of probabilities per mode
%     mag        table: bin_lo_dBi, prob, exceed (P(MAG >= bin_lo))
%     pct        table: pct, upper_edge_dBi, interp_dBi
%     mean_dB, sigma_dB   mean and standard deviation of the MAG bins
%                (bin centers), as DetermineMeanAndSigma_dB reports
%     n_samples  samples per mode, per antenna
%     assumptions  string array
%
%   Example:
%     a = struct("name","A","g_max_dBi",30,"bw_deg",2, ...
%         "gain_fn",@(az,t) max(30 - 12*(az/2).^2, 0),"backlobe_dBi",0, ...
%         "scan_width_deg",360);
%     o = macd_rev1(a, a);   % MAG histogram of two identical antennas
%
%   DRAFT//PRE-DECISIONAL
%
%   See also simple_mag_nexrad_asr11_rev1, report_figures_macd_rev1.

arguments
    ant1 (1,1) struct
    ant2 (1,1) struct
    opts.step_div (1,1) double {mustBePositive} = 4
    opts.bin_dB (1,1) double {mustBePositive} = 1
    opts.pcts (:,1) double {mustBeInRange(opts.pcts, 0, 100)} = [50 90 99 99.9 99.99 99.999]'
end
ants = {fill_defaults(ant1), fill_defaults(ant2)};
ants{2}.bearing_deg = norm_az(ants{1}.bearing_deg + 180);   % as GetParameters

w = opts.bin_dB;
out.samples = cell(1, 2);
out.hist = cell(1, 2);
out.mode_hist = cell(1, 2);
out.n_samples = zeros(1, 2);
lo = zeros(1, 2);
p = cell(1, 2);
for a = 1:2
    A = ants{a};
    % Step 2 grid: lowest bin at the backlobe gain truncated toward zero
    % ((int) in C#), enough bins to hold the mainbeam gain.
    lo(a) = fix(A.backlobe_dBi/w)*w;
    nb = floor((A.g_max_dBi - lo(a))/w) + 1;
    H = zeros(nb, height(A.modes));
    for m = 1:height(A.modes)
        % Step 1: pointings at a quarter-beamwidth step over the scan.
        inc = A.bw_deg/opts.step_div;
        n = floor(A.scan_width_deg/inc) + 1;
        brg = norm_az(A.mainbeam_az_deg - A.scan_width_deg/2 + (0:n - 1)'*inc);
        srel = norm_az(brg - A.bearing_deg);        % signed, for display
        g = A.gain_fn(abs(srel), A.modes.tilt_deg(m));
        if m == 1
            out.samples{a} = table(srel, g, 'VariableNames', ["rel_az_deg", "gain_dBi"]);
        end
        out.n_samples(a) = n;
        % Step 2: one count per pointing, floor binning, clamped to the range.
        k = min(max(floor((g - lo(a))/w + 1e-9) + 1, 1), nb);
        H(:, m) = accumarray(k, 1, [nb 1])/n;
    end
    % Step 3: mode average weighted by P(mode, el).
    pw = A.modes.prob/sum(A.modes.prob);
    p{a} = H*pw;
    out.mode_hist{a} = H;
    out.hist{a} = table(lo(a) + (0:nb - 1)'*w, p{a}, 'VariableNames', ["bin_lo_dBi", "prob"]);
end

% Step 4: convolution. Lower edges add, so MAG bin k starts at
% lo1 + lo2 + (k - 1) w.
pm = conv(p{1}, p{2});
x = lo(1) + lo(2) + (0:numel(pm) - 1)'*w;
ex = flip(cumsum(flip(pm)));
out.mag = table(x, pm, ex, 'VariableNames', ["bin_lo_dBi", "prob", "exceed"]);

% Step 5: percentage points.
cdf = cumsum(pm);
up = zeros(size(opts.pcts));
ip = zeros(size(opts.pcts));
for i = 1:numel(opts.pcts)
    q = opts.pcts(i)/100;
    k = find(cdf >= q - 1e-12, 1);
    up(i) = x(k) + w;
    below = 0;
    if k > 1
        below = cdf(k - 1);
    end
    ip(i) = x(k) + w*(q - below)/pm(k);
end
out.pct = table(opts.pcts, up, ip, 'VariableNames', ["pct", "upper_edge_dBi", "interp_dBi"]);
c = x + w/2;
out.mean_dB = sum(pm.*c);
out.sigma_dB = sqrt(sum(pm.*(c - out.mean_dB).^2));
out.assumptions = ["Gain goes in the bin whose lower edge is at or below it (floor); DistributionAnalysis.cs isn't in the project"; ...
    "Gains outside the bin range go in the end bins"; ...
    "Percentage point given as the bin upper edge and as a linear interpolation inside the bin"; ...
    "Mode averaging (Step 3) follows the link budget working paper, Appendix B, Eq. (B-1)"];
end

% -------------------------------------------------------------------------
function A = fill_defaults(A)
d = struct("name", "", "scan_width_deg", 360, "mainbeam_az_deg", 0, "bearing_deg", 0, ...
    "modes", table(0, 1, 'VariableNames', ["tilt_deg", "prob"]));
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(A, f{i})
        A.(f{i}) = d.(f{i});
    end
end
need = ["g_max_dBi", "bw_deg", "gain_fn", "backlobe_dBi"];
for f = need
    if ~isfield(A, f)
        error("macd_rev1:missingField", "Antenna struct needs field %s.", f);
    end
end
end

function a = norm_az(a)
% Azimuth to (-180, 180], as NormalizeAzimuth and DetermineRelativeAzimuth.
a = mod(a + 180, 360) - 180;
a(a == -180) = 180;
end
