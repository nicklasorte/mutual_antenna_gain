function files = mag_pair_kml_rev1(filename, pairs, R, opts)
%MAG_PAIR_KML_REV1 Write pair MAG results to a KMZ (or KML) for Google Earth.
%   files = mag_pair_kml_rev1(filename, pairs, R) writes one line per radar
%   pair, colored by its mutual antenna gain (MAG) at one percentile, and
%   one placemark per site. Each line's balloon lists the pair's takeoff
%   angles, every percentile in R, the worst mode pair and the largest MAG
%   on the path. The models are location-agnostic; the coordinates come
%   from the pair list.
%
%   filename  output file, ".kmz" (default if no extension) or ".kml". A
%             KMZ holds doc.kml, which links pairs.kml and sites.kml. A
%             ".kml" name writes <name>_pairs.kml and <name>_sites.kml.
%   pairs     table, one row per pair, in the row order of R:
%               lat_tx, lon_tx, lat_rx, lon_rx   site coordinates (deg)
%               name_tx, name_rx                 optional site names
%   R         mag_pair_cdf_rev1 output for the same pairs
%   opts      struct, all fields optional:
%     .pct         percentile that colors the lines (default 99.99; must be
%                  in R.pct_list)
%     .levels_dBi  color class edges, dBi (default [40 50 60 70]): five
%                  classes, below 40 to 70 and above
%     .title       document name (default "MAG pairs")
%
%   files  string array of the files written
%
%   Example:
%     R = mag_pair_cdf_rev1(T_nexrad, T_atc, el_tx, el_rx);
%     mag_pair_kml_rev1("pairs.kmz", pairs, R);
%
%   Marked DRAFT//PRE-DECISIONAL in the document name and every balloon.
%
%   See also mag_pair_cdf_rev1, kmlwrite, track2.

arguments
    filename (1, 1) string
    pairs table
    R (1, 1) struct
    opts (1, 1) struct = struct()
end
opts = fill_opts(opts);
need = ["lat_tx", "lon_tx", "lat_rx", "lon_rx"];
miss = setdiff(need, pairs.Properties.VariableNames);
if ~isempty(miss)
    error("mag_kml:missingColumns", "pairs needs column(s): %s.", strjoin(miss, ", "));
end
n = height(pairs);
if height(R.pct) ~= n
    error("mag_kml:sizeMismatch", "pairs has %d rows but R has %d pairs.", n, height(R.pct));
end
jp = find(abs(R.pct_list - opts.pct) < 1e-9, 1);
if isempty(jp)
    error("mag_kml:badPct", "opts.pct = %g is not in R.pct_list.", opts.pct);
end
names_tx = col_or(pairs, "name_tx", "Site " + string((1:n)') + "A");
names_rx = col_or(pairs, "name_rx", "Site " + string((1:n)') + "B");
mark = "DRAFT//PRE-DECISIONAL";

% ---- lines: great circle, 20 points each
ell = wgs84Ellipsoid("km");
lat = cell(n, 1);
lon = cell(n, 1);
d_km = distance(pairs.lat_tx, pairs.lon_tx, pairs.lat_rx, pairs.lon_rx, ell);
for i = 1:n
    [lat{i}, lon{i}] = track2(pairs.lat_tx(i), pairs.lon_tx(i), pairs.lat_rx(i), pairs.lon_rx(i), ell, "degrees", 20);
end
L = geoshape(lat, lon, "Geometry", "line");
pv = R.pct{:, 3:2 + numel(R.pct_list)};
v = pv(:, jp);
cls = discretize(v, [-Inf, opts.levels_dBi(:)', Inf]);
pal = palette(numel(opts.levels_dBi) + 1);
lname = names_tx + " - " + names_rx + sprintf(" (%g%%: ", opts.pct) + compose("%.1f", v) + " dBi)";
desc = strings(n, 1);
pct_lab = compose("%g%%", R.pct_list);
for i = 1:n
    rows = "<tr><td>" + pct_lab' + "</td><td>" + compose("%.2f", pv(i, :)') + " dBi</td></tr>";
    desc(i) = "<b>" + mark + "</b><br/>" + ...
        sprintf("Path %s (tx) to %s (rx), %.1f km<br/>", names_tx(i), names_rx(i), d_km(i)) + ...
        sprintf("Takeoff angles: %.3f deg at tx, %.3f deg at rx<br/>", R.pct.el_tx_deg(i), R.pct.el_rx_deg(i)) + ...
        "<table border=""1"">" + "<tr><th>Percentile</th><th>MAG</th></tr>" + strjoin(rows, "") + "</table>" + ...
        sprintf("Worst mode pair at %g%%: %s<br/>", opts.pct, R.worst_mode(i, jp)) + ...
        sprintf("Largest MAG on the path: %.2f dBi<br/>", R.pct.max_on_path_dBi(i)) + ...
        sprintf("Percentiles conservative by up to %.2f dB (grid rounding).", R.bias_dB);
end
lw = 1 + cls;                                    % heavier lines for higher classes

% ---- sites: unique by name and position
st = table([names_tx; names_rx], [pairs.lat_tx; pairs.lat_rx], [pairs.lon_tx; pairs.lon_rx], ...
    'VariableNames', ["name", "lat", "lon"]);
[~, iu] = unique(st(:, ["name", "lat", "lon"]), "rows", "stable");
st = st(iu, :);
P = geopoint(st.lat, st.lon);
sdesc = "<b>" + mark + "</b><br/>" + compose("%.5f, %.5f", st.lat, st.lon);

% ---- write
[fdir, base, ext] = fileparts(filename);
if ext == ""
    ext = ".kmz";
end
if fdir == ""
    fdir = string(pwd);
end
switch lower(ext)
    case ".kml"
        f_pairs = fullfile(fdir, base + "_pairs.kml");
        f_sites = fullfile(fdir, base + "_sites.kml");
        write_parts(f_pairs, f_sites);
        files = [f_pairs; f_sites];
    case ".kmz"
        tmp = tempname;
        mkdir(tmp);
        c = onCleanup(@() rmdir(tmp, "s"));
        write_parts(fullfile(tmp, "pairs.kml"), fullfile(tmp, "sites.kml"));
        write_doc(fullfile(tmp, "doc.kml"), opts.title + " " + mark, opts, pal);
        zfile = fullfile(tmp, "out.zip");
        % doc.kml first: viewers open the first KML in the archive.
        zip(zfile, ["doc.kml", "pairs.kml", "sites.kml"], tmp);
        out = fullfile(fdir, base + ".kmz");
        copyfile(zfile, out);
        files = out;
    otherwise
        error("mag_kml:badExt", "filename must end in .kmz or .kml.");
end

    function write_parts(fp, fs)
        kmlwrite(fp, L, "Name", cellstr(lname), "Description", cellstr(desc), ...
            "Color", pal(cls, :), "LineWidth", lw);
        kmlwrite(fs, P, "Name", cellstr(st.name), "Description", cellstr(sdesc));
    end
end

% =========================================================================
function write_doc(f, name, opts, pal)
% doc.kml of the KMZ: links the pairs and sites layers and a legend.
% Considered: kmlwrite, which writes features only (no NetworkLink).
edges = [-Inf, opts.levels_dBi(:)', Inf];
leg = strings(size(pal, 1), 1);
for k = 1:size(pal, 1)
    hex = sprintf("#%02X%02X%02X", round(255*pal(k, :)));
    if k == 1
        lab = sprintf("below %g dBi", edges(2));
    elseif k == size(pal, 1)
        lab = sprintf("%g dBi and above", edges(k));
    else
        lab = sprintf("%g to %g dBi", edges(k), edges(k + 1));
    end
    leg(k) = sprintf("<font color=""%s"">&#9632;</font> %s<br/>", hex, lab);
end
txt = [
    "<?xml version=""1.0"" encoding=""UTF-8""?>"
    "<kml xmlns=""http://www.opengis.net/kml/2.2"">"
    "<Document>"
    "<name>" + esc(name) + "</name>"
    "<description><![CDATA[" + sprintf("Line color: MAG at %g%%<br/>", opts.pct) + strjoin(leg, "") + "]]></description>"
    "<NetworkLink><name>Pairs</name><Link><href>pairs.kml</href></Link></NetworkLink>"
    "<NetworkLink><name>Sites</name><Link><href>sites.kml</href></Link></NetworkLink>"
    "</Document>"
    "</kml>"];
writelines(txt, f);
end

function s = esc(s)
s = replace(s, ["&", "<", ">"], ["&amp;", "&lt;", "&gt;"]);
end

% -------------------------------------------------------------------------
function pal = palette(n)
% Sequential, light to dark (parula), so order reads in grayscale too.
% Considered: the NTIA gray-plus-accent palette, which has no ordered
% sequence of five.
pal = flipud(parula(n + 1));
pal = pal(2:end, :);
pal = flipud(pal);
end

function v = col_or(T, name, dflt)
if ismember(name, T.Properties.VariableNames)
    v = string(T.(name));
else
    v = dflt;
end
end

function opts = fill_opts(opts)
d = struct("pct", 99.99, "levels_dBi", [40 50 60 70], "title", "MAG pairs");
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
        opts.(f{i}) = d.(f{i});
    end
end
extra = setdiff(fieldnames(opts), f);
if ~isempty(extra)
    error("mag_kml:unknownOption", "Unknown option(s): %s.", strjoin(extra, ", "));
end
opts.title = string(opts.title);
end
