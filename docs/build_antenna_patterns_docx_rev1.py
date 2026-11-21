"""Build docs/antenna_patterns_working_paper_rev1.docx from the Markdown source.

DRAFT//PRE-DECISIONAL

The Markdown file is the source; this script turns it into a Word document
with Word figure and table captions, Word equations, and figure alt text
stored as the image description (not printed). Run from the docs folder:

    python3 build_antenna_patterns_docx_rev1.py

Needs pandoc (3.1 or later) on the path. The reference document
mutual_antenna_gain_methodology_report_rev5.docx supplies the styles and
the DRAFT//PRE-DECISIONAL header and footer.
"""
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

SRC = Path("antenna_patterns_working_paper_rev1.md")
OUT = Path("antenna_patterns_working_paper_rev1.docx")
REF = Path("mutual_antenna_gain_methodology_report_rev5.docx")

# Display equations: Markdown line -> TeX (pandoc writes Word equations).
EQUATIONS = {
    "G(θ) = 45.5 − 12 (θ / 0.925)<sup>2</sup> dBi,":
        r"$$G(\theta) = 45.5 - 12\left(\frac{\theta}{0.925^\circ}\right)^2 \ \text{dBi},$$",
    "cos θ = cos t × cos φ.":
        r"$$\cos\theta = \cos t \, \cos\varphi.$$",
    "G(φ) = 34 + max[g<sub>az</sub>(φ) + g<sub>el</sub>(0°), −30.5] dBi.":
        r"$$G(\varphi) = 34 + \max\left[g_{az}(\varphi) + g_{el}(0^\circ),\ -30.5\right] \ \text{dBi}.$$",
    "g<sub>az</sub> = 20 log<sub>10</sub> |F(μ)|, μ = π K sin φ / 1.45°, K = 68.8 (M.1851-2 Table 4),":
        r"$$g_{az} = 20\log_{10}\left|F(\mu)\right|, \qquad \mu = \frac{\pi K \sin\varphi}{1.45^\circ}, \qquad K = 68.8,$$"
        "\n\nwith K from M.1851-2 Table 4,",
    "g<sub>az</sub> = −17.51 ln(2.33 |φ| / 1.45°).":
        r"$$g_{az} = -17.51 \, \ln\left(\frac{2.33\,|\varphi|}{1.45^\circ}\right).$$",
}
INLINE = {
    "where F(μ) = [S(μ + π/2) + S(μ − π/2)] / (2 F(0)) and S(b) = sin(b)/b.":
        r"where $F(\mu) = \frac{S(\mu + \pi/2) + S(\mu - \pi/2)}{2F(0)}$ and $S(b) = \sin(b)/b$.",
}

IMG = re.compile(r"^!\[(?P<alt>[^\]]*)\]\((?P<path>[^)]+)\)\s*$")
FIGCAP = re.compile(r"^\*(?P<cap>Figure \d+\..*)\*\s*$")
TABCAP = re.compile(r"^\*(?P<cap>Table \d+\..*)\*\s*$")


def transform(text):
    """Return pandoc Markdown and the list of figure alt texts, in order."""
    lines = text.split("\n")
    out, alts = [], []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = IMG.match(line)
        if m:
            # Next non-blank line must be the italic figure caption.
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            c = FIGCAP.match(lines[j]) if j < len(lines) else None
            if not c:
                raise ValueError(f"Image on line {i + 1} has no figure caption")
            alts.append(m["alt"])
            out.append(f"![{c['cap']}]({m['path']}){{width=6.5in}}")
            i = j + 1
            continue
        t = TABCAP.match(line)
        if t:
            out.append(f"Table: {t['cap']}")
            i += 1
            continue
        if line.strip() in EQUATIONS:
            out.append(EQUATIONS[line.strip()])
            i += 1
            continue
        for k, v in INLINE.items():
            if k in line:
                line = line.replace(k, v)
        out.append(line)
        i += 1
    missing = [k for k in list(EQUATIONS) + list(INLINE) if k not in text]
    if missing:
        raise ValueError(f"Equation lines not found in the source: {missing}")
    return "\n".join(out), alts


def set_alt_text(docx, alts):
    """Store each figure's alt text as the image description."""
    tmp = Path(tempfile.mkdtemp())
    with zipfile.ZipFile(docx) as z:
        z.extractall(tmp)
    doc = tmp / "word" / "document.xml"
    xml = doc.read_text(encoding="utf-8")
    pics = list(re.finditer(r"<wp:docPr [^>]*?/>", xml))
    if len(pics) != len(alts):
        raise ValueError(f"{len(pics)} images in the document, {len(alts)} alt texts")
    for m, alt in zip(reversed(pics), reversed(alts)):
        tag = re.sub(r'descr="[^"]*"', f'descr="{escape(alt, {chr(34): "&quot;"})}"', m.group(0))
        xml = xml[:m.start()] + tag + xml[m.end():]
    # The picture element carries a second description; give it the same text.
    pics2 = list(re.finditer(r"<pic:cNvPr [^>]*?/>", xml))
    for m, alt in zip(reversed(pics2), reversed(alts)):
        tag = m.group(0)
        if "descr=" in tag:
            tag = re.sub(r'descr="[^"]*"', f'descr="{escape(alt, {chr(34): "&quot;"})}"', tag)
        else:
            tag = tag.replace("/>", f' descr="{escape(alt, {chr(34): "&quot;"})}"/>')
        xml = xml[:m.start()] + tag + xml[m.end():]
    doc.write_text(xml, encoding="utf-8")
    # Drop media files that no relationship uses (left over from the
    # reference document).
    used = set()
    for rels in (tmp / "word" / "_rels").glob("*.rels"):
        used |= set(re.findall(r'Target="media/([^"]+)"', rels.read_text(encoding="utf-8")))
    for p in (tmp / "word" / "media").glob("*"):
        if p.name not in used:
            p.unlink()
    out = Path(str(docx) + ".tmp")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(tmp.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(tmp).as_posix())
    shutil.move(out, docx)
    shutil.rmtree(tmp)


def main():
    md, alts = transform(SRC.read_text(encoding="utf-8"))
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False, encoding="utf-8") as f:
        f.write(md)
        tmp_md = f.name
    subprocess.run(["pandoc", tmp_md, "-f", "markdown", "-o", str(OUT),
                    f"--reference-doc={REF}", "--resource-path=."], check=True)
    Path(tmp_md).unlink()
    set_alt_text(OUT, alts)
    print(f"Wrote {OUT} with {len(alts)} figures")


if __name__ == "__main__":
    main()
