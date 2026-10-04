#!/usr/bin/env python3
"""Build results/report.html: a single self-contained browser page with the
final-region  paper plots, the comparisons with data and the FASER
predictions, PNGs embedded as data URIs.

The page is READ LOCALLY -- open results/report.html in a browser. Because it
is self-contained (no external CSS, JS, fonts or images) it works straight
from file://, and nothing constrains its size. A second, shareable copy is
written outside the repo at $BENCH_REPORT_COPY (see config.sh); that is the
one to hand to someone else, since report.html itself is gitignored.
"""
import base64
import io
import json
import os
import shutil
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = f"{BASE}/results"
OUT = f"{RESULTS}/report.html"
sys.path.insert(0, f"{BASE}/analysis")
import paths  # noqa: E402  -- config.sh is the one source of every path
import beams  # noqa: E402


# CC neutrino DIS and the FASER and data studies live in results_nu/.
RESULTS_NU = f"{BASE}/results_nu"


# NOTE: Sherpa and Herwig each return half the physical neutrino cross-section
# (they average over two initial lepton helicities when a massless neutrino has
# only one). This is corrected in analyze_nu.py as a named constant and is
# deliberately NOT surfaced on the page -- it is an under-the-hood fix, not a
# result. The evidence and the correction are documented in NEUTRINO_NOTES.md.


CSS = """
:root {
  --ground: #f6f6f3; --card: #fcfcfb; --ink: #1e2126; --muted: #5d6470;
  --line: #e2e2dc; --accent: #2a78d6; --accent-ink: #205da8;
  --warn: #c2601f; --warnbg: #fdf3e9;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground: #15171b; --card: #1e2126; --ink: #e7e8e4; --muted: #9aa1ab;
    --line: #2c2f35; --accent: #7db3ea; --accent-ink: #9cc5f0;
    --warn: #e08a4b; --warnbg: #2a2119;
  }
}
:root[data-theme="dark"] {
  --ground: #15171b; --card: #1e2126; --ink: #e7e8e4; --muted: #9aa1ab;
  --line: #2c2f35; --accent: #7db3ea; --accent-ink: #9cc5f0;
  --warn: #e08a4b; --warnbg: #2a2119;
}
* { box-sizing: border-box; }
body {
  background: var(--ground); color: var(--ink); margin: 0;
  font-family: Charter, "Bitstream Charter", "Sitka Text", Cambria, Georgia, serif;
  font-size: 17px; line-height: 1.55;
}
.wrap { max-width: min(2100px, 97vw); margin: 0 auto;
  padding: 2rem 1rem 3rem; }
.mono {
  font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
  font-variant-numeric: tabular-nums;
}
header { max-width: 122ch; }
.eyebrow {
  font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
  font-size: 0.72rem; letter-spacing: 0.14em; text-transform: uppercase;
  color: var(--accent-ink); margin: 0 0 0.9rem;
}
h1 { font-size: 2.1rem; line-height: 1.15; margin: 0 0 0.8rem;
     font-weight: 700; text-wrap: balance; }
.setup { color: var(--muted); margin: 0; }
h2 { font-size: 1.35rem; margin: 0 0 0.2rem; text-wrap: balance; }
.secnote { color: var(--muted); margin: 0; max-width: 122ch;
  text-align: justify; text-justify: inter-word;
  hyphens: auto; -webkit-hyphens: auto; }
section { margin-top: 3rem; }
/* TABLES ARE CENTRED LIKE THE FIGURES (user, 2026-08-28).  The wrapper is
   capped at the text measure and centred in it, and the table is then centred
   inside the wrapper, so a narrow table sits under the middle of the
   paragraph instead of hugging the left edge.
   THE CAP IS THE FULL 122ch, NOT the figures' 92ch, and that difference is
   deliberate: a figure is an image that scales, whereas a table has an
   intrinsic width, so squeezing the wide ones to 92ch would replace a
   readable table with a horizontal scrollbar.  `overflow-x: auto` still
   catches the genuinely too-wide ones.
   NOT `display: flex` with `justify-content: center` here, which is the
   obvious way to centre and is a known trap: once the table is wider than the
   wrapper, centring an overflowing flex item pushes its left edge out of the
   scrollable area, where no scrollbar can reach it.  A block wrapper with
   `margin-inline: auto` on the table centres it when it fits and degrades to
   ordinary left-anchored scrolling when it does not. */
.tablewrap { overflow-x: auto; margin: 1.4rem auto 0; max-width: 122ch; }
table { border-collapse: collapse; min-width: 34rem; margin-inline: auto; }
/* a small table INSIDE a prose block: no minimum width, left aligned with
   the paragraph it explains, so it reads as part of the sentence rather
   than as one of the page's own result tables */
table.mini { min-width: 0; margin-inline: 0; font-size: 0.94em; }
table.mini th, table.mini td { padding: 0.2rem 1.1rem 0.2rem 0; }
table.mini td + td, table.mini th + th { text-align: right; }
th, td { text-align: left; padding: 0.45rem 1.4rem 0.45rem 0;
         border-bottom: 1px solid var(--line); white-space: nowrap; }
th { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
     font-size: 0.72rem; letter-spacing: 0.1em; text-transform: uppercase;
     font-weight: 500; color: var(--muted); }
td.num { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
         font-variant-numeric: tabular-nums; font-size: 0.92rem; }
tr.analytic td:first-child { color: var(--accent-ink); }
/* a result the sample cannot support: shown, but never as a number */
tr.insufficient td { opacity: 0.68; }
.nostat { font-style: italic; letter-spacing: 0.01em; cursor: help;
          border-bottom: 1px dotted currentColor; }
tr.prelim td { color: var(--muted); }
tr.sechead td { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
                font-size: 0.72rem; letter-spacing: 0.1em; text-transform: uppercase;
                color: var(--accent-ink); padding-top: 0.9rem; border-bottom: 0; }
.tag { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
       font-size: 0.62rem; letter-spacing: 0.08em; text-transform: uppercase;
       border: 1px solid var(--line); border-radius: 3px;
       padding: 0.05rem 0.35rem; margin-left: 0.5rem; vertical-align: 0.08em; }
.tabs { display: flex; flex-wrap: wrap; gap: 0.4rem; margin: 1.3rem 0 0;
        border-bottom: 1px solid var(--line); }
.tabs button {
  font: inherit; font-size: 0.92rem; color: var(--muted); cursor: pointer;
  background: none; border: 0; border-bottom: 2px solid transparent;
  padding: 0.5rem 0.9rem; margin-bottom: -1px; border-radius: 4px 4px 0 0;
}
.tabs button:hover { color: var(--ink); background: var(--card); }
.tabs button[aria-selected="true"] {
  color: var(--accent-ink); border-bottom-color: var(--accent); font-weight: 600;
}
.tabs button:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
[role="tabpanel"] { padding-top: 1rem; }
.missing { color: var(--muted); max-width: 122ch; text-align: justify;
  hyphens: auto; border-left: 2px solid var(--line);
           padding: 0.1rem 0 0.1rem 1rem; margin: 0.6rem 0 1.4rem; font-style: italic; }
.tabs.top { gap: 0.1rem; margin: 2rem 0 0; border-bottom-width: 2px; }
.tabs.top button { font-size: 1.02rem; padding: 0.6rem 1.2rem; }
/* outermost level: which event selection */
.tabs.sel { gap: 0.35rem; margin: 2.2rem 0 0; border-bottom: 0; }
.tabs.sel button { font-size: 1.05rem; font-weight: 600; padding: 0.55rem 1.3rem;
                   border: 1px solid var(--line); border-radius: 999px;
                   background: var(--card); }
.tabs.sel button[aria-selected="true"] {
  background: var(--accent); border-color: var(--accent); color: #fff; }
.tabs.sel button[aria-selected="true"]:hover { background: var(--accent-ink); }
.selnote { color: var(--muted); max-width: 122ch; margin: 1.4rem 0 0;
  text-align: justify; hyphens: auto;
           border-left: 2px solid var(--line); padding: 0.1rem 0 0.1rem 1rem; }
.atpoint { max-width: 122ch; text-align: justify; hyphens: auto;
  border-left: 3px solid var(--accent);
           background: var(--card); padding: 0.85rem 1.1rem;
           margin: 1.2rem 0 0; border-radius: 0 4px 4px 0;
           color: var(--muted); }
.atpoint-tag { display: block; font-family: ui-monospace, "SF Mono", Menlo,
               Consolas, monospace; font-size: 0.68rem; letter-spacing: 0.1em;
               text-transform: uppercase; color: var(--accent-ink);
               margin-bottom: 0.45rem; }
.warn { max-width: 122ch; text-align: justify; hyphens: auto;
  border-left: 3px solid var(--warn); background: var(--warnbg);
        padding: 0.85rem 1.1rem; margin: 1.4rem 0 0; border-radius: 0 4px 4px 0; }
.warn .mono { font-size: 0.86em; }
[role="tabpanel"][hidden] { display: none; }
/* EVERY figure block is one centred column at 75% of the text measure.
   The default grid used to auto-fill columns of at least 460px, so most
   figures rendered two-up and packed from the left -- 126 of the 149 blocks
   on the page are this class, which is why restyling `.grid.one` alone
   changed almost nothing. */
.grid { display: grid; max-width: 122ch;
  grid-template-columns: minmax(0, 92ch); justify-content: center;
        gap: 1.1rem; margin-top: 1.4rem; }
/* Every figure block is ONE CENTRED COLUMN at 75% of the text measure.
   `.grid.one` is kept only so existing markup still matches; it adds
   nothing, so the two cannot drift apart.
   NOTE FOR THE NEXT EDITOR: this is a CSS string inside a Python file.
   A `#` comment here is NOT a comment -- CSS reads it as a selector and
   swallows the rule that follows, which is exactly how an earlier version
   of this block silently disabled itself. Use / * * / only. */
.grid.one { }
figure { margin: 0; background: var(--card); border: 1px solid var(--line);
         border-radius: 6px; padding: 0.7rem; }
figure button { display: block; width: 100%; border: 0; padding: 0;
                background: none; cursor: zoom-in; border-radius: 3px; }
figure button:focus-visible { outline: 2px solid var(--accent); outline-offset: 3px; }
figure img { display: block; width: 100%; height: auto;
             background: #fcfcfb; border-radius: 3px; }
figcaption { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
             font-size: 0.78rem; color: var(--muted); padding: 0.55rem 0.2rem 0.05rem; }
dialog { border: 0; border-radius: 8px; padding: 0.6rem; max-width: min(97vw, 1900px);
         background: #fcfcfb; }
dialog::backdrop { background: rgba(10, 12, 16, 0.72); }
dialog img { display: block; width: 100%; height: auto; }
footer { margin-top: 3.5rem; color: var(--muted); font-size: 0.85rem;
         border-top: 1px solid var(--line); padding-top: 1rem; }

.byline{margin:.55rem 0 .2rem;font-size:1.02rem;font-weight:600;color:#1e2430}
.byline .affil{display:block;font-weight:400;font-size:.86rem;color:#5d6470;
  margin-top:.18rem;max-width:86ch;line-height:1.45}
.auxnote{margin:.35rem 0 0;font-size:.88rem;color:#5d6470;max-width:86ch}
"""

JS = """
// each .tabs list is an independent group (hadron level, ME level)
for (const list of document.querySelectorAll('.tabs')) {
  const tabs = [...list.querySelectorAll('button')];
  const select = btn => {
    for (const t of tabs) {
      const on = t === btn;
      t.setAttribute('aria-selected', on);
      t.tabIndex = on ? 0 : -1;
      document.getElementById(t.getAttribute('aria-controls')).hidden = !on;
    }
  };
  tabs.forEach((t, i) => {
    t.addEventListener('click', () => select(t));
    t.addEventListener('keydown', e => {
      const d = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
      if (!d) return;
      e.preventDefault();
      const next = tabs[(i + d + tabs.length) % tabs.length];
      select(next); next.focus();
    });
  });
}

const dlg = document.getElementById('zoom');
const dimg = document.getElementById('zoomimg');
for (const btn of document.querySelectorAll('figure button')) {
  btn.addEventListener('click', () => {
    dimg.src = btn.querySelector('img').src;
    dimg.alt = btn.querySelector('img').alt;
    dlg.showModal();
  });
}
dlg.addEventListener('click', e => { if (e.target === dlg) dlg.close(); });
"""


# Figures are embedded as LOSSLESS WebP, re-encoded here at build time; the
# files on disk stay PNG, because the paper's LaTeX cannot read WebP.
#
# Measured on this report's own 143 figures: lossless WebP is ~37% of the PNG
# size, taking the page from 15.6 MB to about 5.8 MB with pixel-identical
# output.  Two things worth not re-deriving:
#   - LOSSY WebP is far WORSE here (q=95 gives 74%).  These are flat-colour
#     line plots, which lossless compresses well and lossy does not.
#   - method=6 buys one further point (35.8%) for 23x the encode time
#     (165 s vs 7 s for the full set).  Not worth it.
WEBP_OPTS = dict(lossless=True, quality=100, method=4)


def b64(png):
    """(mime_subtype, base64 payload) for one figure, WebP where possible."""
    try:
        from PIL import Image
    except ImportError:      # Pillow absent: embed the PNG unchanged
        with open(png, "rb") as f:
            return "png", base64.b64encode(f.read()).decode()

    im = Image.open(png)
    # matplotlib writes RGBA with a fully opaque alpha channel; dropping it is
    # lossless and saves a plane.  Checked rather than assumed, so a future
    # figure that really is transparent keeps its alpha.
    if im.mode == "RGBA" and im.getchannel("A").getextrema() == (255, 255):
        im = im.convert("RGB")
    buf = io.BytesIO()
    im.save(buf, "WEBP", **WEBP_OPTS)
    return "webp", base64.b64encode(buf.getvalue()).decode()


# Every figure is embedded ONCE and referenced by key, so a plot shown under
# more than one selection tab carries a single payload rather than one per tab.
# This began as a way under the 16 MB artifact ceiling; the report is read
# locally now and has no ceiling, but the page still loads faster for it.
IMAGES = {}


def _check_css(css):
    """Refuse to ship a stylesheet containing a Python-style comment.

    CSS has no `#` comment.  A line beginning with one is read as a SELECTOR,
    and the parser then swallows the rule that follows it -- so a stray `#`
    silently disables the next rule and the page renders with the old layout
    while the source says otherwise.  That happened here: five `#` lines
    inserted while editing this string disabled the figure-grid rule, and the
    only symptom was that a layout change "did nothing".
    """
    bad = [(i + 1, l.strip()) for i, l in enumerate(css.splitlines())
           if l.strip().startswith("#")]
    if bad:
        raise SystemExit(
            "CSS contains Python-style '#' comments, which disable the rule "
            "that follows them:\n"
            + "\n".join(f"  line {n}: {t[:70]}" for n, t in bad[:6])
            + "\nUse /* ... */ instead.")
    # every rule should balance; an unbalanced brace eats the rest of the file
    if css.count("{") != css.count("}"):
        raise SystemExit(f"CSS braces do not balance: {css.count('{')} open "
                         f"vs {css.count('}')} close")
    return css


def fig(prefix, key, cap, resdir=RESULTS):
    path = f"{resdir}/{prefix}{key}.png"
    # A missing figure is said so in place rather than embedded as a broken
    # image or allowed to abort the whole page.
    if not os.path.exists(path):
        return (f'<figure><p class="missing">no {prefix}{key} figure'
                f'</p><figcaption>{cap}</figcaption></figure>')
    ref = IMAGES.get(path)
    if ref is None:
        ref = f"i{len(IMAGES)}"
        IMAGES[path] = ref
    return (f'<figure><button type="button" aria-label="enlarge">'
            f'<img data-img="{ref}" alt="{key} comparison"></button>'
            f'<figcaption>{cap}</figcaption></figure>')


def image_script():
    """The one copy of each PNG, plus the loader that fills in every <img>.

    ALSO WRITES A MANIFEST of what was embedded.  The figures are converted to
    WebP on the way in, so nothing in the finished page identifies the PNG it
    came from -- and twice a figure has been built, committed and never
    rendered because nothing complained.  tools/check_report_complete.py reads
    this manifest and compares it against what is on disk.
    """
    payloads = {ref: b64(path) for path, ref in IMAGES.items()}
    try:
        with open(f"{RESULTS}/report_figures.json", "w") as _f:
            json.dump(sorted(os.path.relpath(p, BASE) for p in IMAGES),
                      _f, indent=1)
    except Exception as _exc:                                 # noqa: BLE001
        print(f"  (could not write the figure manifest: {_exc})")
    entries = ",\n".join(f'"{ref}":"data:image/{mime};base64,{data}"'
                          for ref, (mime, data) in payloads.items())
    return ("const IMG={\n" + entries + "\n};\n"
            "for (const el of document.querySelectorAll('img[data-img]')) {\n"
            "  el.src = IMG[el.dataset.img];\n}\n")


def tabs(views, aria, ns="", cls="", bodies=None):
    """A tab switcher: one panel per view, its content from `bodies`.

    `ns` namespaces the element ids so several tab groups can coexist on the
    page.
    """
    bar, panels = [], []
    for i, entry in enumerate(views):
        view, label = entry[:2]
        vid = f"{ns}{view}"
        on = "true" if i == 0 else "false"
        bar.append(f'<button role="tab" id="t-{vid}" aria-controls="p-{vid}" '
                   f'aria-selected="{on}" tabindex="{0 if i == 0 else -1}">'
                   f'{label}</button>')
        body = (bodies or {}).get(view, "")
        panels.append(
            f'<div role="tabpanel" id="p-{vid}" aria-labelledby="t-{vid}"'
            f'{"" if i == 0 else " hidden"}>{body}</div>')
    return (f'<div class="tabs {cls}" role="tablist" aria-label="{aria}">'
            f'{"".join(bar)}</div>{"".join(panels)}')


NOMAD_NOTE = (
    "<b>Dimuon production against a measurement.</b> NOMAD "
    "measured the fraction of charged-current interactions producing a "
    "second, opposite-sign muon, in neutrino-energy bins from 6 to 300 "
    "GeV on an iron target. A second muon in charged-current DIS comes "
    "from a semileptonic charm decay, so the ratio is a charm fraction "
    "times a branching ratio. The calculation is integrated over the region "
    "the measurement&rsquo;s own theory used, on an iron nucleon built from "
    "free-proton parton distributions.<br><br>"
    "<b>The shape failure was the charm threshold, and the figure now shows "
    "it</b> (user, 2026-09-08: &ldquo;are you sure you are computing the "
    "right observable?&rdquo;). Reading the reference implementation of "
    "<a href=\"https://arxiv.org/abs/2009.00014\">arXiv:2009.00014</a> "
    "against ours settled the observable: the integration region, the W "
    "propagator, the target-mass term in the coefficient and the effective "
    "branching ratio all agree, and the charm structure function really is "
    "the one with charm in the <i>final</i> state. What differs is the "
    "treatment of the charm mass.<br><br>"
    "<b>The benchmark&rsquo;s FONLL convention is nearly massless here.</b> "
    "FONLL is the massless result plus a <i>damped</i> mass correction, and "
    "over NOMAD&rsquo;s kinematics the damping factor is 0.03 at "
    "Q&sup2; = 2.7 GeV&sup2; and 0.19 at 4 &mdash; so the charm threshold is "
    "almost entirely switched off. Recomputing the same ratio with three "
    "numerators and one denominator makes the consequence plain:<br><br>"
    "&nbsp;&nbsp;massless charm &mdash; &chi;&sup2;/N = 216, theory/data 2.54 "
    "falling to 1.45<br>"
    "&nbsp;&nbsp;FONLL (this benchmark) &mdash; 126, 2.16 falling to 1.37<br>"
    "&nbsp;&nbsp;<b>massive, n<sub>f</sub> = 3 &mdash; 53, and 1.36 to "
    "1.32</b><br><br>"
    "<b>With the massive calculation the energy dependence is reproduced and "
    "what is left is a flat normalisation.</b> The ratio to data stops "
    "sloping and sits at 1.31 to 1.36 across sixteen to two hundred GeV; the "
    "shape was the threshold and nothing else. A constant third is not a "
    "shape problem, and the leading suspect is the strange density: NOMAD "
    "dimuons are the classic constraint on it, the reference calculation "
    "fits strangeness <i>to these very data</i>, and NNPDF4.0&rsquo;s "
    "strange sea is not fitted to them. The reference also runs at NNLO "
    "with target-mass corrections, an intrinsic-charm parton set and a "
    "Q&sup2; floor of 1 GeV&sup2;, none of which this calculation has.<br>"
    "<br>"
    "<b>Two suspects were ruled out rather than assumed away.</b> The "
    "Q&sup2; floor moves the ratio the wrong way &mdash; raising it "
    "<i>lowers</i> the prediction, so the reference&rsquo;s lower floor "
    "would make its own job harder &mdash; and the free-proton target was "
    "already corrected to iron, which was itself worth a factor of two.<br>"
    "<br>"
    "No PDF uncertainty is included; the &chi;&sup2; treats the points as "
    "independent, since the correlations were never released, and is an "
    "indication of size only; and the theory is evaluated at each bin "
    "centre, roughest in the first bin (6 to 22 GeV represented by 15.9).")


def figure_or_note_nu(prefix, key, cap, resdir=None):
    """A figure, or a note saying it has not been generated.

    Defaults to results_nu/, where the data and FASER studies live.
    """
    resdir = resdir or RESULTS_NU
    path = f"{resdir}/{prefix}{key}.png"
    if not os.path.exists(path):
        return (f'<p class="missing">{os.path.basename(path)} not yet '
                f'generated</p>')
    return fig(prefix, key, cap, resdir)


def paper_fig(path, cap):
    """Embed one paper figure by exact path.

    A paper plot names its own output file, so it is embedded by that path
    rather than through `fig()`'s prefix + key.
    """
    ref = IMAGES.get(path)
    if ref is None:
        ref = f"i{len(IMAGES)}"
        IMAGES[path] = ref
    return (f'<figure><button type="button" aria-label="enlarge">'
            f'<img data-img="{ref}" alt="{cap}"></button>'
            f'<figcaption>{cap}</figcaption></figure>')


def _load_paper_plots(subdir="paper_plots"):
    """Import every script in analysis/<subdir>/, in filename order.

    `subdir` is "paper_plots", the final-region figures (user,
    2026-09-13); the earlier directory analysis/paper_plots/ was deleted with the
    0.2 < y < 0.9 results on 2026-09-19.

    DISCOVERED, NOT LISTED.  A hand-maintained list of figures is a second
    place to update, and the one that gets forgotten; the directory is the
    list.  Each script names its own output, title, message and -- since
    2026-09-08 -- the SECTION of the write-up it belongs to.

    Returns [(filename, module_or_None, error_string_or_None)].
    """
    import importlib.util
    ppdir = f"{BASE}/analysis/{subdir}"
    if not os.path.isdir(ppdir):
        return []
    out = []
    for fn in sorted(f for f in os.listdir(ppdir)
                     if f.endswith(".py") and not f.startswith("_")):
        path = os.path.join(ppdir, fn)
        try:
            spec = importlib.util.spec_from_file_location(
                f"_{subdir}_{os.path.splitext(fn)[0]}", path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            out.append((fn, mod, None))
        except Exception as exc:                              # noqa: BLE001
            out.append((fn, None, f"{type(exc).__name__}: {exc}"))
    return out


def paper_plot_block(fn, mod, err, heading=True, subdir="paper_plots"):
    """One figure: its message, the figure itself, and the script that made it.

    Split out on 2026-09-08 so that one rendering serves every per-figure
    sub-tab.
    """
    slug = os.path.splitext(fn)[0]
    if mod is None:
        return (f'<h2>{slug}</h2><p class="missing">script failed to '
                f'import: {err}</p>')
    title = getattr(mod, "TITLE", slug)
    msg = getattr(mod, "MESSAGE", "").strip()
    resdir = RESULTS_NU if getattr(mod, "RESULTS", "results_nu") \
        == "results_nu" else RESULTS
    png = getattr(mod, "OUTPUT", "")
    out = [f'<h2>{title}</h2>'] if heading else []
    if msg:
        out.append('<p class="secnote"><b>What this figure shows.</b> '
                   + msg.replace("\n\n", "<br><br>") + '</p>')
    out.append('<div class="grid one">')
    full = os.path.join(resdir, png) if png else ""
    if full and os.path.exists(full):
        out.append(paper_fig(full, getattr(mod, "CAPTION", title)))
    else:
        out.append(f'<p class="missing">{png or "no OUTPUT declared"} '
                   f'not generated yet &mdash; run '
                   f'<code>analysis/{subdir}/{fn}</code>.</p>')
    out.append('</div>')
    # PROVENANCE, by user instruction: the reader is told which script made
    # the figure and where it lives.
    nclaims = len(getattr(mod, "CLAIMS", []) or [])
    out.append(
        f'<p class="secnote" style="font-size:0.85em">Script: '
        f'<code>analysis/{subdir}/{fn}</code>'
        + (f' &middot; {nclaims} checked claim'
           f'{"s" if nclaims != 1 else ""} '
           f'(<code>tools/check_paper_plots.py'
           f'{" --dir " + subdir if subdir != "paper_plots" else ""} '
           f'{slug}</code>)'
           if nclaims else
           ' &middot; <b>no checked claims yet</b>')
        + '</p>')
    return "".join(out)


PAPER_PLOTS_V2_NOTE = (
    "The paper plots, all of them in the <b>final benchmark region</b> "
    "(user, 2026-09-13): a <b>tungsten</b> target built from 74 protons and "
    "110 neutrons with NNPDF4.0 NNLO and no nuclear effects, "
    "<b>Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no cut on y</b>, on "
    "both currents, at 400, 700, 1000, 2000 and 4000 GeV. Cross-sections "
    "are per nucleon. Every sample behind these figures is generated inside "
    "the region, on the proton and on the neutron separately; the analytic "
    "references include target-mass corrections. <b>Each figure has its own "
    "script</b>, named beneath it, in <code>analysis/paper_plots/</code>; "
    "each figure&rsquo;s physics message is re-derived by executable claims "
    "from the result files it plots, checked by "
    "<code>tools/check_paper_plots.py</code>.")


def paper_plots_panel():
    """The paper plots in the final region (user, 2026-09-13), fed from
    analysis/paper_plots/.  The landing tab since 2026-09-19, when the the earlier production
    "Paper plots" tab was removed with the 0.2 < y < 0.9 results.  A figure
    keeps the number it had in earlier by declaring V1_NUMBER.
    """
    sub = "paper_plots"
    loaded = _load_paper_plots(sub)
    if not loaded:
        return (f'<h2>Paper plots</h2>'
                f'<p class="secnote">{PAPER_PLOTS_V2_NOTE}</p>'
                '<p class="missing">No figure yet.</p>')
    views, bodies = [], {}
    for fn, mod, err in loaded:
        n = getattr(mod, "V1_NUMBER", None) if mod else None
        key = f"v{n}" if n else os.path.splitext(fn)[0]
        views.append((key, f"Plot #{n}" if n else key, ""))
        bodies[key] = paper_plot_block(fn, mod, err, subdir=sub)
    return ('<h2>Paper plots</h2>'
            f'<p class="secnote">{PAPER_PLOTS_V2_NOTE}</p>'
            + tabs(views, "paper figure", ns="pp2-", cls="sel",
                   bodies=bodies))


def nomad_table():
    """The NOMAD points against the calculation, energy by energy."""
    p = f"{RESULTS_NU}/nomad_dimuon_fonll.json"
    if not os.path.exists(p):
        return ('<p class="missing">nomad_dimuon_fonll.json not yet '
                'generated &mdash; run analysis/nomad_dimuon.py</p>')
    with open(p) as f:
        d = json.load(f)
    out = ['<div class="tablewrap"><table><thead><tr>'
           '<th>E<sub>&nu;</sub> bin [GeV]</th><th>centre</th>'
           '<th>data [10<sup>&minus;3</sup>]</th>'
           '<th>YADISM NLO (FONLL) [10<sup>&minus;3</sup>]</th>'
           '<th>theory / data</th></tr></thead><tbody>']
    for r in d["points"]:
        err = (r["stat"] ** 2 + r["sys"] ** 2) ** 0.5
        out.append(
            f'<tr><td>{r["e_lo"]:g}&ndash;{r["e_hi"]:g}</td>'
            f'<td>{r["e"]:g}</td>'
            f'<td>{1e3*r["r"]:.3f} &plusmn; {1e3*err:.3f}</td>'
            f'<td>{1e3*r["r_th"]:.3f}</td>'
            f'<td>{r["r_th"]/r["r"]:.2f}</td></tr>')
    out.append('</tbody></table></div>')
    out.append(f'<p class="secnote">&chi;&sup2; per point '
               f'{d["chi2_per_point"]:.0f} over {d["ndat"]} points, '
               f'uncorrelated. At Q&sup2; &gt; 4 GeV&sup2; it is '
               f'{d["chi2_q2min_4"]/d["ndat"]:.0f} and at 6 GeV&sup2; '
               f'{d["chi2_q2min_6"]/d["ndat"]:.0f}, so the floor is not '
               f'what drives it.</p>')
    return "".join(out)


FASER_RATES_NOTE = (
    "<b>Applications of the benchmark to FASER&nu;.</b> The neutrino flux is "
    "the fast neutrino flux simulation of Kling and Nevay (<a "
    "href=\"https://arxiv.org/abs/2105.08270\">2105.08270</a>, <a "
    "href=\"https://arxiv.org/abs/2402.13318\">2402.13318</a>), vendored in "
    "<code>data/faser_flux/</code>. The target is the 1.1 t of Run-3 "
    "tungsten, 74 protons and 110 neutrons, behind the flux file&rsquo;s own "
    "aperture. Both studies on this tab are at 300 fb<sup>&minus;1</sup>, "
    "the luminosity at which the FASER dimuon estimate the first of them is "
    "set against is quoted.")


# The public note; the full one, which also compares with a FASER-internal
# estimate, is read from data/faser_internal/ where it exists (not distributed).
FASER_DIMUON_NOTE = (
    '<b>The opposite-sign dimuon signal at FASER</b>, as a cut flow at '
    '300 fb<sup>&minus;1</sup> built with this benchmark&rsquo;s '
    'machinery: the vendored FASER&nu; flux, 1.1 t of tungsten as 74 '
    'protons and 110 neutrons behind the flux&rsquo;s own aperture, and '
    'event samples on a ladder of beam energies from our GENIE (default '
    'tune, every CC process, Pythia 8 decays) and from POWHEG-V2 at NLO '
    'showered with Pythia 8. Every cut is applied to the generated muons; '
    'the steps that need the hadronic final state ride on the '
    'ladder&rsquo;s events, and GENIE is also run on the flux histogram '
    'directly, which reproduces the ladder to within 1.02 on every '
    'row.<br><br><b>POWHEG-V2 at NLO gives 1.6 times the GENIE dimuon '
    'rate, and it is the charm fraction.</b> Restricted to the same '
    'Q&sup2; &gt; 4 GeV&sup2; region the two generators agree on the CC '
    'rate to 0.4%, but POWHEG-V2 makes charm in 13.3% of its events '
    'against GENIE&rsquo;s 8.9%: massless charm at NLO against '
    'GENIE&rsquo;s slow-rescaled massive charm, the same difference the '
    'charm paper plots measure across the energy range. The steps after '
    'it agree: the charm-to-muon step is 9.31% against 9.04%, since both '
    'decay charm with Pythia 8, and the momentum-and-angle acceptances '
    'differ by under 10%. Under the spectrometer cuts, both muons above '
    '100 GeV and inside 25 mrad, the expected signal at 300 '
    'fb<sup>&minus;1</sup> is 7.9 events from GENIE and 13.1 from '
    'POWHEG-V2 for &nu;<sub>&mu;</sub>, 10.5 and 17.3 with the '
    'antineutrino included, a ratio of 1.64. The absolute rates, row by '
    'row, are the second figure.<br><br><b>The PDF dependence of the '
    'POWHEG-V2 rates</b> (third figure): every point of the ladder is '
    'reweighted a posteriori to every member of six PDF sets, on the '
    'proton and, through isospin-mirrored copies of each set, on the '
    'neutron, and each set&rsquo;s 68% band is LHAPDF&rsquo;s own '
    'prescription over its members; GRV98 is a single LO member and has '
    'none. The nominal member reproduces the unweighted cut flow to '
    '1.0000 on every row. On the CC row the sets agree to within 6.2% of '
    'one another (0.956 to 1.015 of NNPDF4.0, GRV98 the lowest) with '
    'NNPDF4.0&rsquo;s own band 0.5%; on the charm rows the same sets '
    'spread over 22% (CT18 at 0.78, ATLASpdf21 0.80, ABMP16 0.86, MSHT20 '
    '0.91) and GRV98 falls to 0.56, because charm production in CC DIS is '
    'the strange PDF. The spectrometer row, both muons above 100 GeV '
    'inside 25 mrad, reads 17.3 &plusmn; 0.7 events for NNPDF4.0, 14.3 '
    'for CT18, 17.1 for MSHT20, 15.7 for ATLASpdf21, 16.0 for ABMP16 and '
    '9.9 for GRV98, with the antineutrino included: the sets disagree by '
    '17% where the largest band is 15% and NNPDF4.0&rsquo;s is 4%, and '
    'the 7-point scale envelope on that row is &minus;3.6% to +4.0%. A '
    'dimuon measurement at FASER constrains the strange sea at the level '
    'where the global fits disagree.'
)
_DIMUON_NOTE_PRIVATE = os.path.join(BASE, "data", "faser_internal", "dimuon_report_note.html")
if os.path.exists(_DIMUON_NOTE_PRIVATE):
    with open(_DIMUON_NOTE_PRIVATE) as _f:
        FASER_DIMUON_NOTE = _f.read()


def dimuon_cutflow_section():
    """The FASER dimuon cut flow (user, 2026-09-07): the collaboration's own
    GENIE estimate of the opposite-sign dimuon signal at 300 fb^-1, against
    our GENIE and POWHEG-V2, from analysis/faser_dimuon_cutflow.py."""
    pc = f"{RESULTS_NU}/faser_dimuon_cutflow.json"
    if not os.path.exists(pc):
        return ('<h2>The dimuon cut flow, against the FASER estimate</h2>'
                '<p class="missing">not run yet &mdash; tools/genie_dimuon_ladder.sh, '
                'tools/powheg_v2_dimuon_ladder.sh, then '
                'analysis/faser_dimuon_cutflow.py rates</p>')
    with open(pc) as f:
        d = json.load(f)
    # the FASER-internal reference (not distributed): present only where
    # data/faser_internal/ is, and redacted from the public release's result
    _r = d.get("reference") or {}
    ref, RL = _r.get("rows"), _r.get("label", "")
    cols = [("genie_total", "GENIE, all CC"),
            ("genie_q2", "GENIE, Q&sup2; &gt; 4 GeV&sup2;"),
            ("powheg_q2", "POWHEG-V2, Q&sup2; &gt; 4 GeV&sup2;")]
    cols = [(c, l) for c, l in cols if c in d["columns"]]
    nu = "&nu;<sub>&mu;</sub>"
    nub = "&nu;<sub>&mu;</sub>&nbsp;+&nbsp;&nu;&#772;<sub>&mu;</sub>"
    title = ('The dimuon cut flow, against the FASER estimate' if ref
             else 'The dimuon cut flow')
    out = [f'<h2>{title}</h2>',
           f'<p class="secnote">{FASER_DIMUON_NOTE}</p>',
           '<div class="grid one">' + figure_or_note_nu(
               "cmp_faser_dimuon_cutflow", "",
               "Left: the cumulative efficiency of each step of the cut flow "
               f"for &nu;<sub>&mu;</sub>: {RL + ', ' if ref else ''}our GENIE (squares "
               "all CC, diamonds Q&sup2; &gt; 4 GeV&sup2;, crosses the same "
               "GENIE run on the flux histogram directly rather than through "
               "the energy ladder) and POWHEG-V2 at NLO; filled markers have "
               "both muons inside 25 mrad as the spectrometer tier requires, "
               "open markers drop the angular cut. Centre: the efficiencies "
               "that build the rate, against neutrino energy on tungsten, GENIE "
               "in red and POWHEG-V2 in blue. Right: where the events come "
               "from in the flux, for the CC total, the charm-to-muon step and "
               "the tightest spectrometer row, counted per bin of neutrino "
               "energy.") + '</div>',
           '<div class="grid two">' + figure_or_note_nu(
               "cmp_faser_dimuon_rates", "",
               "Absolute event rates at 300 fb<sup>&minus;1</sup>, row by row, "
               f"&nu;<sub>&mu;</sub> alone: {RL + ' in their fiducial volume, ' if ref else ''}"
               "our GENIE and POWHEG-V2 in the full 1.1 t (filled) and "
               "scaled to their first row (open), which puts the three on one "
               "volume.") + figure_or_note_nu(
               "cmp_faser_dimuon_pdf", "",
               "The same POWHEG-V2 rates under six PDF sets, "
               "&nu;<sub>&mu;</sub> + &nu;&#772;<sub>&mu;</sub>, with each "
               "set&rsquo;s own PDF uncertainty at 68% confidence level "
               "(GRV98 has no error members); lower panel: the ratio to "
               "NNPDF4.0.") + '</div>']
    # the cut-flow table
    hdr = ('<div class="tablewrap"><table><thead><tr><th>selection</th>'
           + (f'<th>{RL}</th>' if ref else ''))
    for c, l in cols:
        hdr += f'<th>{l}, {nu}</th><th>{l}, {nub}</th>'
    hdr += '</tr></thead><tbody>'
    out.append(hdr)

    def cell(r):
        return (f'<td class="num">{r["events"]:.1f} &plusmn; {r["err"]:.1f}'
                f'<br><small>{r["step_eff_pct"]:.3g}% / '
                f'{r["cum_eff_pct"]:.3g}%</small></td>')
    # THE "bench" STEP IS NOT SHOWN (user, 2026-09-19): it is the old
    # benchmark dimuon tier with its 0.2 < y < 0.9 window, and only the final
    # region survives.  The JSON still carries it; the table skips it.
    steps = [st for st in d["steps"] if st["key"] != "bench"]
    for st in steps:
        k, lab = st["key"], st["label"]
        lab = (lab.replace("->", "&rarr;").replace("theta", "&theta;")
               .replace("Q2", "Q&sup2;").replace("mu ", "&mu; ").replace(" mu", " &mu;"))
        row = f'<tr><td>{lab}</td>'
        if not ref:
            pass
        elif k in ref:
            nev, step, cum = ref[k]
            row += (f'<td class="num">{nev:g}<br><small>{step:g}% / '
                    f'{cum:g}%</small></td>')
        else:
            row += '<td class="num">&mdash;</td>'
        for c, _l in cols:
            row += cell(d["columns"][c]["14"][k]) + cell(d["columns"][c]["sum"][k])
        out.append(row + '</tr>')
    out.append('</tbody></table></div>'
               '<p class="secnote">Each cell: events at 300 fb<sup>&minus;1</sup> '
               'on 1.1 t of tungsten with the Monte Carlo error, then the step '
               'efficiency and the cumulative efficiency relative to the first row. '
               + (f'The {RL} first row is in an unstated fiducial volume, so '
                  'its events are not directly comparable; its efficiencies are. '
                  if ref else '') +
               'The '
               'momentum rows require both muons inside 25 mrad, as the '
               'spectrometer tier does; the rows without the angular cut show '
               'what that acceptance costs.</p>')
    # the flux-mode closure, GENIE
    fm = d.get("flux_mode", {})
    if fm:
        out.append('<div class="tablewrap"><table><thead><tr><th>step</th>')
        for c in fm:
            l = dict(cols).get(c, c)
            out.append(f'<th>{l}: ladder, {nub}</th><th>flux mode, {nub}</th>'
                       f'<th>ladder / flux</th>')
        out.append('</tr></thead><tbody>')
        for st in steps:
            k = st["key"]
            row = f'<tr><td>{st["label"].replace("->", "&rarr;").replace("theta", "&theta;").replace("Q2", "Q&sup2;")}</td>'
            for c, rec in fm.items():
                a = d["columns"][c]["sum"][k]
                b = rec["sum"][k]
                lf = rec.get("ladder_over_flux", {}).get(k)
                row += (f'<td class="num">{a["events"]:.1f} &plusmn; {a["err"]:.1f}</td>'
                        f'<td class="num">{b["events"]:.1f} &plusmn; {b["err"]:.1f}</td>'
                        f'<td class="num">{lf:.3f}</td>' if lf else
                        f'<td class="num">{a["events"]:.1f}</td><td class="num">{b["events"]:.1f}</td><td>&mdash;</td>')
            out.append(row + '</tr>')
        n_ev = f"{sum(next(iter(fm.values()))['n_events'].values()):,}".replace(",", "&thinsp;")
        out.append('</tbody></table></div>'
                   f'<p class="secnote">The flux-mode check: the same GENIE handed '
                   f'the flux histogram itself ({n_ev} events over the four '
                   f'beam&ndash;nucleon combinations), so that every efficiency is a '
                   f'plain count with no interpolation in energy. The ladder route, '
                   f'which is the one POWHEG-V2 can take, is what the table above '
                   f'uses.</p>')
    return "".join(out)


def pions_section():
    """Single-inclusive pion and kaon production at FASERnu (2026-09-07),
    from analysis/faser_pions.py."""
    p = f"{RESULTS_NU}/faser_pions.json"
    if not os.path.exists(p):
        return ('<h2>Single-inclusive pion production</h2>'
                '<p class="missing">not run yet &mdash; tools/faser_pions_passes.sh, '
                'then analysis/faser_pions.py rates</p>')
    return pions_section_body(p) + sidis_section()


# ---------------------------------------------------------- SIDIS at FASER
SIDIS_TAB_NOTE = (
    "<b>Everything this benchmark has on semi-inclusive DIS at FASER</b>, "
    "collected on one tab (user, 2026-09-21) so it can be handed to theory "
    "colleagues as a single reference. The process is "
    "&nu;N&nbsp;&rarr;&nbsp;&mu;hX with an identified charged pion or kaon, "
    "the observable semi-inclusive calculations deliver and fragmentation "
    "functions are fitted to.<br><br>"
    "<b>One region throughout</b>: Q<sup>2</sup> &gt; 4 GeV<sup>2</sup> and "
    "W &gt; 3 GeV, no cut on y and none on x, on a tungsten target of 74 "
    "protons and 110 neutrons with nuclear effects assumed to vanish. "
    "z = E<sub>h</sub>/&nu; with both energies in the <b>target rest "
    "frame</b>. Yields are quoted at 300 fb<sup>&minus;1</sup> and carry a "
    "60% hadron selection efficiency; multiplicities and the comparison "
    "against the NNLO calculation do not, because those are properties of "
    "the events rather than of the measurement.<br><br>"
    "<b>The neutrino current alone</b> (user, 2026-09-07), the one deliberate "
    "exception to carrying the two currents together; the muon-side spectra "
    "are extracted and on disk.")

SIDIS_EFF_NOTE = (
    "<b>What a calculation needs to become a FASER yield</b> (user, "
    "2026-09-21). A semi-inclusive calculation can be done in the benchmark "
    "region &mdash; Q<sup>2</sup> &gt; 4 GeV<sup>2</sup>, W &gt; 3 GeV and "
    "nothing else &mdash; but not with FASER&rsquo;s detector selection in "
    "it. These are the factors that bridge the two, bin by bin in z:"
    "<br><br><code>dN<sub>FASER</sub>/dz = "
    "&epsilon;(z) &times; dN<sub>region</sub>/dz</code><br><br>"
    "with the right-hand side counting <b>all</b> charged pions (or kaons) "
    "of the final state, no angular or energy cut on the hadron, in events "
    "passing the region cuts and no others.<br><br>"
    "<b>&epsilon; factorises into three pieces</b>, reported separately so "
    "that any one of them can be replaced: "
    "&epsilon;<sub>tier</sub>(z), the FASER Tier E event selection nested in "
    "the region (a scattered lepton above 200 GeV at more than 5 mrad, at "
    "least 5 charged tracks above 1 GeV inside tan&thinsp;&theta; &lt; 0.5 "
    "of which 4 inside 0.1, and &Delta;&phi; &gt; &pi;/2 between the lepton "
    "and the hadronic system); &epsilon;<sub>acc</sub>(z), the emulsion "
    "track acceptance on the hadron itself, tan&thinsp;&theta; &lt; 0.5 and "
    "E<sub>h</sub> &gt; 1 GeV; and a flat &epsilon;<sub>sel</sub> = 0.6 for "
    "selecting the hadron.<br><br>"
    "<b>&epsilon;<sub>tier</sub> depends on z although the cut is not made "
    "on the hadron.</b> Tier E asks for at least five charged tracks, and an "
    "event in which one hadron carries most of the hadronic energy has "
    "fewer of them, so the factor falls from 0.40 at z &rarr; 0 to 0.18 at "
    "z = 0.8 on the &nu;<sub>&mu;</sub> flux. The emulsion acceptance is the "
    "small piece: it is exactly one above z = 0.05, since a hadron carrying "
    "that much of a multi-hundred-GeV &nu; is always forward and always "
    "above 1 GeV, and only the lowest bin loses anything (0.88 for "
    "pions).<br><br>"
    "<b>Use the value at the beam energy, not the flux average, if the "
    "calculation is done at one energy.</b> Tier E needs a scattered lepton "
    "above 200 GeV, so its event efficiency runs (in POWHEG-V2) from 0.03 at 300 GeV to "
    "0.48 at 1 TeV and 0.57 at 2 TeV; the flux-averaged numbers in the table "
    "below are the ones to use for a calculation already folded over the "
    "FASER flux. The &nu;<sub>e</sub> column is larger throughout for the "
    "same reason: that flux is eight times smaller but considerably harder, "
    "a mean energy of 733 GeV against 312 GeV, so Tier E keeps 0.40 of it "
    "against 0.30.<br><br>"
    "<b>The factor is the mean of the five generators and the band their "
    "envelope</b>; the tables below quote the same mean. The spread is the model "
    "uncertainty to carry: &epsilon; is a ratio of two rates from the same "
    "sample, so the cross-section, the flux and the target cancel in it and "
    "what is left is the hadronisation model. It is under 10% of the value "
    "up to z &asymp; 0.5 and widens beyond, where the region&rsquo;s own "
    "yield is running out.<br><br>"
    "<b>Not included in &epsilon;</b>: nuclear corrections (absent on both "
    "sides of the ratio), detector smearing in z, and any vertex- or "
    "lepton-identification efficiency beyond the Tier E kinematics above.")


def sidis_efficiency_section():
    """The Tier E efficiency factors in z (user, 2026-09-21), from
    analysis/faser_sidis_efficiency.py."""
    p = f"{RESULTS_NU}/faser_sidis_efficiency.json"
    if not os.path.exists(p):
        return ('<h2>Tier E efficiency factors</h2>'
                '<p class="missing">not run yet &mdash; '
                'analysis/faser_sidis_efficiency.py</p>')
    with open(p) as f:
        d = json.load(f)
    edges = d["z_edges"]
    out = ['<h2>Tier E efficiency factors in z</h2>',
           f'<p class="secnote">{SIDIS_EFF_NOTE}</p>',
           '<div class="grid one">',
           figure_or_note_nu("cmp_faser_sidis_", "efficiency",
                             "The Tier E efficiency factor in z for charged "
                             "pions and kaons with the FASER selection, "
                             "flux-averaged, for &nu;<sub>e</sub> (left) and "
                             "&nu;<sub>&mu;</sub> (right). The line is the "
                             "mean of the five generators and the band their "
                             "envelope. The shaded strip is the z &lt; 0.1 "
                             "the yield study cuts. The factor already "
                             "contains the 60% hadron selection efficiency."),
           figure_or_note_nu("cmp_faser_sidis_", "efficiency_energy",
                             "The energy dependence of the same factors in "
                             "POWHEG-V2 + Pythia 8: the Tier E event "
                             "efficiency and &epsilon; for pions and kaons "
                             "at two values of z (top), and each divided by "
                             "the Tier E event efficiency (bottom). A single "
                             "generator is shown because GENIE's default "
                             "tune stops at 1 TeV, so a five-generator mean "
                             "would change composition along the axis."),
           '</div>']
    # the table a collaborator copies from
    for h, hlab in (("pi", "&pi;<sup>&plusmn;</sup>"),
                    ("K", "K<sup>&plusmn;</sup>")):
        out.append(f'<h3>{hlab}: &epsilon;(z), flux-averaged</h3>')
        out.append('<div class="tablewrap"><table><thead><tr>'
                   '<th>z bin</th>'
                   '<th>&epsilon;<sub>tier</sub> (&nu;<sub>&mu;</sub>)</th>'
                   '<th>&epsilon;<sub>acc</sub></th>'
                   '<th>&epsilon; (&nu;<sub>&mu;</sub>)</th>'
                   '<th>generator spread</th>'
                   '<th>&epsilon; (&nu;<sub>e</sub>)</th>'
                   '<th>generator spread</th></tr></thead><tbody>')
        for i in range(len(edges) - 1):
            row = [f'<td class="mono">{edges[i]:.2f}&ndash;{edges[i+1]:.2f}</td>']
            for fl in ("nu_mu", "nu_e"):
                # the mean of the five generators, and their envelope
                r = d["flux_averaged_mean"][fl][h]
                vals = ([r["eps_total_min"][i], r["eps_total_max"][i]]
                        if r["n_generators"][i] else [])
                if fl == "nu_mu":
                    row.append(f'<td class="mono">{r["eps_tier"][i]:.3f}</td>')
                    row.append(f'<td class="mono">{r["eps_acc"][i]:.3f}</td>')
                row.append(f'<td class="mono">{r["eps_total"][i]:.4f}</td>')
                row.append('<td class="mono">%s</td>'
                           % (f'{min(vals):.4f}&ndash;{max(vals):.4f}'
                              if vals else '&mdash;'))
            out.append('<tr>' + "".join(row) + '</tr>')
        out.append('</tbody></table></div>')
    out.append('<p class="secnote" style="font-size:0.85em">Script: '
               '<code>analysis/faser_sidis_efficiency.py</code> '
               '&middot; figure <code>analysis/plot_faser_sidis_efficiency.py</code> '
               '&middot; the per-energy factors, and the yields the ratio is '
               'built from, are in '
               '<code>results_nu/faser_sidis_efficiency.json</code>.</p>')
    return "".join(out)


def sidis_region_yields_section():
    """POWHEG-V2 in the benchmark region alone (user, 2026-09-21): the
    numbers an analytic NLO calculation with fragmentation functions is
    benchmarked against, before the Tier E factors of the previous view turn
    it into a FASERnu yield.  analysis/faser_sidis_region_yields.py; every
    number here is read from its JSON, so the prose cannot go stale."""
    p = f"{RESULTS_NU}/faser_sidis_region_yields.json"
    if not os.path.exists(p):
        return ('<h2>Region yields</h2><p class="missing">not run yet &mdash; '
                'analysis/faser_sidis_region_yields.py</p>')
    with open(p) as f:
        d = json.load(f)
    z = d["z_edges"]
    ff = d["flux_folded"]
    es = sorted(d["fixed_energy"], key=float)
    lo_e = float(es[0])
    cut = [i for i in range(len(z) - 1) if z[i] >= 0.1 - 1e-9]

    def tot(fl, h, idx=None):
        v = ff[fl][h]
        return sum(v[i] for i in (idx if idx is not None else range(len(v))))

    note = (
        "<b>The benchmark for an analytic calculation</b> (user, 2026-09-21). "
        "POWHEG-V2 + Pythia 8 with the benchmark settings (NNPDF4.0 NNLO, "
        "&mu;<sub>F</sub> = &mu;<sub>R</sub> = Q, massless charm and bottom, "
        "no QED) in the benchmark region alone: Q<sup>2</sup> &gt; 4 "
        "GeV<sup>2</sup> and W &gt; 3 GeV and <b>no other cut</b> &mdash; "
        "no Tier E, no cut on the hadron and no hadron selection efficiency. "
        "Every charged pion or kaon of the final state is counted (stable if "
        "c&tau; &gt; 10 mm), with z = E<sub>h</sub>/&nu; in the target rest "
        "frame, on tungsten per nucleon, (74&sigma;<sub>p</sub> + "
        "110&sigma;<sub>n</sub>)/184, with nuclear effects neglected. This "
        "is the denominator of the efficiency factors in the previous view, "
        "so the FASERnu yield is "
        "<code>dN<sub>FASER</sub>/dz = &epsilon;(z) &times; "
        "dN<sub>region</sub>/dz</code>.<br><br>"
        "<b>Two comparisons.</b> At fixed beam energy, d&sigma;/dz per nucleon "
        f"at the {len(es)} energies of the ladder "
        f"({', '.join(f'{float(e):g}' for e in es)} GeV): the like-for-like "
        "test of a fragmentation-function calculation, before any flux. "
        "Folded over the FASERnu flux, the expected hadron counts per z bin "
        f"at {d['lumi_fb']:g} fb<sup>&minus;1</sup> on 1.1 t of tungsten: "
        f"{ff['nu_mu']['events_region']:.0f} &nu;<sub>&mu;</sub> and "
        f"{ff['nu_e']['events_region']:.0f} &nu;<sub>e</sub> charged-current "
        f"events in the region, carrying {tot('nu_mu', 'pi', cut):.0f} and "
        f"{tot('nu_e', 'pi', cut):.0f} charged pions above z = 0.1. The fold "
        "evaluates the cross-section at the flux file's own energies, with "
        "the z spectrum interpolated in log E between the ladder points; "
        f"{ff['nu_mu']['events_below_lowest_sample'] / ff['nu_mu']['events_region']:.1%} "
        "of the &nu;<sub>&mu;</sub> region events "
        f"({ff['nu_e']['events_below_lowest_sample'] / ff['nu_e']['events_region']:.1%} "
        f"for &nu;<sub>e</sub>) sit below the lowest sampled energy, "
        f"{lo_e:g} GeV, where the rate is extrapolated and the z shape held "
        "at that point.<br><br>"
        "<b>The fluxes, the tables and a driver</b> are packaged in "
        "<code>share/faser-sidis/</code>: the FASERnu fluxes as CSV, the "
        "numbers of this view and the efficiency factors, and a standalone "
        "Python driver that folds any calculation tabulated in E and z over "
        "the flux, compares it with POWHEG-V2 and applies &epsilon;(z).")
    # every pion, not the hardest one (user, 2026-09-21: "how can the pion
    # yields be larger than the number of events?"); analysis/sidis_leading_vs_all.py
    lp = f"{RESULTS_NU}/sidis_leading_vs_all.json"
    if os.path.exists(lp):
        with open(lp) as f:
            lv = json.load(f)
        zz = lv["z_edges"]
        al, le = lv["tungsten"]["per_event_all"], lv["tungsten"]["per_event_leading"]
        ic = [i for i in range(len(zz) - 1) if zz[i] >= lv["z_tag"] - 1e-9]
        first = ic[0]
        agree = next(zz[i] for i in range(len(zz) - 1)
                     if all(abs(le[j] / al[j] - 1) < 0.005 for j in range(i, len(zz) - 1)
                            if al[j] > 0))
        m = lv["tungsten"]["n_above_ztag"]
        note += (
            "<br><br><b>Every charged pion is counted, not only the hardest one "
            "in the event</b>, which is why the pion counts exceed the number of "
            "events. This is the observable collinear factorisation describes, "
            "d&sigma;<sup>h</sup>/dz = &Sigma;<sub>q</sub> &sigma;<sub>q</sub> "
            "&otimes; D<sub>q</sub><sup>h</sup>(z), with the fragmentation "
            "function a number density, and the one HERMES and COMPASS publish as "
            "multiplicities. The hardest pion per event is a different "
            "observable, not given by a fragmentation function. On the 1 TeV "
            f"POWHEG-V2 sample a region event carries {sum(al[i] for i in ic):.2f} "
            f"charged pions above z = {lv['z_tag']:g}, and "
            f"{100 * m[0]:.0f}% / {100 * m[1]:.0f}% / {100 * m[2]:.0f}% / "
            f"{100 * sum(m[3:]):.0f}% of events have 0 / 1 / 2 / 3 or more of them. "
            f"Keeping only the hardest pion would leave {sum(le[i] for i in ic):.2f} per "
            f"event: {100 * le[first] / al[first]:.0f}% of the pions in the "
            f"{zz[first]:.2f}&ndash;{zz[first + 1]:.2f} bin, and all of them from "
            f"z = {agree:.2f} up, where no event in the sample has two.")
    out = ['<h2>Region yields: POWHEG-V2 with Q<sup>2</sup> &gt; 4 GeV<sup>2</sup> '
           'and W &gt; 3 GeV only</h2>',
           f'<p class="secnote">{note}</p>',
           '<div class="grid one">',
           figure_or_note_nu("cmp_faser_sidis_", "region_yields",
                             "POWHEG-V2 hadron counts per z bin in the region "
                             f"alone, folded over the FASERnu flux at {d['lumi_fb']:g} "
                             "fb<sup>&minus;1</sup>, for &nu;<sub>e</sub> "
                             "(left) and &nu;<sub>&mu;</sub> (right): charged "
                             "pions and kaons (solid), and the same times the "
                             "efficiency factor &epsilon;(z) of the previous "
                             "view, i.e. the FASERnu yield (dashed; the band is "
                             "the five-generator envelope of &epsilon;)."),
           figure_or_note_nu("cmp_faser_sidis_", "region_dsigma",
                             "POWHEG-V2 d&sigma;/dz per nucleon of tungsten in "
                             "the region alone, &nu;<sub>&mu;</sub> CC at each "
                             "energy of the ladder, for charged pions (left) "
                             "and kaons (right). &nu;<sub>e</sub> CC is the "
                             "same calculation."),
           '</div>']
    # the flux-folded counts
    out.append(f'<h3>Flux-folded counts per z bin, {d["lumi_fb"]:g} '
               'fb<sup>&minus;1</sup>, region only</h3>')
    out.append('<div class="tablewrap"><table><thead><tr><th>z bin</th>'
               '<th>&pi;<sup>&plusmn;</sup> (&nu;<sub>e</sub>)</th>'
               '<th>K<sup>&plusmn;</sup> (&nu;<sub>e</sub>)</th>'
               '<th>&pi;<sup>&plusmn;</sup> (&nu;<sub>&mu;</sub>)</th>'
               '<th>K<sup>&plusmn;</sup> (&nu;<sub>&mu;</sub>)</th>'
               '</tr></thead><tbody>')
    cells = [("nu_e", "pi"), ("nu_e", "K"), ("nu_mu", "pi"), ("nu_mu", "K")]
    for i in range(len(z) - 1):
        out.append(f'<tr><td class="mono">{z[i]:.2f}&ndash;{z[i+1]:.2f}</td>'
                   + "".join(f'<td class="mono">{ff[fl][h][i]:.1f}</td>'
                             for fl, h in cells) + '</tr>')
    for lab, idx in (("z &gt; 0.1", cut), ("all z", None)):
        out.append(f'<tr><td><b>{lab}</b></td>'
                   + "".join(f'<td class="mono"><b>{tot(fl, h, idx):.0f}</b></td>'
                             for fl, h in cells) + '</tr>')
    out.append('</tbody></table></div>')
    # the fixed-energy summary
    out.append('<h3>Fixed energy: the region cross-section and d&sigma;/dz '
               'per nucleon</h3>')
    out.append('<div class="tablewrap"><table><thead><tr><th>E<sub>&nu;</sub> '
               '[GeV]</th><th>&sigma;<sub>region</sub> [pb]</th>'
               '<th>&lang;n<sub>&pi;</sub>&rang;</th>'
               '<th>&lang;n<sub>K</sub>&rang;</th>'
               '<th>&sigma;<sub>&pi;</sub>, z &gt; 0.1 [pb]</th>'
               '<th>&sigma;<sub>K</sub>, z &gt; 0.1 [pb]</th>'
               '</tr></thead><tbody>')
    for e in es:
        r = d["fixed_energy"][e]
        out.append(f'<tr><td class="mono">{float(e):g}</td>'
                   f'<td class="mono">{r["sigma_region_pb"]:.4f}</td>'
                   f'<td class="mono">{r["pi"]["multiplicity"]:.3f}</td>'
                   f'<td class="mono">{r["K"]["multiplicity"]:.3f}</td>'
                   f'<td class="mono">{sum(r["pi"]["dsigma_pb"][i] for i in cut):.4f}</td>'
                   f'<td class="mono">{sum(r["K"]["dsigma_pb"][i] for i in cut):.4f}</td></tr>')
    out.append('</tbody></table></div>')
    out.append('<p class="secnote" style="font-size:0.85em">Every bin, species '
               '(&pi;<sup>+</sup>, &pi;<sup>&minus;</sup>, K<sup>+</sup>, '
               'K<sup>&minus;</sup> separately too) and energy is in '
               '<code>share/faser-sidis/data/powheg_v2/</code>. Script: '
               '<code>analysis/faser_sidis_region_yields.py</code> &middot; '
               'figures <code>analysis/plot_faser_sidis_region_yields.py</code> '
               '&middot; package <code>tools/export_sidis_share.py</code>.</p>')
    return "".join(out)


def sidis_panel():
    """SIDIS at FASER, a top-level tab (user, 2026-09-21).

    "In the report, please produce a separate tab 'SIDIS at FASER' and
    collect there all our SIDIS results, since I will use them to start a
    collaboration with some theory colleagues."  So the single-inclusive
    pion study and the NNLO comparison MOVE here out of "Predictions for
    FASER", the four paper figures are repeated here so the tab stands on
    its own, and the Tier E efficiency factors are new.
    """
    pions = f"{RESULTS_NU}/faser_pions.json"
    plots = {os.path.splitext(fn)[0]: (fn, mod, err)
             for fn, mod, err in _load_paper_plots()}
    figs = []
    # pp13_sidis_pions and pp14_pions_emulsion were dropped from the paper
    # plots (user, 2026-09-21: "pp15 and pp16 are sufficient"); the same
    # comparisons stay on this tab through the "gen" and "nnlo" sections.
    for slug in ("pp15_sidis_yields", "pp16_sidis_yields_kaon"):
        if slug in plots:
            figs.append(paper_plot_block(*plots[slug]))
    bodies = {
        "eff": sidis_efficiency_section(),
        "region": sidis_region_yields_section(),
        "gen": (pions_section_body(pions) if os.path.exists(pions) else
                '<p class="missing">analysis/faser_pions.py rates not run</p>'),
        "nnlo": sidis_section(),
        "figs": ("<h2>The SIDIS paper figures</h2>"
                 '<p class="secnote">The same figures as on the '
                 '&ldquo;Paper plots&rdquo; tab, repeated here so this tab '
                 'stands on its own.</p>' + "".join(figs)),
    }
    views = [("eff", "Tier E efficiency factors", ""),
             ("region", "Region yields (POWHEG-V2)", ""),
             ("gen", "Generator comparison", ""),
             ("nnlo", "Comparison with NNLO", ""),
             ("figs", "Paper figures", "")]
    return ('<h2>SIDIS at FASER</h2>'
            f'<p class="secnote">{SIDIS_TAB_NOTE}</p>'
            + tabs(views, "SIDIS view", ns="sidis-", cls="sel", bodies=bodies))


# FFs that are not in the comparison, and why (user, 2026-09-21: "a deep dive
# in the literature").  Facts as the papers and the LHAPDF index state them.
FF_NOT_INCLUDED = [
    ("HAPS-PiFF1.0, HAPS-KaFF1.0", "2606.16754",
     "SIA + HERMES and COMPASS SIDIS, NLO and NNLO (June 2026). Released in LHAPDF "
     "format according to the paper, but not yet in the LHAPDF index, and only "
     "charge-separated sets are described."),
    ("DSS07, DSS14 (&pi;), DSS17 (K)", "hep-ph/0703242, 1410.6027, 1702.06353",
     "Global NLO fits (SIA, SIDIS, pp). Distributed by the authors as their own "
     "grids and code, not in LHAPDF."),
    ("NPC23 NLO, main set", "2401.02781, 2407.04422",
     "In LHAPDF, but its grid starts at Q = 4 GeV, above the Q = 2 GeV of the "
     "region; its low-Q companion from 2502.17837 is used instead."),
    ("AKK08, HKNS07, KKP, Kretzer", "0803.2768, hep-ph/0702250, hep-ph/0010289, hep-ph/0003177",
     "Earlier fits, mostly to SIA; code or tables from the authors, not in LHAPDF."),
    ("NNFF1.1h, SHK22.h", "1807.03310, 2202.10779",
     "Unidentified charged hadrons, not pions or kaons."),
]


def ff_panel():
    """Fragmentation functions for the SIDIS study, a top-level tab (user,
    2026-09-21).  analysis/ff_comparison.py and plot_ff_comparison.py; every
    number in the prose is read from the JSON."""
    p = f"{RESULTS_NU}/ff_comparison.json"
    if not os.path.exists(p):
        return ('<h2>Fragmentation functions</h2><p class="missing">not run yet '
                '&mdash; analysis/ff_comparison.py</p>')
    with open(p) as f:
        d = json.load(f)
    import numpy as np
    z = np.array(d["z"])
    zlo, zhi = d["z_window_faser"]
    win = (z >= zlo - 1e-9) & (z <= zhi + 1e-9)
    w10 = d["nucc_weights"]["10"]

    def span(h, order, q="10"):
        sets = {r["family"]: r for r in d["sets"].values()
                if r["hadron"] == h and r["order"] == order}
        ref = np.array(sets[d["reference"][order]]["values"][q]["nucc"]["central"])[win]
        out = {}
        for fam, r in sets.items():
            v = r["values"][q]["nucc"]
            c = np.array(v["central"])[win]
            unc = (np.array(v["hi"]) - np.array(v["lo"]))[win] / 2 / c
            out[fam] = (float(np.min(c / ref)), float(np.max(c / ref)),
                        float(np.median(unc)))
        return out
    pin, kn = span("pi", "nlo"), span("K", "nlo")
    lo_pi = min(v[0] for v in pin.values())
    hi_pi = max(v[1] for v in pin.values())
    unc_pi = (min(v[2] for v in pin.values()), max(v[2] for v in pin.values()))
    rest_k = {f: v for f, v in kn.items() if f != "NNFF1.0"}
    note = (
        "<b>Which fragmentation functions a calculation of the FASER SIDIS "
        "yields could use</b> (user, 2026-09-21). FASER&nu; has no magnet, so "
        "it counts &pi;<sup>+</sup> + &pi;<sup>&minus;</sup> and "
        "K<sup>+</sup> + K<sup>&minus;</sup>: only the <b>charge-summed</b> "
        "fragmentation functions enter, D<sub>q</sub> = D<sub>q&#772;</sub>. "
        "Every family in LHAPDF with identified charged pions and kaons is "
        "compared, at each order it is released at, in the z range of the "
        f"study ({zlo:g} &lt; z &lt; {zhi:g}, unshaded) and at Q<sup>2</sup> = "
        f"{d['q2'][0]:g} and {d['q2'][1]:g} GeV<sup>2</sup>, from the "
        "bottom of the region to about its mean. Where a set stores only the "
        "positive hadron (JAM) the sum is built as D<sup>h+</sup><sub>q</sub> "
        "+ D<sup>h+</sup><sub>q&#772;</sub> replica by replica. Bands are 68% "
        "CL with each set's own prescription.<br><br>"
        "<b>What the measurement probes.</b> The first column combines the "
        "flavours with the leading-order weights of the quark that fragments "
        "in &nu;<sub>&mu;</sub> CC on tungsten at the region's typical "
        f"kinematics ({d['nucc_note'].split(', ')[2]}, {d['nucc_note'].split(', ')[3]}): "
        f"u {100 * w10['u']:.0f}%, c {100 * w10['c']:.0f}%, "
        f"d {100 * w10['d']:.1f}%, s {100 * w10['s']:.1f}%. The FASER yield is "
        "therefore a measurement of u &rarr; h, with a tenth of charm "
        "fragmentation; the gluon enters only beyond leading order.<br><br>"
        "<b>The sets disagree by more than they say.</b> On that combination "
        f"at Q<sup>2</sup> = {d['q2'][0]:g} GeV<sup>2</sup> the NLO pion sets "
        f"span {lo_pi:.2f} to {hi_pi:.2f} of {d['reference']['nlo']} across "
        f"{zlo:g} &lt; z &lt; {zhi:g}, while each quotes an uncertainty of "
        f"{100 * unc_pi[0]:.0f}&ndash;{100 * unc_pi[1]:.0f}%: the choice of set, "
        "not any one band, is the fragmentation uncertainty of a FASER "
        "prediction. For kaons NNFF1.0, a fit to e<sup>+</sup>e<sup>&minus;</sup> "
        f"data alone, is {kn['NNFF1.0'][0]:.1f}&ndash;{kn['NNFF1.0'][1]:.1f} times "
        f"{d['reference']['nlo']} on this u-dominated combination, a flavour "
        "separation that data cannot make; the remaining sets span "
        f"{min(v[0] for v in rest_k.values()):.2f} to "
        f"{max(v[1] for v in rest_k.values()):.2f} of it. For an NLO calculation "
        "the NLO sets are the consistent choice; MAPFF1.0, NPC23 and JAM24 are "
        "the most recent fits that include SIDIS.")
    figs = {}
    for h, hl in (("pi", "&pi;<sup>+</sup> + &pi;<sup>&minus;</sup>"),
                  ("K", "K<sup>+</sup> + K<sup>&minus;</sup>")):
        body = ['<div class="grid one">']
        for order, ol in (("nlo", "NLO"), ("nnlo", "NNLO")):
            body.append(figure_or_note_nu(
                "cmp_ff_", f"{h}_{order}",
                f"Charge-summed {hl} fragmentation functions z D(z) at {ol}: "
                "the leading-order &nu;<sub>&mu;</sub> CC combination and the "
                "flavours u, d, s, c, g (columns), at Q<sup>2</sup> = "
                f"{d['q2'][0]:g} and {d['q2'][1]:g} GeV<sup>2</sup> (row "
                f"pairs), with the ratio to {d['reference'][order]} of the "
                f"same order. The shaded strips are outside {zlo:g} &lt; z "
                f"&lt; {zhi:g}."))
        body.append('</div>')
        figs[h] = "".join(body)
    rows = []
    for fam, m in d["families"].items():
        names = [n for n, r in d["sets"].items() if r["family"] == fam]
        qmin = min(d["sets"][n]["qmin"] for n in names)
        rows.append(f'<tr><td><b>{fam}</b></td><td>{m["who"]}</td>'
                    f'<td class="mono">{m["arxiv"]}</td><td>{m["data"]}</td>'
                    f'<td>{m["orders"]}</td><td>{m["errors"]}; {m["method"]}</td>'
                    f'<td class="mono">{qmin:.2f}</td>'
                    f'<td class="mono">{"<br>".join(names)}</td></tr>')
    sets_html = (
        '<h3>The sets compared</h3><div class="tablewrap"><table><thead><tr>'
        '<th>family</th><th>authors</th><th>arXiv</th><th>data</th>'
        '<th>orders released</th><th>uncertainties</th>'
        '<th>Q<sub>min</sub> [GeV]</th><th>LHAPDF sets used</th>'
        '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
        '<h3>Not included, and why</h3><div class="tablewrap"><table><thead>'
        '<tr><th>set</th><th>arXiv</th><th>reason</th></tr></thead><tbody>'
        + "".join(f'<tr><td><b>{a}</b></td><td class="mono">{b}</td><td>{c}</td></tr>'
                  for a, b, c in FF_NOT_INCLUDED)
        + '</tbody></table></div>'
        '<p class="secnote" style="font-size:0.85em">Scripts: '
        '<code>analysis/ff_comparison.py</code> &middot; '
        '<code>analysis/plot_ff_comparison.py</code> &middot; data '
        '<code>results_nu/ff_comparison.json</code>. The LHAPDF sets are '
        'installed from lhapdfsets.web.cern.ch.</p>')
    views = [("pi", "Pions", ""), ("K", "Kaons", ""), ("sets", "The sets", "")]
    return ('<h2>Fragmentation functions for the SIDIS study</h2>'
            f'<p class="secnote">{note}</p>'
            + tabs(views, "FF view", ns="ff-", cls="sel",
                   bodies={"pi": figs["pi"], "K": figs["K"], "sets": sets_html}))


def faser_rates_panel():
    """Applications of the benchmark to FASER, as a top-level tab.

    A TOP-LEVEL TAB (user, 2026-09-01), one application per sub-tab (user,
    2026-09-07).  Since 2026-09-19 only the applications built on the final
    region or on region-free inputs remain: the dimuon cut flow and the
    single-inclusive pion study (NEUTRINO ONLY, the rule-2b carve-out the
    user granted on 2026-09-07).  The inclusive CC and NC rates, the charm
    application and the Run 3/Run 4 projections were computed in the
    0.2 < y < 0.9 region and were removed with it (user, 2026-09-19).
    """
    # THE SIDIS STUDY MOVED OUT on 2026-09-21 (user: "produce a separate tab
    # 'SIDIS at FASER' and collect there all our SIDIS results"), so this tab
    # holds the dimuon application alone.
    bodies = {"dimuon": dimuon_cutflow_section()}
    views = [("dimuon", "Dimuon production in CC DIS", "")]
    return ('<h2>Predictions for FASER</h2>'
            f'<p class="secnote">{FASER_RATES_NOTE}</p>'
            + tabs(views, "FASER application", ns="faser-", cls="sel",
                   bodies=bodies))


FASER_PIONS_NOTE = (
    "<b>Single-inclusive charged pions in the emulsion</b> (user, 2026-09-07): "
    "the feasibility of a FASER measurement of &nu;N &rarr; &mu;&pi;X, the "
    "process the NNLO semi-inclusive DIS calculations of "
    "<a href=\"https://arxiv.org/abs/2510.00100\">2510.00100</a> (neutral and "
    "charged current) describe and fragmentation functions are fitted to. For "
    "every DIS event of each sample, every charged pion and kaon in the final "
    "state &mdash; stable under the ctau &gt; 10 mm convention, so the primary "
    "hadrons the emulsion sees at the vertex &mdash; is recorded with its lab "
    "energy, its fraction z = E<sub>h</sub>/&nu; of the hadronic energy and its "
    "transverse momentum, inside the tan&thinsp;&theta; &lt; 0.5 track acceptance "
    "of the emulsion.<br><br>"
    "<b>In the benchmark region, on the tungsten samples (user, 2026-09-18).</b> The "
    "region is the benchmark region, Q<sup>2</sup> &gt; 4 GeV<sup>2</sup> "
    "and W &gt; 3 GeV with <b>no y window</b>, and the FASER Tier E selection "
    "nested in it; the samples are those of the paper plots, tungsten by "
    "isospin with the proton and the neutron generated separately and "
    "combined per region with the weights 74 &sigma;<sub>p</sub> : 110 "
    "&sigma;<sub>n</sub>. Every generator is here, Herwig and Sherpa included. "
    "The spectra are per event; the rate is the "
    "region&rsquo;s tungsten cross-section on the SIDIS ladder of 300, 400, "
    "700, 1000, 2000 and 4000 GeV, carried across the flux on the YADISM NLO "
    "(ZM-VFNS, with target-mass corrections) shape of the region. The 300 GeV point exists because the "
    "tier asks for a scattered lepton above 200 GeV: its efficiency is 0.48 at "
    "1 TeV, 0.14 at 400 GeV, 0.03 at 300 GeV and exactly zero at and below "
    "200 GeV, while a third of the flux-weighted rate of the region sits below "
    "400 GeV &mdash; holding the 400 GeV efficiency constant below it, as the "
    "earlier version of this study did, overestimated the yield by several "
    "per cent.<br><br>"
    "<b>This study shows the neutrino current only</b> (user, 2026-09-07): "
    "\u201cto simplify the discussion &hellip; show only neutrino DIS here, "
    "remove the muon DIS in this specific study\u201d. It is a deliberate exception to "
    "carrying the two currents together; the muon-side spectra are still "
    "extracted and still on disk.<br><br>"
    "<b>The numbers.</b> At 300 fb<sup>&minus;1</sup> POWHEG-V2 expects "
    "2313 Tier E events with 6.88 charged pions each inside the emulsion "
    "acceptance; at a 60% hadron selection efficiency that is 9.5 thousand "
    "selected pions. Herwig gives 7.37 per event, "
    "Sherpa 7.11, GENIE&rsquo;s default tune 5.80 and HEDIS 5.85, so the "
    "multiplicity spreads by a factor 1.27 across the matched generators and "
    "FASER&rsquo;s own, and GENIE delivers 6.7 thousand selected pions. The shape "
    "differs more than the normalisation, and it is the hadronisation model: "
    "below z = 0.1 GENIE has 0.67 of POWHEG-V2&rsquo;s pions and above z = 0.5 "
    "it has 1.47 times as many &mdash; a harder fragmentation spectrum from "
    "AGKY than from the Lund string &mdash; while Sherpa&rsquo;s cluster model "
    "goes the other way, 1.16 and 0.72, and Herwig&rsquo;s cluster model sits "
    "between, 1.10 and 0.94. Kaons are 12% of the pions for POWHEG-V2 and 15% "
    "for GENIE, 1.0 to 1.4 thousand selected. A measurement binned in z with these "
    "statistics resolves the models at every z, and it is a quantity "
    "on which the generators disagree at the level of a factor "
    "rather than of per cent. No analytic reference exists for a hadron "
    "spectrum, so this comparison is generator against generator; the section "
    "below sets the same spectra against a calculation. The centre panel's "
    "GENIE default-tune curve ends at E<sub>h</sub> = 1 TeV because that tune "
    "is declared valid to 1 TeV and its spectrum is held at the 1 TeV one "
    "above it."
)


def pions_section_body(p):
    with open(p) as f:
        d = json.load(f)
    out = ['<h2>Single-inclusive pion production in the emulsion</h2>',
           f'<p class="secnote">{FASER_PIONS_NOTE}</p>']
    # NEUTRINO ONLY (user, 2026-09-07), the one carve-out from rule 2b on this
    # page; the muon-side results are still computed and still on disk.
    for cur, beam in (("nu", "Neutrino charged current"),):
        out.append(f'<h2>{beam}, FASER Tier E in the benchmark region, '
                   f'300 fb<sup>&minus;1</sup></h2>')
        out.append('<div class="grid one">' + figure_or_note_nu(
            f"cmp_faser_pions_{cur}_pi", "",
            "Charged pions inside the emulsion track acceptance, per bin of "
            "z = E<sub>h</sub>/&nu; (left), lab energy (centre) and transverse "
            "momentum (right), in events passing FASER Tier E at "
            "300 fb<sup>&minus;1</sup>; lower panels the ratio to POWHEG-V2. "
            "The vertical axis is a count in the bin, not a density.")
                   + '</div>')
        out.append('<div class="grid one">' + figure_or_note_nu(
            f"cmp_faser_pions_{cur}_K", "",
            "The same for charged kaons.") + '</div>')
        gens = d["currents"][cur]["generators"]
        out.append('<div class="tablewrap"><table><thead><tr><th>generator</th>'
                   '<th>Tier E events</th><th>&pi;<sup>&plusmn;</sup> per event, all</th>'
                   '<th>&pi;<sup>&plusmn;</sup> per event, emulsion acceptance</th>'
                   '<th>&pi;<sup>&plusmn;</sup> in the emulsion</th>'
                   '<th>K<sup>&plusmn;</sup> in the emulsion</th>'
                   '<th>&pi;<sup>&plusmn;</sup> per event, benchmark region (no tier)</th>'
                   '<th>rate normalisation</th></tr></thead><tbody>')
        for key, g in gens.items():
            te = g["selections"].get("sidis_e")
            inc = g["selections"].get("sidis_w")
            if not te:
                continue
            out.append(f'<tr><td>{g["label"]}</td><td class="num">{te["events"]:.0f}</td>'
                       f'<td class="num">{te["pi_all_per_event"]:.2f}</td>'
                       f'<td class="num">{te["pi_emul_per_event"]:.2f}</td>'
                       f'<td class="num">{te["pi_emul_total"]:.0f}</td>'
                       f'<td class="num">{te["K_emul_total"]:.0f}</td>'
                       f'<td class="num">{inc["pi_emul_per_event"] if inc else float("nan"):.2f}</td>'
                       f'<td>{te["normalisation"]}</td></tr>')
        out.append('</tbody></table></div>')
    out.append('<p class="secnote">Both columns of pions per event count pions '
               'in events passing Tier E; the last numeric column is the same '
               'in the benchmark region without the tier, so the effect of the '
               'multiplicity requirement on the multiplicity itself is '
               'visible. The per-event numbers are from each sample&rsquo;s '
               'own events, tungsten by isospin; the totals multiply them by '
               'the Tier E rate on the energy ladder (the YADISM NLO (ZM-VFNS) shape of '
               'the benchmark region times each generator&rsquo;s own ratio, with the tier '
               'efficiency measured at 300, 400, 700, 1000, 2000 and 4000 GeV '
               'and exactly zero at 200 GeV).</p>')
    return "".join(out)


# The public note.  The full one adds the comparison with the NNLO calculation
# of arXiv:2504.05376, whose curves were provided by its authors and are not
# distributed; it is read from data/faser_internal/ where it exists.
FASER_SIDIS_NOTE = (
    '<b>Charged pions in z at FASER&nu;</b>: the yields, their '
    'statistical reach, and how often a DIS event carries a pion, in the '
    'benchmark region on tungsten.<br><br><b>The region.</b> The '
    'benchmark region, Q<sup>2</sup> &gt; 4 GeV<sup>2</sup> and W &gt; 3 '
    'GeV, with <b>no y window</b> and <b>no cut on x</b> on any yield: '
    'the paper cuts x &gt; 0.1 and we account for it only where we '
    'compare with their numbers; at these energies it would keep 0.72 of '
    'the Tier E rate and is used in no yield here. The yields come from '
    'the samples on tungsten by isospin, the comparison with the '
    'calculation from the proton sample at 300 GeV, since the calculation '
    'is on a proton.<br><br><b>The yields.</b> Every yield is neutrino '
    'plus antineutrino charged-current scattering (since 2026-10-04), '
    'each generator with its own antineutrino samples and the '
    'antineutrino flux of the same flavour. At 300 fb<sup>&minus;1</sup> '
    'POWHEG-V2 expects 2813 &nu;<sub>&mu;</sub> + '
    '&nu;&#772;<sub>&mu;</sub> charged-current Tier E events, 18% of them '
    'from the antineutrino, carrying, at a 60% hadron selection '
    'efficiency, 1513 &pi;<sup>+</sup> and 1167 &pi;<sup>&minus;</sup> '
    'above z = 0.1 inside the emulsion track acceptance &mdash; 2680 '
    'pions, with a statistical error of 2.1%. The &nu;<sub>e</sub> flux '
    'is eight times smaller but considerably harder &mdash; a mean energy '
    'of 733 GeV against 312 GeV, since it is fed by charm decay rather '
    'than by pions and kaons &mdash; so Tier E keeps 0.42 of it against '
    '0.31 (neutrinos), and the same cross-section gives 1154 events and '
    '1074 pions. Across the five generators the &nu;<sub>&mu;</sub> + '
    '&nu;&#772;<sub>&mu;</sub> yields span 2292 to 3156 events and 2091 '
    'to 2925 pions. The charge ratio '
    '&pi;<sup>+</sup>/&pi;<sup>&minus;</sup> is 1.30 for POWHEG-V2, '
    'between 1.25 and 1.31 for every generator, and it falls with z in '
    'all of them, because &nu;p &rarr; '
    '&mu;<sup>&minus;</sup>&pi;<sup>+</sup>X runs on the valence d &rarr; '
    'u transition already at Born level while a &pi;<sup>&minus;</sup> '
    'needs a sea combination; the antineutrino, which favours '
    '&pi;<sup>&minus;</sup>, dilutes the asymmetry. A charge-summed '
    'spectrum averages that away, which is why the yields are '
    'split.<br><br><b>How far in z a measurement reaches.</b> Pions from '
    'one event are not independent counts, so the statistical error is '
    'the compound-Poisson one and exceeds the square root of the pion '
    'count. The selection efficiency works the other way, since keeping '
    'each pion independently decorrelates the pions within an event: the '
    'inflation is 1.15 rather than the 1.24 of the generated sample, so '
    'the error grows by less than 1/&radic;0.6. Per bin of 0.05 in z it '
    'stays below 10% up to z = 0.40 on the &nu;<sub>&mu;</sub> charged '
    'current and z = 0.25 for &nu;<sub>e</sub> &mdash; which is where the '
    'generators and the calculation are already tens of per cent apart. '
    'Sherpa&rsquo;s MC@NLO weights turn its highest z bin negative in a '
    'few of the yield spectra; those bins are clipped at zero and the '
    'JSON records which.<br><br><b>How often a DIS event has a pion at '
    'all</b> (user, 2026-09-07). In the region with no other cut, 98.6% '
    'of &nu;<sub>&mu;</sub> charged-current events contain at least one '
    'charged pion by POWHEG-V2 and 98.9% by GENIE&rsquo;s default tune; '
    'requiring one above z = 0.1 inside the emulsion acceptance takes '
    'that to 78.5% and 75.3%.'
)
SIDIS_NNLO_PRIVATE = False
_SIDIS_NOTE_PRIVATE = os.path.join(BASE, "data", "faser_internal", "sidis_nnlo_report_note.html")
if os.path.exists(_SIDIS_NOTE_PRIVATE):
    with open(_SIDIS_NOTE_PRIVATE) as _f:
        FASER_SIDIS_NOTE = _f.read()
    SIDIS_NNLO_PRIVATE = True


# ---------------------------------------------------------------- SIDIS vs NNLO
# NEUTRINO ONLY (user, 2026-09-07): "to simplify the discussion (I know is a
# rule violation, but accepted) show only neutrino DIS here, remove the muon
# DIS in this specific study".  The muon-side numbers are still computed and
# still in results_nu/faser_sidis.json; they are not displayed.
SIDIS_FLAVOURS = [("nu", "nu_mu", "&nu;<sub>&mu;</sub> CC"),
                  ("nu", "nu_e", "&nu;<sub>e</sub> CC")]


def _sidis_excluded(d):
    """Which generators cannot appear here, and why. Never silent."""
    seen = {}
    for cur, fname, _lab in SIDIS_FLAVOURS:
        for key, why in d["currents"][cur]["flavours"][fname]["excluded"].items():
            seen.setdefault(why, key)
    if not seen:
        return ""
    items = "".join(f"<li>{why}</li>" for why in seen)
    return ('<p class="secnote"><b>Not shown here, and why.</b><ul>'
            + items + '</ul></p>')


def _sidis_cutflow_table(d):
    """The nested chain of cuts, and what each one costs."""
    rows = [k for k, _ in [(r["key"], r["label"]) for r in d["regions"]]]
    labels = {r["key"]: r["label"] for r in d["regions"]}
    out = ['<div class="tablewrap"><table><thead><tr><th>beam</th>'
           '<th>generator</th>'
           + "".join(f'<th>{labels[k]}</th>' for k in rows)
           + '</tr></thead><tbody>']
    for cur, fname, flabel in SIDIS_FLAVOURS:
        for key, g in d["currents"][cur]["flavours"][fname]["generators"].items():
            r = g["regions"]
            cells = []
            for k in rows:
                cells.append(f'<td class="num">{r[k]["events"]:.0f}</td>'
                             if k in r else '<td class="num">&mdash;</td>')
            out.append(f'<tr><td>{flabel}</td><td>{g["label"]}</td>'
                       + "".join(cells) + '</tr>')
    out.append('</tbody></table></div>')
    return "".join(out)


def _sidis_yield_table(d):
    """Charged pion yields with z > 0.1, per beam and generator."""
    reg = d["yield_region"]
    zc = d["z_min"]
    out = [f'<div class="tablewrap"><table><thead><tr><th>beam</th>'
           f'<th>generator</th><th>events</th>'
           f'<th>&pi;<sup>+</sup>, z &gt; {zc:g}</th>'
           f'<th>&pi;<sup>&minus;</sup>, z &gt; {zc:g}</th>'
           f'<th>&pi;<sup>&plusmn;</sup> total</th>'
           f'<th>statistical error</th>'
           f'<th>&pi;<sup>&minus;</sup>/&pi;<sup>+</sup></th>'
           f'<th>K<sup>&plusmn;</sup></th>'
           f'<th>rate normalisation</th></tr></thead><tbody>']
    for cur, fname, flabel in SIDIS_FLAVOURS:
        for key, g in d["currents"][cur]["flavours"][fname]["generators"].items():
            r = g["regions"].get(reg)
            if not r:
                continue
            pp, pm = r["pip_total_zcut"], r["pim_total_zcut"]
            tot, err = r["pi_total_zcut"], r["pi_total_zcut_err"]
            out.append(f'<tr><td>{flabel}</td><td>{g["label"]}</td>'
                       f'<td class="num">{r["events"]:.0f}</td>'
                       f'<td class="num">{pp:.0f}</td>'
                       f'<td class="num">{pm:.0f}</td>'
                       f'<td class="num">{tot:.0f}</td>'
                       f'<td class="num">&plusmn;{err:.0f} '
                       f'({100 * err / tot if tot else float("nan"):.1f}%)</td>'
                       f'<td class="num">{(pm / pp if pp else float("nan")):.3f}</td>'
                       f'<td class="num">{r["K_total_zcut"]:.0f}</td>'
                       f'<td>{r["normalisation"]}</td></tr>')
    out.append('</tbody></table></div>')
    return "".join(out)


def _sidis_zbin_table(d):
    """The reference generator's yields bin by bin in z, with their errors."""
    reg = d["yield_region"]
    edges = d["z_edges"]
    zc = d["z_min"]
    cols = []
    for cur, fname, flabel in SIDIS_FLAVOURS:
        key = d["reference"][cur]
        g = d["currents"][cur]["flavours"][fname]["generators"].get(key)
        if not g or reg not in g["regions"]:
            continue
        cols.append((flabel, g["label"], g["regions"][reg]))
    if not cols:
        return ""
    head = "".join(f'<th>{lab} &pi;<sup>+</sup></th>'
                   f'<th>{lab} &pi;<sup>&minus;</sup></th>'
                   f'<th>{lab} &pi;<sup>&plusmn;</sup> &plusmn; stat</th>'
                   for lab, _, _ in cols)
    out = [f'<div class="tablewrap"><table><thead><tr><th>z bin</th>{head}</tr>'
           '</thead><tbody>']
    tot = [0.0] * (3 * len(cols))
    var = [0.0] * len(cols)
    for i in range(len(edges) - 1):
        if edges[i] < zc - 1e-12:
            continue
        cells = []
        for j, (_lab, _gl, r) in enumerate(cols):
            a, b = r["pip_z"][i], r["pim_z"][i]
            t, e = r["pi_z"][i], r["pi_z_err"][i]
            tot[3 * j] += a
            tot[3 * j + 1] += b
            tot[3 * j + 2] += t
            var[j] += e * e
            cells.append(f'<td class="num">{a:.0f}</td>'
                         f'<td class="num">{b:.0f}</td>'
                         f'<td class="num">{t:.0f} &plusmn; {e:.0f}</td>')
        out.append(f'<tr><td>{edges[i]:.2f}&ndash;{edges[i + 1]:.2f}</td>'
                   + "".join(cells) + '</tr>')
    cells = []
    for j in range(len(cols)):
        cells.append(f'<td class="num"><b>{tot[3 * j]:.0f}</b></td>'
                     f'<td class="num"><b>{tot[3 * j + 1]:.0f}</b></td>'
                     f'<td class="num"><b>{tot[3 * j + 2]:.0f} &plusmn; '
                     f'{var[j] ** 0.5:.0f}</b></td>')
    out.append('<tr><td><b>total</b></td>' + "".join(cells)
               + '</tr></tbody></table></div>')
    return "".join(out)


def _sidis_presence_table(d):
    """What fraction of DIS events has a charged pion in it."""
    zc = d["z_min"]
    out = ['<div class="tablewrap"><table><thead><tr><th>beam</th>'
           '<th>generator</th><th>region</th>'
           '<th>any charged pion</th>'
           f'<th>one with z &gt; {zc:g}</th>'
           f'<th>one with z &gt; {zc:g} in the emulsion acceptance</th>'
           '</tr></thead><tbody>']
    labels = {r["key"]: r["label"] for r in d["regions"]}
    for cur, fname, flabel in SIDIS_FLAVOURS:
        if fname == "nu_e":
            continue          # identical to nu_mu: same events, other flux
        for key, g in d["currents"][cur]["flavours"][fname]["generators"].items():
            for region in ("sidis_w", d["yield_region"]):
                r = g["regions"].get(region)
                if not r:
                    continue
                p = r["pion_presence"]
                out.append(f'<tr><td>{flabel}</td><td>{g["label"]}</td>'
                           f'<td>{labels.get(region, region)}</td>'
                           f'<td class="num">{100 * p["any_z"]:.2f}%</td>'
                           f'<td class="num">{100 * p["zcut"]:.2f}%</td>'
                           f'<td class="num">{100 * p["zcut_emulsion"]:.2f}%</td>'
                           '</tr>')
    out.append('</tbody></table></div>')
    return "".join(out)


def _sidis_mult_table(d):
    """Our multiplicity against the published NNLO curve, at their energy."""
    c = d.get("comparison")
    if not c:
        return ""
    zc = d["z_min"]
    out = [f'<div class="tablewrap"><table><thead><tr><th>calculation or '
           f'generator</th><th>beam energy</th>'
           f'<th>&pi;<sup>+</sup> per DIS event, z &gt; {zc:g}</th>'
           f'<th>ratio to the NNLO</th>'
           f'<th>mean dM/dz over the NNLO, z &gt; 0.2</th>'
           f'<th>the same without their x &gt; 0.1</th>'
           f'<th>selected events in the sample</th></tr></thead><tbody>']
    out.append(f'<tr><td>Bonino et al. NNLO</td>'
               f'<td class="num">{c["energy_gev"]:.0f} GeV</td>'
               f'<td class="num">{c["reference_total_zcut"]:.3f}</td>'
               f'<td class="num">1</td><td class="num">1</td>'
               f'<td class="num">&mdash;</td><td class="num">&mdash;</td></tr>')
    for key, g in c["generators"].items():
        flag = "" if g["matched_energy"] else " (not matched)"
        out.append(f'<tr><td>{g["label"]}</td>'
                   f'<td class="num">{g["energy_used_gev"]:.0f} GeV{flag}</td>'
                   f'<td class="num">{g["total_zcut"]:.3f}</td>'
                   f'<td class="num">{g["ratio_total_zcut"]:.3f}</td>'
                   f'<td class="num">{g["mean_ratio_z_gt_0p2"]:.2f}</td>'
                   f'<td class="num">{g["total_zcut_no_x"]:.3f}</td>'
                   f'<td class="num">{g["n_selected"]:,}</td></tr>')
    out.append('</tbody></table></div>')
    return "".join(out)


def sidis_section():
    """Charged-pion yields in z against arXiv:2504.05376 (2026-09-07)."""
    p = f"{RESULTS_NU}/faser_sidis.json"
    if not os.path.exists(p):
        return ('<h2>Comparison with the NNLO SIDIS calculation</h2>'
                '<p class="missing">not run yet &mdash; analysis/faser_sidis.py</p>')
    with open(p) as f:
        d = json.load(f)
    nnlo = bool(d.get("comparison")) and SIDIS_NNLO_PRIVATE
    zc = d["z_min"]
    e_lab = (d.get("comparison") or {}).get("energy_gev", d["paper_energy_gev"])
    out = ['<h2>Charged-pion yields in z, against the NNLO SIDIS '
           'calculation</h2>' if nnlo else '<h2>Charged-pion yields in z</h2>',
           f'<p class="secnote">{FASER_SIDIS_NOTE}</p>']
    if nnlo:
        out += [
           '<div class="grid one">' + figure_or_note_nu(
               "cmp_faser_sidis_mult", "",
               "The pion multiplicity dM/dz per DIS event in the region of "
               "arXiv:2504.05376, x &gt; 0.1 and W &gt; 3 GeV, at their beam "
               f"energy of {e_lab:.0f} GeV. Left: our generators under their "
               "published LO, NLO and NNLO curves with seven-point scale "
               "bands; the box gives the integral above "
               f"z = {zc:g} for each. Right: the ratio to their NNLO, with "
               "its own band about one.") + '</div>',
           _sidis_mult_table(d),
           '<p class="secnote">The last numeric column drops the x &gt; 0.1 '
           'of the calculation and keeps everything else, so what their x '
           'floor does to a multiplicity is a measured number rather than an '
           'assumption. The comparison is made at their beam energy: 300 GeV '
           'is not one of the benchmark energies, so a sample was generated '
           'for it.</p>']
    else:
        out.append('<p class="missing">The comparison with the NNLO calculation '
                   'of arXiv:2504.05376 uses curves provided by its authors, '
                   'which are not distributed with this repository.</p>')
    out += [
           f'<h2>Event yields in z at {d["lumi_fb"]:.0f} '
           'fb<sup>&minus;1</sup></h2>',
           '<div class="grid one">' + figure_or_note_nu(
               "cmp_faser_sidis_yields", "",
               "Charged pions inside the emulsion track acceptance, "
               "counted per bin of z, in events passing FASER Tier E with "
               f"W &gt; 3 GeV at {d['lumi_fb']:.0f} fb<sup>&minus;1</sup>. "
               "The vertical axis is a count in the bin, not a density. "
               "There is no y "
               "window and no cut on x. Solid pi<sup>+</sup>, dotted "
               "pi<sup>&minus;</sup>; the shaded band at low z is the "
               f"z &lt; {zc:g} that is cut. Lower panels: the "
               "pi<sup>&minus;</sup>/pi<sup>+</sup> ratio, which is where the "
               "flavour sensitivity of neutrino SIDIS sits.") + '</div>',
           _sidis_yield_table(d),
           '<p class="secnote">Yields are for the tungsten target of this '
           'tab, 74 protons and 110 neutrons per nucleus. The two neutrino '
           'rows of each generator share one cross-section and differ only '
           'through the flux: charged-current DIS is the same calculation for '
           'an electron and for a muon in the final state. What is not '
           'modelled is the detector side of that difference &mdash; an '
           'electron showers in the emulsion where a muon leaves a track, and '
           'Tier E&rsquo;s cuts on the scattered lepton are applied to it as '
           'if it were reconstructed as cleanly.</p>',
           '<h2>How big the statistical error would be</h2>',
           '<div class="grid one">' + figure_or_note_nu(
               "cmp_faser_sidis_stat", "",
               "Left: the expected yield per bin of z for POWHEG-V2, the "
               "reference generator, with the statistical error a "
               "measurement would have. Right: that error as a fraction of "
               "the yield, with the naive square root of the pion count "
               "dotted beneath it. Pions from one event are not independent "
               "counts, so the error is the compound-Poisson one, "
               "sqrt(N<sub>ev</sub> &lang;n<sup>2</sup>&rang;), and it is "
               "larger.") + '</div>',
           '<h2>What fraction of DIS events has a pion in it</h2>',
           _sidis_presence_table(d),
           '<p class="secnote">1 &minus; P(no charged pion), per DIS event in '
           'the region named. The last column is the one a measurement sees: '
           'a charged pion above the z cut and inside the emulsion track '
           'acceptance. The two neutrino flavours share these numbers, since '
           'they share the events.</p>',
           '<h2>The chain of cuts</h2>',
           _sidis_cutflow_table(d),
           '<p class="secnote">Expected events at '
           f'{d["lumi_fb"]:.0f} fb<sup>&minus;1</sup>, each column the one to '
           'its left plus one cut, starting from the benchmark region '
           '(Q&sup2; &gt; 4 GeV&sup2;, W &gt; 3 GeV, no y window, tungsten by '
           'isospin). The x &gt; 0.1 column is a diagnostic: no '
           'yield quoted on this page uses it.</p>',
           '<h2>Yields bin by bin in z</h2>',
           _sidis_zbin_table(d),
           '<p class="secnote">POWHEG-V2, the reference generator, in events '
           'passing FASER Tier E with W &gt; 3 GeV, counting pions with '
           f'z &gt; {zc:g} inside the emulsion track acceptance.</p>',
           _sidis_excluded(d)]
    return "".join(out)


FASER_ELEC_NOTE = (
    "<b>FASER measures muon neutrinos twice over with its electronic "
    "detector</b> &mdash; the muon is reconstructed in the spectrometer while "
    "the FASER&nu; tungsten is only the target &mdash; and both measurements "
    "are here (user, 2026-09-08): "
    "<a href=\"https://repository.cern/records/eg0mj-1ep53\">"
    "CERN-FASER-CONF-2026-005</a> at 186 fb<sup>&minus;1</sup>, with 766.8 "
    "observed events unfolding to 3052.3 interactions in the fiducial "
    "volume, and <a href=\"https://arxiv.org/abs/2412.03186\">"
    "arXiv:2412.03186</a> at 65.6 fb<sup>&minus;1</sup> with 338.1.<br><br>"
    "<b>The observable is &minus;L/E<sub>&nu;</sub>, and that is why no "
    "rapidity information is needed.</b> L is the lepton number, so the "
    "<b>sign</b> of the axis is the neutrino charge and its magnitude is "
    "1/E: one variable carrying energy and charge together, which our "
    "one-dimensional per-charge flux files predict directly. The bin above "
    "1 TeV holds both charges, because the spectrometer cannot sign a track "
    "that stiff. The note&rsquo;s other figures <i>are</i> differential in "
    "rapidity; the flux files carry none, so those are not reproduced "
    "here.<br><br>"
    "<b>The left panel is a cross-section against a cross-section.</b> "
    "FASER publishes the flux-weighted mean charged-current cross-section "
    "per nucleon in each bin (Table III of the paper), so no detector, no "
    "geometry and no luminosity enter either side. Our GENIE in FASER&rsquo;s "
    "own tune lands within 4 to 9% of it in every bin &mdash; a check of our "
    "machinery rather than of the physics, since their simulation is GENIE "
    "too &mdash; and POWHEG-V2 sits between 0.83 and 1.07 of it, low in the "
    "two 100&ndash;300 GeV bins because our ladder is generated above "
    "Q&sup2; &gt; 4 GeV&sup2; and 4 to 8% of the rate lies below that "
    "floor.<br><br>"
    "<b>The right panel needs the geometry, and one approximation with "
    "it.</b> The measured fiducial volume is a cylinder of 100 mm radius; "
    "the flux files are a <i>count</i> through their authors&rsquo; "
    "25 &times; 25 cm aperture. Converting between them assumes the flux is "
    "uniform across that aperture, and it is not &mdash; it peaks on axis. "
    "Rather than hide the assumption, the size of the miss is measured: our "
    "folded flux is 0.79 to 0.93 of FASER&rsquo;s own simulated flux in five "
    "of the six bins, which <i>is</i> the non-uniformity.<br><br>"
    "<b>The scale band is POWHEG-V2&rsquo;s own</b> (user, 2026-09-08), from "
    "reweighting every ladder point at seven scales and redoing the whole "
    "fold at each rather than propagating a band by hand. It is 1 to 2% on "
    "these yields &mdash; far below the flux uncertainty, which is the same "
    "ordering the Run 3 + Run 4 projections show.")


def faser_electronic_table():
    """The per-bin numbers behind the figure."""
    # the CURRENT POWHEG-V2 ladder (ubexcess_correct 1, genuine neutrons),
    # the one paper plot 12b draws (2026-09-21)
    p = f"{RESULTS_NU}/faser_electronic.json"
    q = f"{BASE}/data/faser_electronic/faser_electronic.json"
    if not (os.path.exists(p) and os.path.exists(q)):
        return ('<p class="missing">results_nu/faser_electronic.json not '
                'yet generated &mdash; run analysis/faser_powheg_rates.py '
                'and analysis/faser_electronic.py</p>')
    with open(p) as f:
        d = json.load(f)
    with open(q) as f:
        dat = json.load(f)
    # the CONF-note block is not part of the public release: its two columns
    # are left out without it
    prl, conf = dat["prl"], dat.get("conf")
    out = ['<div class="tablewrap"><table><thead><tr>'
           '<th>bin</th><th>E<sub>&nu;</sub> [GeV]</th>'
           '<th>&lang;&sigma;&rang; POWHEG-V2</th>'
           '<th>&lang;&sigma;&rang; GENIE</th>'
           '<th>&lang;&sigma;&rang; FASER</th>'
           '<th>our flux / FASER</th>'
           '<th>N<sub>fid</sub> POWHEG-V2</th>'
           '<th>N<sub>fid</sub> GENIE</th>'
           + ('<th>FASER sim.</th><th>FASER data</th>' if conf else '') +
           '<th>acceptance, ours / FASER</th>'
           '</tr></thead><tbody>']
    for i, r in enumerate(d["rows"]):
        rng = (f"{r['e_lo']:.0f}&ndash;{r['e_hi']:.0f}" if r["e_hi"]
               else f"&gt;{r['e_lo']:.0f}")
        sign = {"nu": "&nu;", "nubar": "&nu;&#772;",
                "both": "&nu;+&nu;&#772;"}[r["charge"]]
        out.append(
            f'<tr><td>{sign}</td><td>{rng}</td>'
            f'<td class="num">{r["sigma_1e38_cm2_per_nucleon"]:.0f}</td>'
            f'<td class="num">{r["genie_sigma_1e38_cm2_per_nucleon"]:.0f}</td>'
            f'<td class="num">{prl["sigma_1e38_cm2_per_nucleon"][i]:.1f} '
            f'&plusmn; {prl["sigma_err"][i]:.1f}</td>'
            f'<td class="num">{r["flux_1e6_fb_cm2"]/prl["flux_1e6_fb_cm2"][i]:.2f}</td>'
            f'<td class="num">{r["n_fid_conf"]:.0f}</td>'
            f'<td class="num">{r["genie_n_fid_conf"]:.0f}</td>'
            + (f'<td class="num">{conf["sim_nominal"][i]:.0f} &plusmn; '
               f'{conf["sim_total_unc"][i]:.0f}</td>'
               f'<td class="num">{conf["fig5"][i]:.0f}</td>' if conf else '') +
            f'<td class="num">{r["acceptance_pct"]:.1f}% / '
            f'{prl["acceptance"][i]:.1f}%</td></tr>')
    out.append('</tbody></table></div>')
    out.append(
        '<p class="secnote">&lang;&sigma;&rang; is the flux-weighted mean '
        'charged-current cross-section per nucleon on tungsten, in '
        '10<sup>&minus;38</sup> cm&sup2;; N<sub>fid</sub> is the number of '
        'interactions in the 100 mm-radius fiducial cylinder at '
        '186 fb<sup>&minus;1</sup>. <b>The acceptance column compares two '
        'different things on purpose.</b> Ours is the kinematic acceptance '
        'of FASER&rsquo;s own muon cuts &mdash; E&prime; &gt; 100 GeV and '
        '&theta; &lt; 25 mrad, which is exactly this benchmark&rsquo;s '
        'Tier S &mdash; with no y window, so that it is a fraction of the '
        'same thing theirs is. Theirs additionally carries the radial cuts '
        'and the detector geometry, which this benchmark does not model: '
        'above 300 GeV the ratio settles at a nearly energy-independent '
        '0.57 to 0.63, which is what that geometry costs.</p>')
    return "".join(out)


def faser_electronic_section():
    return ('<h2>Comparison with FASER&rsquo;s electronic detector</h2>'
            + f'<p class="secnote">{FASER_ELEC_NOTE}</p>'
            # THE PAPER PLOT ITSELF (user, 2026-09-21: promoted to paper
            # plot 12b), so the sub-tab cannot drift from it
            + '<div class="grid one">' + figure_or_note_nu(
                "pp_faser_electronic", "",
                "Left: the flux-weighted charged-current cross-section per "
                "nucleon in each &minus;L/E<sub>&nu;</sub> bin, ours against "
                "the value FASER publishes. Right: the number of "
                "&nu;<sub>&mu;</sub> interactions in the fiducial volume at "
                "186 fb<sup>&minus;1</sup>, ours against FASER&rsquo;s "
                "simulation and against their unfolded data. The shaded band "
                "is POWHEG-V2&rsquo;s own 7-point scale envelope. Lower "
                "panels: ratios. The neutrino half of the axis is on the "
                "left of the dotted line and the antineutrino half on the "
                "right. This is paper plot 12b "
                "(analysis/paper_plots/pp12b_faser_electronic.py).") + '</div>'
            + faser_electronic_table())


# ---- DONUT (user, 2026-09-10) ------------------------------------------
DONUT_NOTE = (
    "<b>DONUT is the experiment that first observed the tau neutrino</b> "
    "(arXiv:0711.0728), in a beam-dump beam at Fermilab on a target of 1&nbsp;mm "
    "stainless-steel sheets interleaved with nuclear emulsion. Its Fig.&nbsp;9 "
    "gives the number of charged particles at the primary vertex of all 578 "
    "neutrino interactions it located, beside its own LEPTO-based simulation, "
    "and that is what is reproduced here. The two figures the comparison needs "
    "are read off the paper&rsquo;s vector graphics rather than digitised by "
    "eye: the arXiv PDF stores the histogram outlines and the marker centres "
    "as the numbers that were plotted."
    "<br><br>"
    "<b>The multiplicity at the primary vertex is flavour blind, which is what "
    "makes this possible at all.</b> DONUT&rsquo;s beam is 57% "
    "&nu;<sub>&mu;</sub>, 38% &nu;<sub>e</sub> and 5% &nu;<sub>&tau;</sub>, "
    "and this benchmark has no electron- or tau-neutrino samples. It does not "
    "need them. Every charged-current event puts exactly <i>one</i> lepton "
    "track at the vertex &mdash; an electron, a muon or a tau, and the tau "
    "decays microns downstream so the emulsion counts the parent &mdash; and "
    "the hadronic side does not know which lepton was made. The "
    "user&rsquo;s own argument for the cross-sections, that they are flavour "
    "independent, therefore carries over to the multiplicity for a reason of "
    "its own. What does not carry over is the tau mass, which suppresses the "
    "&nu;<sub>&tau;</sub> cross-section at low y and so tilts the hadronic "
    "energy of the 5% of the sample that is &nu;<sub>&tau;</sub>."
    "<br><br>"
    "<b>The band on each generator is the emulsion acceptance, not a theory "
    "uncertainty.</b> DONUT counts tracks reconstructed in emulsion: a track "
    "has to leave the steel plate it was born in and be seen in the next "
    "emulsion layer, and the automated scanning takes a limited angular "
    "range. The paper states no acceptance, so the prediction is computed "
    "with three of them &mdash; a momentum threshold of 0.1 or 0.3 GeV/c and "
    "tan&thinsp;&theta; below 0.5 or 1.0 &mdash; and the envelope is drawn. "
    "It is worth about one track in the mean, more than the difference "
    "between the two generators, so a comparison that ignored it would read "
    "the acceptance as physics.")

DONUT_LIMITS = (
    "<b>What each generator can and cannot say here</b>, stated in place "
    "rather than left to look like a discrepancy. GENIE runs on an iron "
    "nucleus with its full channel list &mdash; charged and neutral current, "
    "quasi-elastic, resonance, deep inelastic and coherent &mdash; and with "
    "final-state interactions, so it is the entry that can populate the "
    "low-multiplicity bins at all. POWHEG-V2 is a perturbative "
    "charged-current DIS calculation above Q<sup>2</sup> = 2.25 GeV<sup>2</sup> "
    "on a free nucleon with no final-state interaction; it has no "
    "quasi-elastic channel, no resonances and no neutral current, and those "
    "three are exactly what the n<sub>ch</sub> = 1 and 2 bins are made of. "
    "Its curve is a statement about the hadronisation of the DIS part and "
    "should be read as one."
    "<br><br>"
    "<b>The GENIE multiplicity predates the 2026-10-01 charm fix</b> "
    "(patches/genie-aivazis-charm-propagator.diff); the cross-section table "
    "below is recomputed with it and moves by 1e-5. The rerun of the events "
    "stopped on the HepMC converter, which now refuses nuclear events without "
    "a struck nucleon; at DONUT&rsquo;s energies the fix lowers the charm "
    "cross-section by under 1%, which is under 0.05% of the events."
    "<br><br>"
    "<b>And neither prediction has a detector in it.</b> DONUT&rsquo;s own "
    "simulation carries the trigger, the scan and a location efficiency that "
    "is 0.31 for shower-like events against 0.77 for the rest &mdash; so it "
    "depends on the topology and therefore on the multiplicity itself. Ours "
    "carry none of that, and every curve is normalised to the same number of "
    "events, so only shapes are compared. Where our curve and "
    "DONUT&rsquo;s simulation differ, that is physics; where their simulation "
    "and their data differ, it is either physics or their detector and we "
    "cannot tell which from the outside.")


def donut_result():
    p = f"{RESULTS_NU}/donut_nch.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def donut_table():
    """Mean multiplicity and shape chi-square, FORMATTED from the result."""
    d = donut_result()
    if d is None:
        return ('<p class="missing">donut_nch.json not yet generated '
                '&mdash; run analysis/donut_nch.py --merge</p>')
    rows = [("DONUT data", f"{d['data_mean']:.2f}", "&mdash;", "&mdash;"),
            ("DONUT simulation (LEPTO, with their detector)",
             f"{d['donut_mc_mean']:.2f}",
             f"{d['donut_mc_chi2']/d['ndf']:.1f}", "&mdash;")]
    for key, g in d["generators"].items():
        lab = {"genie": ("GENIE G18_02a on iron, all channels, with "
                         "final-state interactions"),
               "powheg": "POWHEG-V2 + Pythia&nbsp;8, free-nucleon iron, CC DIS"}
        v, raw = g["variants"]["default"], g["variants"]["raw"]
        rows.append((lab.get(key, key), f"{v['mean']:.2f}",
                     f"{v['chi2']/d['ndf']:.1f}", f"{raw['mean']:.2f}"))
    out = ['<div class="tablewrap"><table><thead><tr><th>sample</th>'
           '<th>&lang;n<sub>ch</sub>&rang;</th>'
           '<th>&chi;<sup>2</sup>/N vs the data</th>'
           '<th>&lang;n<sub>ch</sub>&rang; with no acceptance</th>'
           '</tr></thead><tbody>']
    for r in rows:
        out.append('<tr><td>%s</td><td class="num">%s</td>'
                   '<td class="num">%s</td><td class="num">%s</td></tr>' % r)
    out.append('</tbody></table></div>')
    out.append(
        f'<p class="secnote">Means are over the drawn range '
        f'n<sub>ch</sub> = 1&ndash;13. The &chi;<sup>2</sup> uses the '
        f'data&rsquo;s Poisson error alone on {d["ndf"]} degrees of freedom '
        f'(one is spent on the common normalisation), so it orders the curves '
        f'rather than measuring a fit quality &mdash; no theory or acceptance '
        f'uncertainty enters it. The last column is the generator&rsquo;s bare '
        f'charged multiplicity, which no emulsion measures; the difference '
        f'from the previous column is what the acceptance is worth.</p>')
    return "".join(out)


DONUT_XSEC_NOTE = (
    "<b>The cross-section DONUT actually published, on the same target.</b> "
    "Writing the charged-current cross-section per nucleon as "
    "&sigma;<sub>&nu;</sub>(E) = &sigma;<sub>const</sub> &times; E &times; "
    "K(E), with K the kinematic suppression from the tau-lepton mass and "
    "safely 1 for electrons and muons over their energy range, DONUT measured "
    "&sigma;<sub>const</sub> for the tau neutrino. <b>The tau mass is divided "
    "out of their number</b>, which is what makes the user&rsquo;s point &mdash; "
    "that the interaction cross-section is flavour independent &mdash; enough "
    "to compare it against a muon-neutrino calculation. Ours is GENIE on "
    "iron, charged current, averaged over neutrino and antineutrino because "
    "DONUT&rsquo;s beam is: they measure the ratio of the two fluxes as "
    "1.05 &plusmn; 0.13 and take it as one.")


def donut_xsec_table():
    """The cross-section comparison, FORMATTED from the result file."""
    p = f"{RESULTS_NU}/donut_xsec.json"
    if not os.path.exists(p):
        return ('<p class="missing">donut_xsec.json not yet generated '
                '&mdash; run analysis/donut_xsec.py</p>')
    with open(p) as f:
        d = json.load(f)
    dn = d["donut"]
    rows = [
        ("DONUT, &nu;<sub>&tau;</sub> (their Eq. 16)",
         f"{dn['sigma_const_tau']:.2f} &plusmn; {dn['stat']:.2f} "
         f"&plusmn; {dn['syst']:.2f}", "&mdash;"),
        ("the &nu;<sub>&mu;</sub> + &nu;&#773;<sub>&mu;</sub> average they "
         "compare it against",
         f"{dn['reference_numu_average']:.2f}", "&mdash;"),
        ("<b>this benchmark: GENIE on iron, flux averaged over their "
         "Fig.&nbsp;2</b>",
         f"{d['sigma_const_flux_averaged']:.3f}",
         f"{d['ratio_to_reference']:.3f}"),
        ("the same at 100 GeV", f"{d['sigma_const_at_100gev']:.3f}",
         f"{d['sigma_const_at_100gev']/dn['reference_numu_average']:.3f}"),
    ]
    out = ['<div class="tablewrap"><table><thead><tr><th>&sigma;<sub>const</sub>'
           ' [10<sup>&minus;38</sup> cm<sup>2</sup>/GeV per nucleon]</th>'
           '<th>value</th><th>ratio to the reference</th>'
           '</tr></thead><tbody>']
    for r in rows:
        out.append('<tr><td>%s</td><td class="num">%s</td>'
                   '<td class="num">%s</td></tr>' % r)
    out.append('</tbody></table></div>')
    out.append(
        f'<p class="secnote">Our flux-averaged value sits '
        f'{d["ratio_to_donut_tau"]:.2f} times DONUT&rsquo;s tau-neutrino '
        f'measurement, {abs(d["pull_vs_donut_tau"]):.1f} standard deviations '
        f'from it once their statistical and systematic errors are combined '
        f'&mdash; their measurement carries a 60% uncertainty, so this is a '
        f'consistency statement rather than a test. Against the muon-neutrino '
        f'number they compare it to it agrees to '
        f'{100*abs(d["ratio_to_reference"]-1):.1f}%, which is also the check '
        f'that the nuclear splines were summed per nucleus and not per '
        f'nucleon: no wrong convention lands within a factor of two.</p>')
    return "".join(out)


def donut_section():
    return ('<h2>DONUT: charged-particle multiplicity in the emulsion</h2>'
            + f'<p class="secnote">{DONUT_NOTE}</p>'
            + '<div class="grid one">' + figure_or_note_nu(
                "cmp_donut_nch", "",
                "Charged particles at the primary vertex of DONUT&rsquo;s 578 "
                "located events. Points are the data with their Poisson "
                "errors, the dashed grey line is DONUT&rsquo;s own LEPTO "
                "simulation with their detector in it, and the coloured lines "
                "are this benchmark&rsquo;s generators folded with "
                "DONUT&rsquo;s Fig.&nbsp;2 spectrum of interacting neutrinos "
                "(&ldquo;FSI&rdquo; in the legend is final-state "
                "interactions). "
                "Bands are the emulsion acceptance, not a theory uncertainty. "
                "Every curve carries the same number of events. "
                "<span class=\"mono\">analysis/donut_nch.py</span>.")
            + '</div>'
            + donut_table()
            + f'<p class="secnote">{DONUT_LIMITS}</p>'
            + '<h2>The tau-neutrino cross-section</h2>'
            + f'<p class="secnote">{DONUT_XSEC_NOTE}</p>'
            + donut_xsec_table())


def data_panel():
    """Comparison against a MEASUREMENT, as its own top-level tab.

    Moved out of the "Physics studies" sub-tabs on 2026-09-01 (user): it was
    the ninth of nine there, and it is the only place on the page where the
    benchmark meets data rather than another calculation.

    SPLIT INTO SUB-TABS on 2026-09-08 (user: "the Comparison with Data tab
    can be split in subtabs, now we have two, one for NOMAD, one for FASER
    emulsion detector").  One measurement per sub-tab, so a third is one
    entry rather than a rewrite.
    """
    nomad = ('<h2>NOMAD dimuon production</h2>'
            + f'<p class="secnote">{NOMAD_NOTE}</p>'
            + '<div class="grid one">' + figure_or_note_nu(
                "cmp_nomad_dimuon", "",
                "The dimuon fraction against neutrino energy: NOMAD data "
                "with statistical and systematic errors combined, against "
                "YADISM at NLO in FONLL with its seven-point scale band. "
                "Horizontal bars are the bin widths, since the calculation "
                "is evaluated at the bin centre. Lower panel: data over "
                "theory.") + '</div>'
            + nomad_table())
    # >>> THE PUBLISHED-RATE REPRODUCTIONS ARE NOT HERE ANY MORE (user,
    # 2026-09-08: "the text after 'Published rate predictions, reproduced
    # with our own machinery' is repeated from other tabs and should be
    # deleted I guess?  In this subtab we care only about NOMAD"). <<<  They
    # were rendered twice, here and on the "Predictions for FASER" tab, on the
    # argument that reproducing a published prediction is also a comparison
    # with data.  It is not: those two check our machinery against another
    # CALCULATION, which is what the rest of this page does, while this tab is
    # the one place the benchmark meets a MEASUREMENT.  They stay on the FASER
    # tab, beside the rates they bear on, and only there.
    # THE FASER EMULSION SUB-TAB WENT WITH THE earlier RESULTS (user, 2026-09-19):
    # it read the 0.2 < y < 0.9 rates; its successor is the FASER emulsion
    # figure on the "Paper plots" tab.
    return tabs([("nomad", "NOMAD dimuons", ""),
                 ("faserelec", "FASER electronic detector", ""),
                 # the tau-neutrino experiment, and the only measurement here
                 # that is a hadronic FINAL STATE rather than a rate
                 ("donut", "DONUT multiplicity", "")],
                "measurement", ns="data-", cls="sel",
                bodies={"nomad": nomad,
                        "faserelec": faser_electronic_section(),
                        "donut": donut_section()})


def audit_stamps():
    """Check every result's `energy_gev` against the energy its NAME claims.

    The report is the one place that reads all ~340 result JSONs at once, so
    it is the cheapest place to run this.  A file called `..._400GeV.json`
    holding a result stamped 1000.0 is the signature of an input path being
    made energy-aware while its output path was left fixed -- which has
    happened here twice, once corrupting the 1 TeV anchor and once (in
    combine_nlo.py) mislabelling all six Herwig NLO scan points.  Neither
    errored; both were caught by this stamp.  Printing the audit on every
    build means a recurrence is visible immediately rather than whenever
    someone next happens to look.
    """
    import re
    bad = []
    # Both result families, because both carry an energy tag in the FILENAME
    # and an energy stamp INSIDE, and it is the disagreement between those two
    # that this catches.  nlo_unc_* was added later and uses exactly the same
    # naming rule, so leaving it out would have exempted the newest files from
    # the audit that exists because this trap recurs.
    prefixes = ("histos_", "nlo_unc_")
    for resdir in (RESULTS, RESULTS_NU):
        for fn in sorted(os.listdir(resdir)):
            pre = next((p for p in prefixes if fn.startswith(p)), None)
            if pre is None or not fn.endswith(".json"):
                continue
            name = fn[len(pre):-len(".json")]
            # THE Q2 FLOOR SUFFIX SITS AFTER THE ENERGY TAG (2026-09-01), so
            # a plain endswith() reads histos_sherpa_400GeV_q2min11.json as
            # untagged and reports twenty perfectly good files as
            # mislabelled.  Twenty false alarms is how an audit stops being
            # read -- the same reason a NaN-vs-NaN diff is filtered --
            # so the suffix is stripped before the tag is looked for rather
            # than the audit being narrowed.
            m = re.search(r"_q2min[0-9.]+$", name)
            if m:
                name = name[:m.start()]
            # ANY energy tag, not only the three scan energies (2026-09-02):
            # the analytic sigma_fid(E) ladder behind the FASER flux
            # convolution writes _15GeV, _632GeV, _3TeV and a dozen more.
            # Matching only the scan energies made every one of them read as
            # untagged, i.e. as the anchor, and the audit reported 26 correct
            # files as mislabelled -- which is how an audit stops being read.
            claimed = beams.ANCHOR_ENERGY
            m = re.search(r"_(\d+(?:\.\d+)?)(GeV|TeV)$", name)
            if m:
                claimed = float(m.group(1)) * (1000.0 if m.group(2) == "TeV"
                                               else 1.0)
            try:
                with open(f"{resdir}/{fn}") as f:
                    st = json.load(f).get("energy_gev")
            except ValueError:
                continue
            if st is not None and abs(float(st) - claimed) > 1e-6:
                bad.append(f"{resdir}/{fn}: name says {claimed:g} GeV, "
                           f"stamp says {st:g} GeV")
    if bad:
        print("*** ENERGY STAMP MISMATCH -- results are mislabelled ***",
              file=sys.stderr)
        for b in bad:
            print("   " + b, file=sys.stderr)
    else:
        print("energy stamps: all results agree with their filenames")
    return bad


def check_claims():
    """Re-verify the numbers written into the blurbs, on every build.

    A standalone tool nobody runs is the same as no tool -- the same lesson
    the insufficient-statistics flag taught when it was stored and then
    ignored by the renderer.  Failures are reported, not fatal: the page is
    still worth having while a sentence is being corrected.
    """
    try:
        sys.path.insert(0, f"{BASE}/tools")
        import check_report_claims
        if check_report_claims.main():
            print("*** report prose disagrees with results/ (above) ***",
                  file=sys.stderr)
        # A GROWN sample leaves every listed mtime untouched, so the manifest
        # tripwire in analyze.py cannot see it; this one counts job
        # directories on disk against the manifest (2026-09-04).
        import check_manifests
        if check_manifests.main([]):
            print("*** a result reads fewer job files than its sample has "
                  "(above) -- the page carries a partial-sample number ***",
                  file=sys.stderr)
        # AND THAT EVERY PAPER PLOT IS WHERE ITS PAPER_SECTION SAYS.  The
        # figures are copied into paper/figures/ by a glob, so one can be
        # produced, synced and tracked while no \includegraphics ever names
        # it -- and the reverse has happened too, a figure the user asked to
        # hold back staying in the document.  Both directions are checked.
        import check_paper_figures_used
        if check_paper_figures_used.main():
            print("*** a paper plot is not where its PAPER_SECTION says "
                  "(above) ***", file=sys.stderr)
    except Exception as exc:            # never fail the report over the check
        print(f"WARNING: could not verify report claims: {exc}",
              file=sys.stderr)


def main():
    audit_stamps()
    check_claims()
    # ONLY THE FINAL REGION SURVIVES (user, 2026-09-19): "all results with the
    # old 0.2 < y < 0.9 cut can be deleted (from file and from report), only
    # survives from now on."  The per-energy tabs, the energy scan, the
    # physics studies, the earlier paper plots and the extra-plots tab went with
    # their results; the paper plots are now the landing tab.
    e_tabs, bodies = [], {}
    e_tabs.append(("paper2", "Paper plots", ""))
    bodies["paper2"] = paper_plots_panel()
    # COMPARISON WITH DATA AND THE FASER PREDICTIONS ARE TOP-LEVEL TABS
    # (user, 2026-09-01).
    e_tabs.append(("data", "Comparison with data", ""))
    bodies["data"] = data_panel()
    e_tabs.append(("faser", "Predictions for FASER", ""))
    bodies["faser"] = faser_rates_panel()
    # SIDIS AT FASER, its own top-level tab (user, 2026-09-21), to be handed
    # to theory collaborators.
    e_tabs.append(("sidis", "SIDIS at FASER", ""))
    bodies["sidis"] = sidis_panel()
    # the fragmentation functions a calculation of those yields could use
    # (user, 2026-09-21)
    e_tabs.append(("ff", "Fragmentation functions", ""))
    bodies["ff"] = ff_panel()
    top = tabs(e_tabs, "report section", ns="e-", cls="sel", bodies=bodies)
    html = f"""<title>FASER DIS Generator Benchmark</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>{_check_css(CSS)}</style>
<div class="wrap">
<header>
  <p class="eyebrow">evgen-benchmark &middot; results</p>
  <h1>FASER DIS Generator Benchmark</h1>
  <p class="setup">Event-generator comparison at FASER kinematics, against
  analytic YADISM structure functions (ZM-VFNS, with target-mass
  corrections). A tungsten target of 74 protons and 110 neutrons with
  NNPDF4.0 NNLO, Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no cut on y,
  &mu;<sub>F</sub> = &mu;<sub>R</sub> = Q, Gmu EW scheme, massless charm,
  no QED radiation, hadrons stable if c&tau; &gt; 10 mm.
  Beam energies from 400 GeV to 4 TeV, two currents.</p>
  <p class="byline">Paulina Hernandez-Sainz, Jelle Koorn and Juan Rojo<br>
  <span class="affil">Department of Physics and Astronomy, Vrije Universiteit
  Amsterdam, and Nikhef Theory Group, Amsterdam, The Netherlands</span></p>
  <p class="auxnote">Auxiliary material for the paper
  <i>Event generation for neutrino and muon DIS at the LHC</i>.</p>
</header>
{top}
<footer class="mono">FASER DIS generator benchmark &middot; Q&sup2; &gt;
4 GeV&sup2;, W &gt; 3 GeV, tungsten &middot; 400 GeV &ndash; 4 TeV &middot;
NNPDF4.0, &mu;<sub>F</sub> = &mu;<sub>R</sub> = Q, ZM-VFNS</footer>
</div>
<dialog id="zoom"><img id="zoomimg" src="" alt=""></dialog>
<script>{image_script()}{JS}</script>
"""
    with open(OUT, "w") as f:
        f.write(html)
    size_mb = os.path.getsize(OUT) / 1e6
    print(f"wrote {OUT} ({size_mb:.1f} MB)")

    # A second copy OUTSIDE the repo, to hand to someone else.  report.html is
    # gitignored (tens of MB of embedded figures, fully regenerable), so the
    # repo is not where a shareable copy can live.  Destination comes from
    # config.sh like every other path -- override with BENCH_REPORT_COPY, or
    # pass --no-copy to skip it.
    if "--no-copy" in sys.argv:
        return
    dest = paths.REPORT_COPY
    try:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copyfile(OUT, dest)
    except OSError as exc:
        # Never fail the report over the copy -- report.html is already written.
        print(f"WARNING: could not write the shareable copy to {dest}: {exc}",
              file=sys.stderr)
    else:
        print(f"wrote {dest} ({size_mb:.1f} MB)  <- open or send this one")


if __name__ == "__main__":
    main()
