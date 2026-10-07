#!/usr/bin/env python3
"""Shared look of the self-contained HTML pages: stylesheet, tab switcher
script, WebP figure embedding.

Extracted on 2026-10-07 from analysis/make_report.py when the benchmark's
HTML report was retired (user: "we don't need it anymore"); the one page
still built with it is the neutrino-telescope study,
telescopes/build_report.py.  The page is read locally from file://, so
everything it needs is inline.
"""
import base64
import io


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
