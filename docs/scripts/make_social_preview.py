"""Render the LinkedIn social-preview PNG (1200x630) from real data.

Standalone Pillow renderer — does not depend on a SVG converter. Re-run any
time the headline numbers change.
"""

from __future__ import annotations
import json
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parents[1] / "social-preview.png"
DATA = Path(__file__).resolve().parents[1] / "data" / "summary.json"

W, H = 1200, 630
PAPER  = (243, 239, 230)
PAPER2 = (235, 230, 216)
INK    = (20, 20, 15)
INK70  = (74, 72, 66)
INK40  = (138, 135, 125)
GREEN  = (45, 90, 61)
GREEN_SOFT = (126, 163, 107)
BLUE   = (30, 74, 122)
BLUE_SOFT = (111, 161, 207)

# Try the best system fonts. Windows ships Times.ttf and Consola.ttf.
# We deliberately use Times for serif so we don't depend on any installed font.
def _font(name_candidates, size):
    for name in name_candidates:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


serif_xl = _font(["timesbi.ttf", "Times New Roman Bold Italic.ttf", "timesi.ttf", "times.ttf"], 76)
serif_lg = _font(["times.ttf", "Times New Roman.ttf"], 76)
serif_md = _font(["times.ttf"], 44)
serif_sm = _font(["times.ttf"], 34)
mono_sm  = _font(["consola.ttf", "Consolas.ttf", "cour.ttf"], 13)
mono_md  = _font(["consola.ttf"], 14)
mono_eyebrow = _font(["consolab.ttf", "consola.ttf"], 22)   # bumped for visibility on social cards
sans_brand = _font(["seguibl.ttf", "segoeuib.ttf", "segoeui.ttf"], 22)


def text(d, xy, s, font, fill):
    d.text(xy, s, font=font, fill=fill)


def main():
    summary = json.loads(DATA.read_text())
    g20 = summary["global"]["2020"]
    change = summary["change_2010_2020"]
    total = g20["total_km3"]
    green = g20["green_km3"]
    blue  = g20["blue_km3"]
    green_pct = green / total * 100
    blue_pct  = blue  / total * 100

    im = Image.new("RGB", (W, H), PAPER)
    d = ImageDraw.Draw(im)

    # ── Top brand bar ──
    # brand mark: two overlapping circles
    d.ellipse((80, 36, 102, 58), fill=GREEN)
    d.ellipse((92, 36, 114, 58), fill=BLUE)
    text(d, (124, 36), "CropsGreenBlueWater", sans_brand, INK)

    # eyebrow rule + label  (longer + thicker rule to match bigger type)
    eyebrow_y = 86
    d.line([(80, eyebrow_y + 12), (134, eyebrow_y + 12)], fill=INK70, width=2)
    text(d, (146, eyebrow_y),
         "ATLAS OF AGRICULTURAL GREEN- AND BLUE WATER USE  ·  2026",
         mono_eyebrow, INK70)

    # ── Editorial headline (3 lines) ──
    # Headline uses the paper's rounded figure ("9%") even though the precise
    # 2010→2020 change is {:.1f}% — the dashboard itself shows full precision.
    y = 200
    # Line 1 — "The world's farms drink"
    text(d, (80, y), "The world’s farms drink", serif_lg, INK)
    y += 92
    # Line 2 — "9% more water than"
    pct_str = f"{round(change['total_pct'])}%"
    text(d, (80, y), pct_str, serif_xl, BLUE)
    pct_w = d.textlength(pct_str, font=serif_xl)
    text(d, (80 + pct_w + 18, y), "more water than", serif_lg, INK)
    y += 92
    # Line 3 — "they did a decade ago." — italic muted
    text(d, (80, y), "they did ", serif_lg, INK)
    they_w = d.textlength("they did ", font=serif_lg)
    text(d, (80 + they_w, y), "a decade ago.", serif_xl, INK40)
    y += 88

    # ── Green / blue split bar ──
    bar_x, bar_w, bar_h = 80, 1040, 90
    split = int(bar_w * green_pct / 100)
    d.rectangle((bar_x, y, bar_x + split - 4, y + bar_h), fill=GREEN)
    d.rectangle((bar_x + split, y, bar_x + bar_w, y + bar_h), fill=BLUE)

    # green-side labels
    d.ellipse((bar_x + 14, y + 14, bar_x + 24, y + 24), fill=GREEN_SOFT)
    text(d, (bar_x + 32, y + 14), "GREEN WATER  ·  RAINFALL",
         mono_eyebrow, PAPER)
    text(d, (bar_x + 14, y + 36),
         f"{green_pct:.1f}%   {green:,.0f} km³/yr",
         serif_md, PAPER)

    # blue-side labels
    bx = bar_x + split + 14
    d.ellipse((bx, y + 14, bx + 10, y + 24), fill=BLUE_SOFT)
    text(d, (bx + 18, y + 14), "BLUE WATER", mono_eyebrow, PAPER)
    text(d, (bx, y + 36), f"{blue_pct:.1f}%", serif_sm, PAPER)

    # ── Bottom strip ──
    d.line([(80, 586), (1120, 586)], fill=(20, 20, 15, 26), width=1)
    text(d, (80, 600),
         f"46 crops  ·  172 countries  ·  2000-2020  ·  0.083° (~10 km) grid",
         mono_md, INK70)
    url = "abebedchukalla.github.io/CropGBWater"
    url_w = d.textlength(url, font=mono_md)
    text(d, (1120 - url_w, 600), url, mono_md, INK70)

    im.save(OUT, "PNG", optimize=True)
    print(f"Wrote {OUT}  ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
