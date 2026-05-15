"""Render the LinkedIn social-preview PNG (1200x630).

Brand: CropGBWater with bold G(green)/B(blue) + stacked vertical Green/Blue labels.
Mark:  Two overlapping circles (green + blue) with radial sheen, light-grey
       intersection, sharper plant icon (thicker strokes, outlined leaves,
       highlight on the water drop) centred in the overlap.
"""

from __future__ import annotations
import json
import math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont, ImageFilter

OUT  = Path(__file__).resolve().parents[1] / "social-preview.png"
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
GREY_OVERLAP = (216, 214, 207)

def _font(names, size):
    for n in names:
        try:
            return ImageFont.truetype(n, size)
        except OSError:
            continue
    return ImageFont.load_default()


serif_xl = _font(["timesbi.ttf", "timesi.ttf", "times.ttf"], 76)
serif_lg = _font(["times.ttf"], 76)
serif_md = _font(["times.ttf"], 44)
serif_sm = _font(["times.ttf"], 34)
mono_md  = _font(["consola.ttf"], 14)
mono_eyebrow = _font(["consolab.ttf", "consola.ttf"], 22)
mono_eyebrow_sm = _font(["consolab.ttf", "consola.ttf"], 14)
sans_brand_bold = _font(["seguibl.ttf", "segoeuib.ttf", "arialbd.ttf"], 30)
mono_stack = _font(["consola.ttf"], 13)


def text(d, xy, s, font, fill):
    d.text(xy, s, font=font, fill=fill)


def draw_brand_mark(d, im, cx, cy, size=72):
    """Two overlapping circles, light-grey intersection, sharp plant in middle.

    Renders on a 4× supersampled canvas then downscales for crisp edges.
    """
    SS = 4  # supersample factor for sharp circles + plant
    W2, H2 = size * 3, size * 3
    layer = Image.new("RGBA", (W2 * SS, H2 * SS), (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)

    # Centre inside the supersampled layer
    px = W2 * SS // 2
    py = H2 * SS // 2
    r = (size * SS) // 2
    gx, bx = px - r // 2, px + r // 2

    # Soft shadow under both circles
    shadow = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.ellipse((gx - r, py - r + 6 * SS, gx + r, py + r + 6 * SS), fill=(20, 20, 15, 55))
    sd.ellipse((bx - r, py - r + 6 * SS, bx + r, py + r + 6 * SS), fill=(20, 20, 15, 55))
    shadow = shadow.filter(ImageFilter.GaussianBlur(radius=6 * SS))
    layer.alpha_composite(shadow)

    # Filled circles
    ld.ellipse((gx - r, py - r, gx + r, py + r), fill=GREEN)
    ld.ellipse((bx - r, py - r, bx + r, py + r), fill=BLUE)
    # Light-grey intersection (oval centred between the two circle centres)
    overlap_w = int(r * 1.0)
    overlap_h = int(r * 1.45)
    ld.ellipse((px - overlap_w // 2, py - overlap_h // 2,
                px + overlap_w // 2, py + overlap_h // 2),
               fill=GREY_OVERLAP)
    # Bright highlight on top of each circle for a "shining" finish
    hl = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    hd = ImageDraw.Draw(hl)
    hd.ellipse((gx - r + int(r * 0.15), py - r + int(r * 0.15),
                gx - r + int(r * 0.70), py - r + int(r * 0.55)),
               fill=(255, 255, 255, 90))
    hd.ellipse((bx - r + int(r * 0.15), py - r + int(r * 0.15),
                bx - r + int(r * 0.70), py - r + int(r * 0.55)),
               fill=(255, 255, 255, 90))
    hl = hl.filter(ImageFilter.GaussianBlur(radius=3 * SS))
    layer.alpha_composite(hl)

    # ── Plant icon (sharper: thicker strokes, outlined leaves, glow on drop)
    # Coordinates in supersampled space; py is the centre.
    stroke = int(2.4 * SS)
    leaf_outline = int(1.6 * SS)

    # Stem
    ld.line([(px, py + int(size * 0.34) * SS // 2),
             (px, py - int(size * 0.34) * SS // 2)],
            fill=(13, 51, 32), width=stroke)

    # Leaves
    def leaf(points, fill, outline):
        ld.polygon(points, fill=fill, outline=outline)
        # second-pass outline for a thicker visible edge
        ld.line(points + [points[0]], fill=outline, width=leaf_outline, joint="curve")

    # Left leaf
    leaf([
        (px,                            py - int(0.02 * size * SS)),
        (px - int(size * 0.20 * SS),    py - int(size * 0.06 * SS)),
        (px - int(size * 0.22 * SS),    py - int(size * 0.20 * SS)),
        (px - int(size * 0.04 * SS),    py - int(size * 0.10 * SS)),
    ], fill=GREEN_SOFT, outline=(13, 51, 32))

    # Right leaf (higher up the stem)
    leaf([
        (px,                            py - int(size * 0.10 * SS)),
        (px + int(size * 0.20 * SS),    py - int(size * 0.14 * SS)),
        (px + int(size * 0.22 * SS),    py - int(size * 0.30 * SS)),
        (px + int(size * 0.04 * SS),    py - int(size * 0.20 * SS)),
    ], fill=GREEN_SOFT, outline=(13, 51, 32))

    # Water drop on top of plant with bright highlight
    drop_cy = py - int(size * 0.38 * SS)
    rx, ry = int(size * 0.05 * SS), int(size * 0.07 * SS)
    ld.ellipse((px - rx, drop_cy - ry, px + rx, drop_cy + ry),
               fill=BLUE_SOFT, outline=(22, 58, 98), width=int(1.6 * SS))
    # Glint
    ld.ellipse((px - rx // 2, drop_cy - ry,
                px - rx // 2 + rx // 2, drop_cy - ry // 4),
               fill=(228, 238, 248, 220))

    # Downsample with LANCZOS for crisp edges
    final = layer.resize((W2, H2), Image.LANCZOS)
    # Crop to actual size centred
    cropped = final.crop(((W2 - size) // 2, (H2 - size) // 2,
                         (W2 + size) // 2, (H2 + size) // 2))
    # Composite onto the base image at (cx, cy)
    im.alpha_composite(cropped.convert("RGBA"), (cx - size // 2, cy - size // 2))


def draw_logo(d, x, y):
    """CropGBWater wordmark with bold G(green) + B(blue) + stacked vertical labels."""
    parts = [("Crop", INK), ("G", GREEN), ("B", BLUE), ("Water", INK)]
    cx = x
    positions = {}
    for s, color in parts:
        w = d.textlength(s, font=sans_brand_bold)
        text(d, (cx, y), s, sans_brand_bold, color)
        positions[s] = (cx, cx + w)
        cx += w

    # Vertical stacked labels under G and B
    gx0, gx1 = positions["G"]
    bx0, bx1 = positions["B"]
    g_center = int((gx0 + gx1) / 2)
    b_center = int((bx0 + bx1) / 2)

    # Rotate text 90° via a transparent overlay
    def stacked(text_str, center_x, color):
        bbox = mono_stack.getbbox(text_str)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        txt = Image.new("RGBA", (tw + 4, th + 4), (0, 0, 0, 0))
        td = ImageDraw.Draw(txt)
        td.text((2, 2), text_str, font=mono_stack, fill=color)
        rotated = txt.rotate(90, expand=True, resample=Image.BICUBIC)
        rw, rh = rotated.size
        return rotated, rw, rh

    g_img, gw, gh = stacked("GREEN", g_center, GREEN)
    b_img, bw, bh = stacked("BLUE",  b_center, BLUE)
    paste_y = y + 38   # under the letter row
    return g_img, g_center, paste_y, b_img, b_center, paste_y


def main():
    summary = json.loads(DATA.read_text())
    g20 = summary["global"]["2020"]
    change = summary["change_2010_2020"]
    total = g20["total_km3"]
    green = g20["green_km3"]
    blue  = g20["blue_km3"]
    green_pct = green / total * 100
    blue_pct  = blue  / total * 100

    im = Image.new("RGBA", (W, H), PAPER + (255,))
    d = ImageDraw.Draw(im)

    # ── Brand mark + logo ──
    draw_brand_mark(d, im, 110, 64, size=72)
    parts = draw_logo(d, 160, 50)
    g_img, g_center, paste_y, b_img, b_center, _ = parts
    # paste vertical labels for G and B
    gw, gh = g_img.size
    bw, bh = b_img.size
    im.paste(g_img, (g_center - gw // 2, paste_y), g_img)
    im.paste(b_img, (b_center - bw // 2, paste_y), b_img)

    # ── Eyebrow ──
    eyebrow_y = 120
    d.line([(80, eyebrow_y + 12), (134, eyebrow_y + 12)], fill=INK70, width=2)
    text(d, (146, eyebrow_y),
         "ATLAS OF AGRICULTURAL GREEN- AND BLUE WATER USE  ·  2026",
         mono_eyebrow, INK70)

    # ── Editorial headline ──
    y = 220
    text(d, (80, y), "The world’s farms drink", serif_lg, INK)
    y += 92
    pct_str = f"{round(change['total_pct'])}%"
    text(d, (80, y), pct_str, serif_xl, BLUE)
    pct_w = d.textlength(pct_str, font=serif_xl)
    text(d, (80 + pct_w + 18, y), "more water than", serif_lg, INK)
    y += 92
    text(d, (80, y), "they did ", serif_lg, INK)
    they_w = d.textlength("they did ", font=serif_lg)
    text(d, (80 + they_w, y), "a decade ago.", serif_xl, INK40)
    y += 86

    # ── Green / Blue use split bar ──
    bar_x, bar_w, bar_h = 80, 1040, 84
    split = int(bar_w * green_pct / 100)
    d.rectangle((bar_x, y, bar_x + split - 4, y + bar_h), fill=GREEN)
    d.rectangle((bar_x + split, y, bar_x + bar_w, y + bar_h), fill=BLUE)

    # green-side labels
    d.ellipse((bar_x + 14, y + 14, bar_x + 24, y + 24), fill=GREEN_SOFT)
    text(d, (bar_x + 32, y + 14), "GREEN WATER USE", mono_eyebrow, PAPER)
    text(d, (bar_x + 14, y + 36),
         f"{green_pct:.1f}%   {green:,.0f} km³/yr",
         serif_md, PAPER)

    # blue-side labels (smaller font + shorter wording to fit narrow strip)
    bx = bar_x + split + 12
    d.ellipse((bx, y + 12, bx + 8, y + 20), fill=BLUE_SOFT)
    text(d, (bx + 14, y + 10), "BLUE USE", mono_eyebrow_sm, PAPER)
    text(d, (bx, y + 32), f"{blue_pct:.1f}%", serif_sm, PAPER)
    text(d, (bx, y + 64), f"{blue:,.0f} km³/yr", mono_eyebrow_sm, (243,239,230,180))

    # ── Bottom strip ──
    d.line([(80, 586), (1120, 586)], fill=(20, 20, 15, 26), width=1)
    text(d, (80, 600),
         f"46 crops  ·  172 countries  ·  2000-2020  ·  0.083° (~10 km) grid",
         mono_md, INK70)
    url = "abebedchukalla.github.io/CropGBWater"
    url_w = d.textlength(url, font=mono_md)
    text(d, (1120 - url_w, 600), url, mono_md, INK70)

    im.convert("RGB").save(OUT, "PNG", optimize=True)
    print(f"Wrote {OUT}  ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
