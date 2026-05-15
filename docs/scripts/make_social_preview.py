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
from PIL import Image, ImageDraw, ImageFont, ImageFilter, ImageChops

OUT  = Path(__file__).resolve().parents[1] / "social-preview.png"
DATA = Path(__file__).resolve().parents[1] / "data" / "summary.json"

W, H = 1200, 630
PAPER  = (243, 239, 230)
PAPER2 = (235, 230, 216)
INK    = (20, 20, 15)
INK70  = (74, 72, 66)
INK40  = (138, 135, 125)
GREEN  = (45, 90, 61)         # wordmark / split-bar colour
GREEN_SOFT = (126, 163, 107)
BLUE   = (30, 74, 122)        # wordmark / split-bar colour
BLUE_SOFT = (111, 161, 207)
GREY_OVERLAP = (216, 214, 207)
# Brighter, more vivid colours just for the globe brand mark
GLOBE_GREEN      = (62, 145, 80)
GLOBE_GREEN_DARK = (22, 80, 40)
GLOBE_BLUE       = (40, 110, 180)
GLOBE_BLUE_DARK  = (16, 56, 110)

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


def draw_brand_mark(d, im, cx, cy, size=88):
    """Two overlapping glossy globes (green left, blue right) with world-map
    grid lines, a white drop at the centre overlap, and a shiny plant inside.

    Rendered at 4x supersampling for crisp edges, then downsampled.
    """
    SS = 4
    # Brand mark aspect is wider than tall (two side-by-side globes)
    W2 = int(size * 1.45)
    H2 = size
    layer = Image.new("RGBA", (W2 * SS, H2 * SS), (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)

    cx_l = W2 * SS // 2
    cy_l = H2 * SS // 2
    r    = int(size * SS * 0.40)
    # Wider separation so each globe shows a clear coloured crescent on its side
    gx   = cx_l - int(r * 0.65)
    bx   = cx_l + int(r * 0.65)

    # Soft drop shadow under each globe
    shadow = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.ellipse((gx - r, cy_l - r + 5 * SS, gx + r, cy_l + r + 5 * SS), fill=(20, 20, 15, 80))
    sd.ellipse((bx - r, cy_l - r + 5 * SS, bx + r, cy_l + r + 5 * SS), fill=(20, 20, 15, 80))
    shadow = shadow.filter(ImageFilter.GaussianBlur(radius=5 * SS))
    layer.alpha_composite(shadow)

    def _gradient_globe(center_x, base_color, dark_color):
        """Saturated filled circle + opaque inner highlight + thin dark ring."""
        # Solid fill in vivid colour
        ld.ellipse((center_x - r, cy_l - r, center_x + r, cy_l + r), fill=base_color)

        # Lighter inner ellipse — gives a 3D shaded sphere feel
        light = tuple(min(255, c + 55) for c in base_color)
        ld.ellipse((center_x - int(r * 0.55), cy_l - int(r * 0.62),
                    center_x + int(r * 0.10), cy_l - int(r * 0.18)),
                    fill=light)

        # Bright specular dot (top-left) — solid white-ish for a "shining" feel
        ld.ellipse((center_x - int(r * 0.34), cy_l - int(r * 0.46),
                    center_x - int(r * 0.10), cy_l - int(r * 0.28)),
                    fill=(255, 255, 255))

        # Thin darker outline ring for the globe edge
        ld.ellipse((center_x - r, cy_l - r, center_x + r, cy_l + r),
                   outline=dark_color, width=max(2, int(1.6 * SS)))

    def _grid_lines(center_x):
        """White world-map grid lines drawn directly with semi-transparent fill.

        We render onto a transparent layer, then alpha_composite (which blends
        instead of replacing).  A circular mask clips lines that leak outside.
        """
        gl = Image.new("RGBA", layer.size, (0, 0, 0, 0))
        gd = ImageDraw.Draw(gl)
        lw = max(1, int(0.8 * SS))
        col = (255, 255, 255, 170)
        # Parallels (horizontal-ish arcs spanning the full globe width)
        for dy in (-int(r * 0.55), -int(r * 0.28), 0, int(r * 0.28), int(r * 0.55)):
            gd.arc((center_x - r, cy_l - r + dy - int(r * 0.05),
                    center_x + r, cy_l + r + dy + int(r * 0.05)),
                   start=0, end=180, fill=col, width=lw)
        # Meridians (vertical ellipses)
        for rx in (int(r * 0.25), int(r * 0.55)):
            gd.ellipse((center_x - rx, cy_l - r,
                        center_x + rx, cy_l + r),
                       outline=col, width=lw)
        # Central vertical line
        gd.line([(center_x, cy_l - r), (center_x, cy_l + r)],
                fill=col, width=lw)

        # Clip grid layer to the globe circle (multiply its alpha by mask)
        mask = Image.new("L", layer.size, 0)
        ImageDraw.Draw(mask).ellipse(
            (center_x - r, cy_l - r, center_x + r, cy_l + r), fill=255)
        rgba = gl.split()
        new_alpha = ImageChops.multiply(rgba[3], mask)
        gl_clipped = Image.merge("RGBA", (rgba[0], rgba[1], rgba[2], new_alpha))
        layer.alpha_composite(gl_clipped)

    _gradient_globe(gx, base_color=GLOBE_GREEN, dark_color=GLOBE_GREEN_DARK)
    _grid_lines(gx)
    _gradient_globe(bx, base_color=GLOBE_BLUE,  dark_color=GLOBE_BLUE_DARK)
    _grid_lines(bx)

    # ── White drop / leaf shape at centre overlap ──
    # Smaller so the green / blue globe crescents stay clearly visible
    drop_w = int(r * 0.62)
    drop_h = int(r * 1.05)
    drop_layer = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    drop_d = ImageDraw.Draw(drop_layer)
    # Drop body: vertical leaf via pieslice + polygon
    drop_d.ellipse((cx_l - drop_w // 2, cy_l - drop_h // 2,
                    cx_l + drop_w // 2, cy_l + drop_h // 2),
                   fill=(248, 248, 246), outline=(0, 0, 0, 38), width=max(1, int(0.7 * SS)))
    # Subtle wet sheen on upper-left of drop
    sheen = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    sh = ImageDraw.Draw(sheen)
    sh.ellipse((cx_l - int(drop_w * 0.45), cy_l - int(drop_h * 0.42),
                cx_l - int(drop_w * 0.05), cy_l - int(drop_h * 0.05)),
               fill=(255, 255, 255, 180))
    sheen = sheen.filter(ImageFilter.GaussianBlur(radius=2 * SS))
    layer.alpha_composite(drop_layer)
    layer.alpha_composite(sheen)

    # ── Plant inside drop ──
    stem_top    = cy_l - int(drop_h * 0.28)
    stem_bottom = cy_l + int(drop_h * 0.25)
    ld.line([(cx_l, stem_top), (cx_l, stem_bottom)],
            fill=(45, 130, 25), width=int(1.8 * SS))

    leaf_outline = (15, 60, 15)
    leaf_fill    = (110, 200, 70)
    leaf_bright  = (200, 240, 160)

    # Left leaf (rounded teardrop tilted up-left)
    left_pts = [
        (cx_l,                              stem_top + int(r * 0.18)),
        (cx_l - int(r * 0.42),              stem_top + int(r * 0.08)),
        (cx_l - int(r * 0.48),              stem_top - int(r * 0.16)),
        (cx_l - int(r * 0.10),              stem_top - int(r * 0.02)),
    ]
    ld.polygon(left_pts, fill=leaf_fill, outline=leaf_outline)
    ld.line(left_pts + [left_pts[0]], fill=leaf_outline,
            width=int(1.0 * SS), joint="curve")
    # left highlight
    ld.line([(cx_l - int(r * 0.32), stem_top - int(r * 0.08)),
             (cx_l - int(r * 0.15), stem_top + int(r * 0.0))],
            fill=leaf_bright, width=int(1.3 * SS))

    # Right leaf (higher, tilted up-right)
    right_pts = [
        (cx_l,                              stem_top - int(r * 0.05)),
        (cx_l + int(r * 0.42),              stem_top - int(r * 0.18)),
        (cx_l + int(r * 0.48),              stem_top - int(r * 0.42)),
        (cx_l + int(r * 0.10),              stem_top - int(r * 0.28)),
    ]
    ld.polygon(right_pts, fill=leaf_fill, outline=leaf_outline)
    ld.line(right_pts + [right_pts[0]], fill=leaf_outline,
            width=int(1.0 * SS), joint="curve")
    # right highlight
    ld.line([(cx_l + int(r * 0.18), stem_top - int(r * 0.30)),
             (cx_l + int(r * 0.36), stem_top - int(r * 0.20))],
            fill=leaf_bright, width=int(1.3 * SS))

    # Downsample and composite
    final = layer.resize((W2, H2), Image.LANCZOS)
    im.alpha_composite(final, (cx - W2 // 2, cy - H2 // 2))


def draw_logo(d, x, y):
    """CropGBWater wordmark: charcoal Crop + green G + blue B + charcoal Water.
    Small "Green" / "Blue" labels in their colors are placed ABOVE the G and B.
    Returns positions so the caller can paint the labels in a second pass."""
    parts = [("Crop", INK), ("G", GREEN), ("B", BLUE), ("Water", INK)]
    cx = x
    positions = {}
    for s, color in parts:
        w = d.textlength(s, font=sans_brand_bold)
        text(d, (cx, y), s, sans_brand_bold, color)
        positions[s] = (cx, cx + w)
        cx += w

    gx0, gx1 = positions["G"]
    bx0, bx1 = positions["B"]
    g_center = int((gx0 + gx1) / 2)
    b_center = int((bx0 + bx1) / 2)
    return g_center, b_center


def draw_logo_labels(d, g_center, b_center, label_y):
    """Paint 'Green' (in green) + 'Blue' (in blue) as a single phrase centred
    above the G/B pair, since the two letters sit adjacent in the wordmark."""
    center_x = (g_center + b_center) // 2
    g_w = d.textlength("Green", font=mono_eyebrow_sm)
    b_w = d.textlength("Blue",  font=mono_eyebrow_sm)
    gap = 8
    total_w = g_w + gap + b_w
    start_x = center_x - total_w // 2
    d.text((start_x,                label_y), "Green", font=mono_eyebrow_sm, fill=GREEN)
    d.text((start_x + g_w + gap,    label_y), "Blue",  font=mono_eyebrow_sm, fill=BLUE)


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
    # Brand mark is wider (two side-by-side globes); place it on the left and
    # shift the wordmark right past its full width.
    mark_size  = 130
    mark_w     = int(mark_size * 1.45)
    mark_cx, mark_cy = 80 + mark_w // 2, 78
    draw_brand_mark(d, im, mark_cx, mark_cy, size=mark_size)
    # Wordmark starts after the brand mark + a small gap
    wordmark_x = mark_cx + mark_w // 2 + 22
    g_center, b_center = draw_logo(d, wordmark_x, 56)
    # "Green" / "Blue" labels centred above the G / B pair
    draw_logo_labels(d, g_center, b_center, label_y=40)

    # ── Eyebrow ──
    eyebrow_y = 124
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
