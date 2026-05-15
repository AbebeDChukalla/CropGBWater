"""Render the animated social-preview GIF (1200x630, 5-frame cycle).

Each frame shows the key finding from one dashboard module:
  1. Overview  — total water +9%
  2. Atlas     — Asia ≈ 52% share
  3. Trends    — green water 84%
  4. Crops     — cereals 47%
  5. Countries — India + China + USA

Reuses the brand-mark drawing routine from make_social_preview.py so the
logo looks identical across PNG (LinkedIn unfurl) and GIF (manual upload).
"""

from __future__ import annotations
import json
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# Import the brand-mark helper from the static-preview generator
import importlib.util, sys
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("sp", HERE / "make_social_preview.py")
sp = importlib.util.module_from_spec(spec); sys.modules["sp"] = sp; spec.loader.exec_module(sp)

OUT  = HERE.parent / "social-preview.gif"
DATA = HERE.parent / "data"
W, H = 1200, 630

PAPER, PAPER2 = sp.PAPER, sp.PAPER2
INK, INK70, INK40 = sp.INK, sp.INK70, sp.INK40
GREEN, GREEN_SOFT, BLUE, BLUE_SOFT = sp.GREEN, sp.GREEN_SOFT, sp.BLUE, sp.BLUE_SOFT

# Fonts — reuse the same fonts as the PNG generator
serif_xl = sp.serif_xl
serif_lg = sp.serif_lg
serif_md = sp.serif_md
mono_md  = sp.mono_md
mono_eyebrow    = sp.mono_eyebrow
mono_eyebrow_sm = sp.mono_eyebrow_sm
sans_brand_bold = sp.sans_brand_bold
mono_stack      = sp.mono_stack


def draw_logo_block(d, im):
    """Brand mark + CropGBWater wordmark + eyebrow at the top of the frame."""
    sp.draw_brand_mark(d, im, 110, 64, size=72)
    g_img, g_center, paste_y, b_img, b_center, _ = sp.draw_logo(d, 160, 50)
    gw, gh = g_img.size
    bw, bh = b_img.size
    im.paste(g_img, (g_center - gw // 2, paste_y), g_img)
    im.paste(b_img, (b_center - bw // 2, paste_y), b_img)
    eyebrow_y = 120
    d.line([(80, eyebrow_y + 12), (134, eyebrow_y + 12)], fill=INK70, width=2)
    d.text((146, eyebrow_y),
           "ATLAS OF AGRICULTURAL GREEN- AND BLUE WATER USE  ·  2026",
           font=mono_eyebrow, fill=INK70)


def draw_footer(d):
    d.line([(80, 586), (1120, 586)], fill=(20, 20, 15, 26), width=1)
    d.text((80, 600),
           "46 crops  ·  172 countries  ·  6 continents  ·  2000-2020",
           font=mono_md, fill=INK70)
    url = "abebedchukalla.github.io/CropGBWater"
    url_w = d.textlength(url, font=mono_md)
    d.text((1120 - url_w, 600), url, font=mono_md, fill=INK70)


def draw_step_dots(d, idx, n):
    """5 progress dots at the very bottom right of the headline area."""
    base_x = 80
    base_y = 200
    for i in range(n):
        on = (i == idx)
        d.rectangle((base_x + i * 36, base_y,
                     base_x + i * 36 + (28 if on else 20), base_y + 4),
                    fill=INK if on else (20, 20, 15, 38))


def draw_dual_stat_bar(d, y, green_pct, blue_pct, green_val_km3, blue_val_km3):
    bar_x, bar_w, bar_h = 80, 1040, 84
    split = int(bar_w * green_pct / 100)
    d.rectangle((bar_x, y, bar_x + split - 4, y + bar_h), fill=GREEN)
    d.rectangle((bar_x + split, y, bar_x + bar_w, y + bar_h), fill=BLUE)
    d.ellipse((bar_x + 14, y + 14, bar_x + 24, y + 24), fill=GREEN_SOFT)
    d.text((bar_x + 32, y + 14), "GREEN WATER USE", font=mono_eyebrow, fill=PAPER)
    d.text((bar_x + 14, y + 36),
           f"{green_pct:.1f}%   {green_val_km3:,.0f} km³/yr",
           font=serif_md, fill=PAPER)
    bx = bar_x + split + 12
    d.ellipse((bx, y + 12, bx + 8, y + 20), fill=BLUE_SOFT)
    d.text((bx + 14, y + 10), "BLUE USE", font=mono_eyebrow_sm, fill=PAPER)
    d.text((bx, y + 32), f"{blue_pct:.1f}%", font=sp.serif_sm, fill=PAPER)
    d.text((bx, y + 64), f"{blue_val_km3:,.0f} km³/yr",
           font=mono_eyebrow_sm, fill=PAPER)


def render_frame(idx, *, summary, continents, groups, countries, change):
    """Build one 1200x630 frame for the given module index (0-4)."""
    im = Image.new("RGBA", (W, H), PAPER + (255,))
    d  = ImageDraw.Draw(im)
    draw_logo_block(d, im)
    draw_step_dots(d, idx, 5)
    draw_footer(d)

    g20 = summary["global"]["2020"]
    g10 = summary["global"]["2010"]

    if idx == 0:
        # OVERVIEW — headline +9%
        y = 230
        d.text((80, y), "The world’s farms drink", font=serif_lg, fill=INK); y += 92
        pct = f"{round(change['total_pct'])}%"
        d.text((80, y), pct, font=serif_xl, fill=BLUE)
        pct_w = d.textlength(pct, font=serif_xl)
        d.text((80 + pct_w + 18, y), "more water than", font=serif_lg, fill=INK); y += 92
        d.text((80, y), "they did ", font=serif_lg, fill=INK)
        w = d.textlength("they did ", font=serif_lg)
        d.text((80 + w, y), "a decade ago.", font=serif_xl, fill=INK40); y += 86
        draw_dual_stat_bar(d, y,
                           g20["green_km3"] / g20["total_km3"] * 100,
                           g20["blue_km3"]  / g20["total_km3"] * 100,
                           g20["green_km3"], g20["blue_km3"])

    elif idx == 1:
        # ATLAS — Asia
        asia = next((c for c in continents if c["name"] == "Asia"), None)
        share = asia["share_pct"] if asia else 52
        countries_count = asia["n_countries"] if asia else 46
        kml = round(asia["total_km3"]) if asia else 3525
        y = 220
        d.text((80, y), "Asia consumes", font=serif_lg, fill=INK); y += 92
        s = f"{share:.0f}%"
        d.text((80, y), s, font=serif_xl, fill=GREEN)
        w = d.textlength(s, font=serif_xl)
        d.text((80 + w + 18, y), "of all crop water,", font=serif_lg, fill=INK); y += 92
        d.text((80, y), "across ", font=serif_lg, fill=INK)
        w = d.textlength("across ", font=serif_lg)
        s2 = f"{countries_count} "
        d.text((80 + w, y), s2, font=serif_xl, fill=BLUE)
        w2 = d.textlength(s2, font=serif_xl)
        d.text((80 + w + w2, y), "countries.", font=serif_lg, fill=INK40); y += 78
        d.text((80, y), f"≈ {kml:,} km³/yr · more than Africa + Europe + Americas combined.",
               font=mono_eyebrow, fill=INK70)

    elif idx == 2:
        # TRENDS — green share + growth
        green_pct = g20["green_km3"] / g20["total_km3"] * 100
        y = 220
        d.text((80, y), "Of every drop of crop water,", font=serif_lg, fill=INK); y += 92
        s = f"{green_pct:.1f}%"
        d.text((80, y), s, font=serif_xl, fill=GREEN)
        w = d.textlength(s, font=serif_xl)
        d.text((80 + w + 18, y), "is rainfall —", font=serif_lg, fill=INK); y += 92
        d.text((80, y), "irrigation barely grew.", font=serif_xl, fill=INK40); y += 86
        d.text((80, y),
               f"Green water consumption grew {change['green_pct']:.1f}%  ·  blue grew only {change['blue_pct']:.1f}%  ·  2010 → 2020",
               font=mono_eyebrow, fill=INK70)

    elif idx == 3:
        # CROPS — cereals share
        cereals = next((g for g in groups if g["id"] == "cereals"), None)
        share = cereals["share_pct"] if cereals else 47
        y = 220
        d.text((80, y), "Cereals alone account for", font=serif_lg, fill=INK); y += 92
        s = f"{share:.0f}%"
        d.text((80, y), s, font=serif_xl, fill=(201, 136, 59))
        w = d.textlength(s, font=serif_xl)
        d.text((80 + w + 18, y), "of crop water.", font=serif_lg, fill=INK); y += 92
        d.text((80, y), "Rice draws ", font=serif_lg, fill=INK)
        w = d.textlength("Rice draws ", font=serif_lg)
        d.text((80 + w, y), "1,000+ km³/yr.", font=serif_xl, fill=BLUE); y += 78
        d.text((80, y),
               "Wheat · Maize · Rice · Sorghum · Barley · Millets — together ≈ half of global crop water.",
               font=mono_eyebrow, fill=INK70)

    elif idx == 4:
        # COUNTRIES — top 3
        top = countries[:3]
        names = ", ".join(c["name"] for c in top)
        tot = sum(c["total_km3"] for c in top)
        share = round(tot / g20["total_km3"] * 100)
        y = 220
        d.text((80, y), "India, China, USA together", font=serif_lg, fill=INK); y += 92
        s = f"~{share}%"
        d.text((80, y), s, font=serif_xl, fill=BLUE)
        w = d.textlength(s, font=serif_xl)
        d.text((80 + w + 18, y), "of all crop water.", font=serif_lg, fill=INK); y += 92
        d.text((80, y), "Each draws ", font=serif_lg, fill=INK)
        w = d.textlength("Each draws ", font=serif_lg)
        d.text((80 + w, y), "400–1,000 km³/yr.", font=serif_xl, fill=INK40); y += 78
        d.text((80, y),
               f"172 countries modelled across 6 continents · click any country in the dashboard for details.",
               font=mono_eyebrow, fill=INK70)

    return im


def main():
    summary    = json.loads((DATA / "summary.json").read_text())
    continents = json.loads((DATA / "continents.json").read_text())["continents"]
    groups     = json.loads((DATA / "crop_groups.json").read_text())["groups"]
    countries  = json.loads((DATA / "countries.json").read_text())["countries"]
    change     = summary["change_2010_2020"]

    frames = [render_frame(i, summary=summary, continents=continents,
                            groups=groups, countries=countries, change=change)
              for i in range(5)]

    # Convert to P mode for smaller GIF and quantize for a clean palette
    frames_p = [f.convert("RGB").quantize(colors=128, dither=Image.FLOYDSTEINBERG) for f in frames]

    # Save as animated GIF.  Each frame: 3.0s. Total cycle: 15s.
    frames_p[0].save(
        OUT,
        save_all=True,
        append_images=frames_p[1:],
        duration=3000,          # ms per frame
        loop=0,                 # infinite loop
        optimize=True,
        disposal=2,
    )
    print(f"Wrote {OUT}  ({OUT.stat().st_size // 1024} KB, {len(frames_p)} frames, "
          f"{sum(3.0 for _ in frames_p):.0f}s cycle)")


if __name__ == "__main__":
    main()
