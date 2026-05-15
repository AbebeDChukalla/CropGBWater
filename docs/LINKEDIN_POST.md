# LinkedIn post drafts

Three drafts in different tones. Pick one, edit lightly, paste into a
new LinkedIn post.

## Two preview images, two ways to use them

- **`social-preview.png`** — used by LinkedIn's automatic unfurl when you
  paste the dashboard URL. You don't have to attach it manually; LinkedIn
  fetches it from the page's `og:image` tag.
- **`social-preview.gif`** — a 5-frame, 15-second animation that cycles
  through the key finding of each module (Overview → Atlas → Trends →
  Crops → Countries). LinkedIn does not animate `og:image` files, so to
  show this animation you need to **attach it as an image to the post**
  manually. The card unfurl will still show, but with the GIF attached
  directly LinkedIn will play the animation in the feed.

The card preview LinkedIn shows is driven by the `<meta og:*>` tags in
`index.html` plus `social-preview.png`. To force a refresh after deploy:
https://www.linkedin.com/post-inspector/

---

## Draft 1 — Editorial / story-led (recommended)

> The world's farms drink 9% more water than they did a decade ago — and
> almost all of that new demand was met by rainfall, not irrigation.
>
> That single finding sits at the heart of our new paper on global crop
> water consumption, modelled across 46 crops at 10-km resolution for
> 2000, 2010, and 2020. To make the dataset easier to explore, I built an
> interactive dashboard around it: country rankings, crop breakdowns,
> trend chart, and a map you can click into for any country's full crop
> mix.
>
> Two numbers that stuck with me while building it:
> • Rice alone accounts for 1,035 km³/yr — about 15% of all crop water.
> • India and China each draw nearly 1,000 km³/yr — together a fifth of
>   the global total.
>
> Try the dashboard → https://abebedchukalla.github.io/CropGBWater/
>
> Code, methodology and the underlying CSV/NetCDF outputs are all open
> (CC BY 4.0). Feedback and ideas for what to add next are welcome.
>
> #WaterFootprint #AgriculturalWater #OpenScience #DataVisualization #Hydrology #FoodSecurity

---

## Draft 2 — Short / punchy

> New dashboard: global crop water consumption across 46 crops and 172
> countries, 2000-2020.
>
> Headline: total crop water use rose 9% in a decade — and 84% of that
> water is green (rainfall absorbed by crops), not blue (irrigation).
>
> Click any country for its crop-by-crop breakdown.
>
> https://abebedchukalla.github.io/CropGBWater/
>
> Built from the dataset behind our paper on spatially explicit crop
> water consumption. All code and aggregated outputs are open.
>
> #WaterFootprint #AgriculturalWater #DataVisualization

---

## Draft 3 — Academic / collaboration-led

> Sharing an interactive companion to our paper "Global spatially explicit
> crop water consumption shows an overall increase of 9% for 46 crops
> from 2010 to 2020."
>
> The dashboard lets researchers and policymakers explore the same model
> outputs without writing a line of code:
>
> – Country-level green / blue water choropleth (2000, 2010, 2020)
> – Per-country detail with top-12 crops and 3-year trend
> – Ranking of all 172 modelled countries
> – Crop composition across 10 groups and 46 individual crops
>
> Dashboard → https://abebedchukalla.github.io/CropGBWater/
> Source code & data → https://github.com/AbebeDChukalla/CropGBWater
>
> Built on the SPAM + GAEZ + PCR-GLOBWB stack. Feedback from the
> hydrology, agronomy, and policy communities is very welcome — especially
> on which additional views would be most useful.
>
> #WaterResources #CropWaterFootprint #PCRGLOBWB #SPAM #OpenScience

---

## Posting checklist

- [ ] Deploy (`DEPLOY.md`) and confirm https://abebedchukalla.github.io/CropGBWater/ loads
- [ ] Open the URL directly: confirm map, donut chart and country click work
- [ ] Run https://www.linkedin.com/post-inspector/ on the URL → confirm
      preview shows the 9%-headline image
- [ ] **For animation:** download `social-preview.gif` from the deployed site
      (`https://abebedchukalla.github.io/CropGBWater/social-preview.gif`) and
      attach it directly to your LinkedIn post — LinkedIn will play the 15-second
      cycle of key findings inline in the feed.
- [ ] Tag co-authors if you want them notified — paste their LinkedIn @handles
- [ ] Post between Tue–Thu, 9–11am local time for best reach
