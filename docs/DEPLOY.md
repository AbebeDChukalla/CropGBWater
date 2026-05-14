# Deploying to GitHub Pages

Goal: serve `webapp/` at **https://abebedchukalla.github.io/CropGBWater/**
from your existing repo, then share that URL on LinkedIn.

GitHub Pages can only serve from the repo root or a `/docs` folder when
"Deploy from a branch" is selected. The simplest path is therefore to
**push the dashboard as a `docs/` folder on `main`**, alongside your
existing Notebooks_v1.0/ code.

---

## One-time setup

Open PowerShell and run these from `D:\Project\G2CWD\CropGBWater\`:

```powershell
# 1. Clone your repo (skip if you already have a local clone)
git clone https://github.com/AbebeDChukalla/CropGBWater.git CropGBWater-repo
cd CropGBWater-repo

# 2. Copy the dashboard into docs/ (note: docs/, not webapp/, for GH Pages)
xcopy /E /I /Y "..\Combined_WF_Outputs\webapp\*" "docs\"

# 3. Stage, commit, push
git add docs
git commit -m "Add interactive CropsGreenBlueWater dashboard"
git push origin main
```

## Enable Pages (one click)

Open https://github.com/AbebeDChukalla/CropGBWater/settings/pages and set:

- **Source:** Deploy from a branch
- **Branch:** `main`  ·  Folder: `/docs`
- **Save**

Wait ~30 seconds, then open:

**https://abebedchukalla.github.io/CropGBWater/**

The first build can take a minute. If the page shows a 404, the build is
still pending — check the *Actions* tab on the repo. Subsequent pushes to
`main` auto-redeploy.

---

## Updating the dashboard later

```powershell
cd D:\Project\G2CWD\CropGBWater\CropGBWater-repo

# If the source CSVs changed, regenerate JSON first:
python ..\Combined_WF_Outputs\webapp\scripts\build_static_data.py

# Then sync the docs/ folder
xcopy /E /I /Y "..\Combined_WF_Outputs\webapp\*" "docs\"
git add docs
git commit -m "Refresh dashboard data"
git push
```

---

## Verifying the LinkedIn unfurl

After the page is live, paste the URL into the
**LinkedIn Post Inspector**: https://www.linkedin.com/post-inspector/

It re-fetches your page and shows what card LinkedIn will render.
If the preview looks stale, click "Inspect" again — LinkedIn caches
unfurls for ~7 days.

The card should show:

- **Title:** CropsGreenBlueWater · Atlas of agricultural water use
- **Description:** Global atlas across 46 crops, 172 countries…
- **Image:** the 8.7%-headline preview (1200×630)

If the image looks broken, confirm `social-preview.png` actually exists at
`https://abebedchukalla.github.io/CropGBWater/social-preview.png` (open
the URL directly in a browser).

---

## Alternative: subdirectory or other host

| Want | Do |
|---|---|
| Keep dashboard on its own branch (no `docs/` folder on `main`) | Push `webapp/` contents to a `gh-pages` branch, set Pages source to `gh-pages` / `/(root)` |
| Custom domain like `cropgbw.org` | After Pages is live, add a `CNAME` file with the domain inside `docs/`, then point DNS A records to GitHub Pages IPs |
| Vercel / Netlify instead | Drag the `webapp/` folder onto netlify.com/drop (instant) or run `npx vercel` from inside `webapp/` |
