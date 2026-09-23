# SANGHA Research Systems — static landing page v1

A dependency-free static site. It can be served directly from GitHub Pages, Cloudflare Pages, Netlify, or any ordinary web host.

## Files
- `index.html` — all HTML, CSS, and motion code in one file
- `assets/sangha-mark.svg` — canonical vertical-proportion mark supplied by the user

## Local preview
From this folder:

```bash
python3 -m http.server 8080
```

Then open `http://localhost:8080`.

## Deployment
Upload the folder as-is. No build step is required.

## Design notes
- Canonical vertical mark proportions are preserved.
- Bauhaus-brutalist grid and typography.
- Fine-line canvas geometry responds subtly to pointer and scroll.
- `prefers-reduced-motion` is respected.
- The two crossbarless A forms are rendered as `Λ`, echoing logical conjunction / recursive “AND”.

## Field Geometry

The self-contained Python/Make source is in [`software/field-geometry`](software/field-geometry).
The public project page is `public/projects/field-geometry.html`; immutable v0.1.0
downloads are under `public/downloads/field-geometry/`.

The deployed Cloudflare Worker serves `public/`, as specified by `wrangler.jsonc`.
Preview the deployed surface with `python3 -m http.server --directory public 8080`.
