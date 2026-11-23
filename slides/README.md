# Presentation slides

Slides for the presentation at BRIN: the current state of the IL–CO₂ solubility model and the proposed computational roadmap (18 slides).

| Path | What it is |
|---|---|
| `index.html` | Standalone version of the deck. Open it in a browser; speaker notes are under each slide. Print to PDF from the browser to get one slide per page. |
| `source/` | The deck source as edited in the Claude Slides viewer: `deck.json` (slide order), `deck.css` (styles) and `slides/*.html` (one file per slide, speaker notes in `<aside>`). |
| `build_slides.py` | Rebuilds `index.html` from `source/`. |

`index.html` loads its two figures from `../images/`, so open it from inside this repository. If the browser blocks local images, serve the repo root:

```bash
python3 -m http.server 8000
```

and open http://localhost:8000/slides/index.html.

To update after editing the slides in the viewer, replace the files in `source/` and run:

```bash
python3 slides/build_slides.py
```
