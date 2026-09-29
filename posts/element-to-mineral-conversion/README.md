# From elements to minerals

`index.qmd` is the Quarto source for the public post. The site build renders this source to its normal public page.

The 29 September 2026 review revises the literature discussion, mathematical qualifications and teaching implementations. Review findings and validation limits are recorded in [`../../reviews/emc-review-2026-09-29.md`](../../reviews/emc-review-2026-09-29.md). From the repository root, run `uv run python scripts/check_emc_article.py` to execute the main numerical examples; the optional NumPyro snippet is deliberately skipped.

[`rendered-report.html`](rendered-report.html) is a retained **historical, pre-review snapshot**. It has not been regenerated from the revised source and must not be treated as its current rendered output. The existing site configuration also copies this snapshot to `posts/element-to-mineral-conversion/rendered-report.html`.

Before publishing the revision, render and inspect `index.qmd` in the locked repository environment. Do not silently replace the historical snapshot with an unvalidated render.
