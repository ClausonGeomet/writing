# Writing

This is the source repository for a public Quarto writing site: essays on geometallurgy, mineral inference, and defensible technical decisions.

The first post, [From elements to minerals: what an assay can tell us](posts/element-to-mineral-conversion/index.qmd), is accompanied by an [archived rendered report](posts/element-to-mineral-conversion/rendered-report.html). The source and archived report are versioned together so the publication record is auditable.

## Local setup and rendering

Install [uv](https://docs.astral.sh/uv/) and [Quarto](https://quarto.org/docs/get-started/). Then create the local environment from the committed lockfile:

```bash
uv sync --frozen
uv run quarto render
```

Preview while writing with:

```bash
uv run quarto preview
```

The Python environment contains only packages imported by executable cells in the QMD: NumPy, SciPy, pandas, plotnine, IPython, and Jupyter execution support. It intentionally excludes the optional sampler examples because their cells use `eval: false`.

For those optional examples, add only the toolchain you plan to use:

```bash
uv add --group bayes arviz pymc
# or
uv add --group bayes arviz numpyro jax
# Bambi is useful for the separate calibration example.
uv add --group bayes bambi
```

JAX hardware-specific installation should follow the [JAX installation guide](https://docs.jax.dev/installation.html).

## Add a post

1. Start with `drafts/<slug>/index.qmd` and add `_metadata.yml` (or front matter) containing `title`, `description`, `author`, `date`, `categories`, and `draft: true`.
2. Preview the draft directly while working. Drafts are excluded from the normal site render and deployment.
3. When ready, move the directory to `posts/<slug>/`, change to `draft: false`, and run `uv run quarto render`.
4. Review `_site/` locally, commit the source and any deliberately retained rendered artefacts, then push `main`.

The home page automatically lists public posts, shows categories, and emits an RSS feed (`index.xml`).

## GitHub Pages

The included workflow installs the locked Python environment, renders Quarto to `_site/`, and deploys that directory. After creating `ClausonGeomet/writing`, set **Settings → Pages → Source** to **GitHub Actions**. A project site will then be served at `https://clausongeomet.github.io/writing/`.

### Custom-domain caveat

GitHub Pages cannot by itself map a project site to the custom-domain subpath `clausongeomet.com/writing`. It supports a custom domain at a site root (or a subdomain), not this independent subpath mapping. To serve that URL, integrate the built `_site/` into the main `clausongeomet.com` build, use a Cloudflare Worker/reverse proxy, or publish this site on a subdomain such as `writing.clausongeomet.com`.

## Publishing hygiene

Only commit material that is safe to make public. The `.gitignore` excludes local environments, execution caches, editor state, secrets, and designated private-data directories; it is not a substitute for reviewing every file before publication.
