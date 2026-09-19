# Draft workspace

Keep unpublished work under `drafts/<slug>/`. These files are explicitly excluded by `_quarto.yml`, so the deployment workflow does not publish them.

To publish a draft, move it to `posts/<slug>/`, set `draft: false` in its `_metadata.yml` or document front matter, add it to version control, and render the site. Public posts are discovered automatically by the home-page listing.
