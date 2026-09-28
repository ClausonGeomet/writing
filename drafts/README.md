# Draft workspace

This repository is public, so **unpublished drafts should not be committed here**.

The repository's `.gitignore` ignores everything under `drafts/` except this README. I normally keep working drafts locally or in a private research/notes repository until they are ready to publish.

When a piece is ready:

1. move the finished article to `posts/<slug>/`;
2. set `draft: false` in its metadata or front matter;
3. add it to version control; and
4. render and check the site before merging.

Quarto also excludes the draft area from the deployed site, but the Git ignore rule is the important protection because the GitHub repository itself is public.
