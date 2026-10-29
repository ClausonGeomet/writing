# EMC Part 2 publication release QA

This directory preserves the fresh technical check of the **local publication candidate**. It is a repository QA record, outside Quarto's site render; it is not a dataset, a scientific run, external peer review, or a deployment record. The independent final review of the complete release diff remains pending.

## Release boundary

The release branch started at `37fca05bc78aaa98a743645cf1d400cfc0d5bde8` (then both `origin/main` and the existing release branch). The eight publication paths were committed as `c1baad3b3cce1539643b10882421cb7e4348e8ef`:

| Allowed path | Purpose |
| --- | --- |
| `_quarto.yml` | One favicon URL line changed to local `assets/logo-mark-final.png` |
| `assets/logo-mark-final.png` | Required favicon |
| `posts/emc-part2-what-assays-can-tell-us/index.qmd` | Approved article, with its sole `draft: true` changed to `draft: false` |
| `posts/emc-part2-what-assays-can-tell-us/publication.css` | Article styling |
| `posts/emc-part2-what-assays-can-tell-us/figures/ambiguous_swap.png` | Synthetic ambiguity figure |
| `posts/emc-part2-what-assays-can-tell-us/figures/direct_vs_transfer.png` | Synthetic measurement comparison figure |
| `posts/emc-part2-what-assays-can-tell-us/supplements/evidence-notes.html` | Reader facing evidence notes and HTML download |
| `posts/emc-part2-what-assays-can-tell-us/supplements/measurement-proposal.html` | Unexecuted measurement proposal and HTML download |

The six post files and logo came from approved source commit `1f9a5345a4f4d4416cc62e6f460071f63b74dc0a`. All support files are byte identical. Restoring the article's one draft flag makes it byte identical to source (source SHA-256 `fbdbb02076f71bf90745bd6477a79d1d5c11997a75b25ccfa142ddb58176ee57`; release SHA-256 `09867a08dda5b997c6d2c71ba04bd32f210af7bafd088557be20637f477816f6`). The root config differs from base only in the favicon line. [Source hashes](release-evidence/source-hashes.json) and the [fresh audit](BUILD_AUDIT.md) support these checks.

## Fresh technical evidence

The [full audit](BUILD_AUDIT.md) records an isolated, clean checkout full production `quarto render --no-execute` that exited 0. Quarto emitted a known nonfatal SCSS diagnostic and warnings for unavailable legacy remote font, polyfill, and logo fetches. This was a successful build **with warnings**. An initial cache database error was resolved by using a new empty XDG cache, without a source change. The successful [build log](release-evidence/build-retry.log), initial [cache failure log](release-evidence/build.log), and [draft guard build log](release-evidence/guard-build.log) are retained.

After adding this QA folder, a separate clean isolated full production `--no-execute` packaging preflight also exited 0: 34 site files, `qa_paths_in_site=[]`, all six candidate output hashes identical to the tested release, and no change in the 34-file global site inventory. This preflight is documented in the run-local `qa-curation.md`; it did not repeat the browser check.

The fresh browser check used Chromium at 320, 390, 768, and 1440 CSS pixels: 11 rendered math items at each width, zero root horizontal excess, working desktop outline and two HTML downloads, and zero page errors. The four full article screenshots are [320](release-evidence/article-320.png), [390](release-evidence/article-390.png), [768](release-evidence/article-768.png), and [1440](release-evidence/article-1440.png); [focused 320 pixel math](release-evidence/math-320.png) documents mobile legibility. See [browser metrics](release-evidence/browser-results.json) and [browser log](release-evidence/browser.log).

The [resource inventory](release-evidence/resource-inventory.json) and [resource check log](release-evidence/resource-check.log) record 34/34 site files served with correct HTTP status, MIME, and SHA-256, 796 local references with zero missing paths or bad fragments, and the article CSS, figures, and two named supplement downloads. A separate draft true guard produced only a 90 byte article stub, zero supplement files, and zero slug hits in home, feed, search, sitemap, and listings. No broad resource glob or QA output was deployed. The article's full GitHub source permalinks remain unchanged; the inventory checked their commit and path objects. The [verification script](release-evidence/verify_release.py) and [browser script](release-evidence/browser_check.py) are retained as reproducibility records.

The article interprets frozen **synthetic and exploratory** studies, with internal review qualifications including a retained initial STOP and retrospective inventory issue. The linked measurement study is **proposed and unexecuted**. No scientific code or experiment ran in this release QA. This technical PASS is not a new scientific signoff. Final release scope review still follows; no merge or deployment has occurred.

[Original evidence manifest](release-evidence/manifest.json) indexes the source proof set. [Copy manifest](COPY_MANIFEST.json) records every packaged source and destination, size, SHA-256, and byte equality. The logs, JSON, scripts, and five PNGs are QA files only; no site build tree, browser profile, secret, user data, raw dataset, archive, or original draft copy is included.
