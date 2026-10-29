# Release verification audit — PASS

**Target:** `/tmp/emc-publication-release`, `release/emc-part2-what-assays-can-tell-us` at `c1baad3b3cce1539643b10882421cb7e4348e8ef`, based on `origin/main` `37fca05bc78aaa98a743645cf1d400cfc0d5bde8`. Working tree clean before and after. This is the Workflow 7 publication release check under `request.md` and `impact.md`; no new algorithm or numerical test specification applies. No source edit, experiment, commit, push, PR, merge or deployment was performed.

## Scope and source identity

`git diff --name-status origin/main...HEAD` showed exactly eight paths: `_quarto.yml`; `assets/logo-mark-final.png`; the post's `index.qmd`, `publication.css`, two figure PNGs and two supplement HTMLs. `_quarto.yml` changes only the favicon from an external URL to `assets/logo-mark-final.png`. No unrelated research, proposal, draft, code, filter, deployment or root style path changed. [`source-hashes.json`](release-evidence/source-hashes.json) records SHA-256 and byte counts for all seven added files: six post files plus favicon. All are byte-identical to approved `1f9a5345a4f4d4416cc62e6f460071f63b74dc0a` except the post's **single** `draft: true` → `draft: false` flag. The source article SHA-256 is `fbdbb02076f71bf90745bd6477a79d1d5c11997a75b25ccfa142ddb58176ee57`; the release variant is `09867a08dda5b997c6d2c71ba04bd32f210af7bafd088557be20637f477816f6`. Historical review of the approved source is inherited as provenance, not claimed as a fresh review of this branch.

## Commands, environment and warnings

- `git archive c1baad3b3cce1539643b10882421cb7e4348e8ef | tar -x -C /tmp/emc-release-test-c1baad3` created the clean isolated production input.
- `(cd /tmp/emc-release-test-c1baad3 && XDG_CACHE_HOME=/tmp/emc-release-cache-c1baad3 quarto render --no-execute)` **exit 0**, created `_site/index.html`; full output: [`build-retry.log`](release-evidence/build-retry.log).
- An initial invocation with `QUARTO_CACHE_DIR=/tmp/emc-release-test-c1baad3/.fresh-quarto-cache` failed before a site was built with `ERROR: unable to open database file`; full output: [`build.log`](release-evidence/build.log). The successful retry used a new, empty XDG cache and the same archive and full production config. This was an environment cache failure, not a source repair.
- A second `git archive` into `/tmp/emc-release-guard-c1baad3`, changing only the candidate flag to `draft: true`, followed by `(cd /tmp/emc-release-guard-c1baad3 && XDG_CACHE_HOME=/tmp/emc-release-guard-cache-c1baad3 quarto render --no-execute)` **exit 0**; full output: [`guard-build.log`](release-evidence/guard-build.log).
- `python3 release-evidence/verify_release.py` and cached Python 3.12.3 Playwright `release-evidence/browser_check.py` each **exit 0** with elevated loopback execution, because the default sandbox returned `EPERM` on localhost socket creation. Scripts and full output: [`verify_release.py`](release-evidence/verify_release.py), [`resource-check.log`](release-evidence/resource-check.log), [`browser_check.py`](release-evidence/browser_check.py), [`browser.log`](release-evidence/browser.log).

Quarto 1.10.18 issued the known nonfatal SCSS parsing/color-variable export diagnostic once in each successful render. Each also emitted four fetch warnings for an existing Google font, polyfill CDN, and existing external logo (twice), with network DNS unavailable in the sandbox. There was no fatal `ERROR:` in either successful log. This is a **clean checkout, successful build with warnings**, not a warning-free build. No global style or external-resource changes were made. Chromium headless shell 153.0.8010.12 and Playwright 1.63.0 ran with `--no-sandbox --disable-setuid-sandbox --disable-dev-shm-usage` against only the trusted loopback preview; analytics beacons were intercepted (4), production source unchanged.

## Test results

| Test / metric | Expected | Actual | Tolerance | Verdict |
| --- | --- | --- | --- | --- |
| Diff scope | Exact eight allowed paths | Exact eight allowed paths | Exact | PASS |
| Approved source bytes | Seven support/article files equal, only candidate draft flag changed | All seven checks true; one flag change | Exact bytes | PASS |
| Production render | Exit 0, full site, no execution | Exit 0, 34 files and 6 HTML pages | Exact | PASS |
| Candidate package | Article, CSS, two PNGs, two HTML downloads | All six present; HTTP 200, SHA and MIME correct | Exact | PASS |
| Recursive site HTTP | All 34 files match disk and MIME | 34/34 | Exact | PASS |
| Local resources and anchors | No missing local path or fragment | 796 local references, 0 missing, 0 bad fragments | Exact | PASS |
| Source permalinks | Referenced commit/path objects exist | All `github.com/ClausonGeomet/writing/blob/...` targets valid in local Git objects | Exact | PASS |
| Draft guard | Only 90 byte stub; no supplement or listing/feed/search/sitemap slug | 90 byte stub, 0 supplements; five surface counts all 0 | Exact | PASS |
| Global resource glob / leaks | No broad `posts/**`, research, dataset, archive, QA or unpublished draft output | No broad glob; no forbidden output path | Exact | PASS |
| Browser widths | 320, 390, 768, 1440; no root excess | `innerWidth = clientWidth = scrollWidth` at each width; DPR 1 | Exact | PASS |
| Rendered math | 11 | 11/11 at every width | Exact | PASS |
| Desktop outline and table | Visible TOC, table present | TOC visible at 1440; one table; no desktop table overflow | Exact | PASS |
| HTML downloads | Two working named downloads | Both succeeded with expected filenames | Exact | PASS |
| Browser runtime errors | Zero | Zero page errors | Exact | PASS |

The full per-file HTTP inventory, path and anchor results, guard counts and permalink checks are in [`resource-inventory.json`](release-evidence/resource-inventory.json). Per-width DOM metrics and download outcomes are in [`browser-results.json`](release-evidence/browser-results.json). Four fresh full article screenshots: [`320`](release-evidence/article-320.png), [`390`](release-evidence/article-390.png), [`768`](release-evidence/article-768.png), [`1440`](release-evidence/article-1440.png); focused mobile math: [`math-320.png`](release-evidence/math-320.png). The 320 math image was inspected visually and is legible. The inherited 67-shot QA belongs to prior approved web assets and was not rerun. No scientific chunk was executed and no old research test spec was invoked. See [`manifest.json`](release-evidence/manifest.json) for evidence sizes and hashes (7,834,015 bytes before manifest).

**Verdict: PASS for publication release verification**, with the existing Quarto SCSS and remote fetch warnings disclosed above. This technical verdict does not authorize merge or deployment.
