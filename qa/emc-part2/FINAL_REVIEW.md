# Final publication release review — PASS WITH NOTE

Reviewed 2026-10-29. **PASS WITH NOTE — the publication scope and fresh technical evidence converge. The known nonfatal SCSS diagnostic and four legacy resource-fetch warnings remain disclosed; they do not block this publication-only PR.** No unresolved release blocker was found.

Exact reviewed HEAD: `69ce8694abe7630f1b2f6669c0a6c850509a7c85` on `release/emc-part2-what-assays-can-tell-us`.
Core publication / browser-tested source: `c1baad3b3cce1539643b10882421cb7e4348e8ef`.
Base main: `37fca05bc78aaa98a743645cf1d400cfc0d5bde8`.
Approved research source: `1f9a5345a4f4d4416cc62e6f460071f63b74dc0a`.

## Clearance and boundary

Explicit clearance for the already authorized **normal fast-forward push of the existing release branch only**, followed by a **publication-only PR against main**. No force push, merge, automerge, deployment, workflow dispatch, research-branch push, or new experiment is cleared. No additional user permission question is required for those authorized ship actions. Shipper must preserve the reviewed scope and use the actual final head SHA in its handoff.

Scriber may subsequently copy this review byte-for-byte to `qa/emc-part2/FINAL_REVIEW.md` and update only the QA README status/provenance. This is an administrative addition, not a new browser-tested source. Preserve the exact reviewed SHA above, all existing evidence, and all eight core file bytes. Any other substantive change falls outside this clearance. This reviewer made no target-repository edits and performed no test, render, browser, scientific, or network execution.

## Scope, ancestry, and scientific preservation

Read request, impact, credentials, release preparation, audit, test handoff, QA curation, docs, log entry, architecture, mailbox, reviewer definition and active reviewer profile. Inspected actual Git diff, tree sizes/modes, commit parents, source comparison, working-tree status, evidence scripts/logs/JSON, all four full article screenshots and the focused mobile math screenshot.

The actual base-to-HEAD diff has 27 paths: exactly eight core paths plus 19 QA paths. Core paths are `_quarto.yml`, `assets/logo-mark-final.png`, and the post's `index.qmd`, `publication.css`, two figure PNGs and two HTML supplements. Config changes only the favicon URL to its local asset. No new broad resources glob, root style/filter/deployment change, research tree, original draft, scientific code, model/API, dataset or numerical result is included. The Python files in QA are the small retained verification scripts, not scientific implementations.

Actual Git comparison against approved `1f9a534` shows just `draft: true` to `draft: false` across the six post files and logo; all other bytes match. The source-hash ledger identifies article source SHA-256 `fbdbb02076f71bf90745bd6477a79d1d5c11997a75b25ccfa142ddb58176ee57` and release SHA-256 `09867a08dda5b997c6d2c71ba04bd32f210af7bafd088557be20637f477816f6`. `execute.enabled: false` remains. Synthetic/exploratory limits, internal rather than external review, prior-harm examples, initial STOP history, retrospective inventory qualification, and unexecuted measurement proposal are preserved. This is not a renewed scientific review.

Parents are exactly `69ce869 -> c1baad3 -> 37fca05`: no merge, rebase or cherry-picked research stack. The release worktree is clean. The research worktree remains at `1f9a534`, with only the previously recorded untracked raw-run directory. Historical project reviews are not overwritten. Credentials/preparation record unchanged old PRs 4/5/6 and remote base; no remote mutation was performed by this reviewer. Live remote state was not re-queried here.

The unchanged deployment workflow triggers on main pushes or explicit workflow dispatch, not release pushes or PR creation. Candidate draft false therefore does not deploy through the cleared ship actions; a later main merge remains unauthorized.

## Evidence challenges and convergence

| Challenge | Evidence inspected and disposition |
| --- | --- |
| Fresh full production build | Clean archive and empty XDG cache commands in audit; successful full three-page render log ends with `_site/index.html`, reported exit 0. Initial cache database failure is retained, with retry explained. No scientific execution. |
| Resources and downloads | Inspected verifier assertions and all 34 inventory records: HTTP 200, expected MIME, matching served/disk bytes and SHA records. Six candidate outputs present; supplements are actual `text/html` downloads, not PDFs. |
| References and guard | 796 local references, zero missing paths/fragments. Separate draft-true archive has a 90-byte stub, zero supplements and zero slug hits on all five index surfaces. The script asserts these outcomes without swallowing assertion failures. |
| Immutable source links | Full accepted commit/file objects recorded valid. This establishes local object provenance, not a fresh external HTTP availability claim. |
| Browser credibility | Script uses actual Chromium/Playwright, four contexts at 320/390/768/1440 CSS pixels, real DOM metrics and download events. JSON/log agree: 11 rendered math items per width, no root excess, zero page errors, desktop TOC visible, two successful named HTML downloads. Desktop table metrics show no excess. |
| Visual evidence | Inspected the four full article PNGs for layout/asset continuity and the focused 320-pixel image for readable math. No visible release-blocking defect found. Full-page images are scaled in review; detailed math legibility is supported by the focused image and DOM evidence. |
| QA packaging | Copy manifest records 17 exact source copies totaling 7,843,095 bytes, plus README/manifest metadata. Git tree confirms the small expected files, regular file modes, and no hidden bulk/raw/private payload. All audit-relative evidence targets are present. |
| Packaging output equivalence | QA curation records a clean no-execute archive preflight, 34 outputs, no QA output paths, and zero changed output bytes across the entire site. Inspected its successful build log. QA-only commit diff leaves every tested core byte unchanged; preflight is not presented as a second browser run. |

Pipeline isolation: this explicitly scoped Workflow 7 run has no new algorithm/specification pipelines. Separate preparation, independent technical audit, and documentation packaging converge on request/impact. No applicable isolation breach was found. Comprehension/spec/test-spec/simulation gates are waived by the dispatch for this publication-only gate; their absence is not represented as completion. No research campaign or simulation is required.

Tolerance integrity: the audit table uses exact bytes/counts/zero-error criteria, consistent with the retained assertions and recorded outputs. No scientific atol/rtol comparison applies, and no tolerance inflation, swallowed failure, or omitted required release scenario was identified. The per-test table is present. Before/after scientific metrics are inapplicable; packaging equivalence is the relevant comparison.

Documentation: docs and process log exist, including filename header, results table, before/after applicability, decisions and handoff. Run-local Mermaid architecture reflects the publication flow. A new root architecture file is correctly excluded by the explicit eight-core-plus-QA allowlist; the generic template requirement does not override user scope. No new root workflow-log artifacts or Brain entries are included.

## Remaining notes

- The build is successful **with warnings**, never warning-free: one known nonfatal Quarto SCSS color-variable export diagnostic and four legacy remote-fetch warnings per successful render. The latter concern the existing font, polyfill and external logo. General site maintenance can address these separately; no unrelated cleanup belongs in this release.
- QA README currently says final review pending. Scriber may make the expressly authorized status update and exact review copy above. Retain the distinction between reviewed metadata HEAD and browser-tested core SHA.
- Source supplements retain approved staging/provenance wording. Their bytes and scientific limits are intentionally preserved by the user's scope. Local object-link checks and loopback MIME checks are not claims of live GitHub/CDN availability or hosted-server MIME behavior.

All applicable release challenges cleared. Route to shipper after any routine Scriber status recording; no STOP routing required.
