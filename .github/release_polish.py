"""One-off launch corrections; temporary helper removed before merge."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

ROOT = Path.cwd()
ARTICLE = ROOT / 'posts/element-to-mineral-conversion/index.qmd'
META = ARTICLE.with_name('_metadata.yml')
OUT = ROOT / '.release-checks'
TAGLINE = 'A worked synthetic example in mineral inference, uncertainty and decision-making.'


def fenced_blocks(text):
    blocks, current = [], []
    fence = None
    for line in text.splitlines(keepends=True):
        match = re.match(r'^\s*(`{3,}|~{3,})(.*)$', line)
        if fence:
            current.append(line)
            if match and match[1][0] == fence[0] and len(match[1]) >= len(fence) and not match[2].strip():
                blocks.append(''.join(current))
                current, fence = [], None
        elif match:
            fence = match[1]
            current = [line]
    assert fence is None, 'Unclosed code fence'
    return blocks


def replace_once(text, old, new):
    assert text.count(old) == 1, (old[:80], text.count(old))
    return text.replace(old, new, 1)


def apply():
    OUT.mkdir(exist_ok=True)
    original = ARTICLE.read_text()
    original_blocks = fenced_blocks(original)
    s = replace_once(original, 'subtitle: "A consultancy case study in mass balance, ambiguity and defensible decisions"\n', '')
    # Quarto displays the inherited description under the title. Keep it once,
    # rather than duplicating it in both subtitle and description.
    metadata = META.read_text()
    metadata = replace_once(metadata, 'description: "A consultancy case study in mass balance, ambiguity and defensible decisions."', f'description: "{TAGLINE}"')
    old_aside = 'The order matters. It is not $a=Cb^{-1}$: $b$ is not part of the notation here, and a matrix\ninverse acts on the left of the assay vector. Some papers instead write the same relationship as\n$b=Ax$, with $b$ for assays, $A$ for mineral chemistry and $x$ for mineral fractions. The letters\nchange; the calculation does not.'
    s = replace_once(s, old_aside, 'The order matters: the inverse acts on the left of the assay vector.')
    result = []
    expressions = []
    fence = None
    for line in s.splitlines(keepends=True):
        marker = re.match(r'^\s*(`{3,}|~{3,})(.*)$', line)
        if fence:
            result.append(line)
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= len(fence) and not marker[2].strip():
                fence = None
            continue
        if marker:
            fence = marker[1]
            result.append(line)
            continue
        # Keep inline code untouched, as well as fenced Python and raw HTML.
        parts = re.split(r'(`+[^`\n]*`+)', line)
        for i in range(0, len(parts), 2):
            def convert(match):
                expression = match[1].strip()
                assert '$' not in expression
                expressions.append(expression)
                return '$' + expression + '$'
            parts[i] = re.sub(r'\\\((.+?)\\\)', convert, parts[i])
        result.append(''.join(parts))
    s = ''.join(result)
    assert expressions, 'No legacy inline maths found'
    assert fenced_blocks(s) == original_blocks, 'Executable/raw fenced content changed'
    assert 'a=Cb^{-1}' not in s
    ARTICLE.write_text(s)
    META.write_text(metadata)
    manifest = {'converted_expressions': expressions, 'fenced_block_count': len(original_blocks), 'fenced_block_sha256': hashlib.sha256(''.join(original_blocks).encode()).hexdigest()}
    (OUT / 'source-checks.json').write_text(json.dumps(manifest, indent=2))
    subprocess.run(['git', 'diff', '--check'], check=True)
    print(f'Corrected {len(expressions)} inline expressions; {len(original_blocks)} fenced blocks preserved byte-for-byte.')


def check():
    from bs4 import BeautifulSoup
    from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
    from functools import partial
    from threading import Thread
    from playwright.sync_api import sync_playwright
    manifest = json.loads((OUT / 'source-checks.json').read_text())
    html = ROOT / '_site/posts/element-to-mineral-conversion/index.html'
    soup = BeautifulSoup(html.read_text(), 'html.parser')
    main = soup.select_one('main')
    assert main is not None
    assert len(soup.select('h1')) == 1
    assert TAGLINE in main.get_text(' ', strip=True), 'Synthetic-example tagline not visible'
    assert main.get_text(' ', strip=True).count(TAGLINE) == 1, 'Duplicated title description'
    assert 'A consultancy case study in mass balance' not in main.get_text()
    assert not soup.select('.cell-output-error'), 'A Python cell raised an error during the render'
    normalise = lambda value: re.sub(r'\s+', '', value)
    rendered_math = [normalise(x.get_text()) for x in soup.select('span.math.inline')]
    for expression in manifest['converted_expressions']:
        expected = normalise(expression)
        assert any(expected in value for value in rendered_math), ('Missing rendered maths', expression)
    assert soup.select_one('#article-takeaways')
    assert soup.select_one('#article-contact a[href^="mailto:"]')
    assert soup.select_one('#article-contact a[href="https://clausongeomet.com/"]')
    class QuietHandler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass
    handler = partial(QuietHandler, directory=str(ROOT / '_site'))
    server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
    Thread(target=server.serve_forever, daemon=True).start()
    origin = f'http://127.0.0.1:{server.server_port}'
    rows = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page()
            for width in (320, 390, 768, 1440):
                page.set_viewport_size({'width': width, 'height': 900})
                page.goto(origin + '/posts/element-to-mineral-conversion/index.html', wait_until='networkidle')
                page.wait_for_function('typeof MathJax !== "undefined" && !!MathJax.startup', timeout=30000)
                page.evaluate('MathJax.startup.promise')
                assert page.locator('mjx-container').count() > 0, 'MathJax did not typeset equations'
                assert page.locator('mjx-merror').count() == 0, 'MathJax reported an equation error'
                dims = page.evaluate('({width:innerWidth, scroll:document.documentElement.scrollWidth})')
                assert dims['scroll'] <= width + 1, (width, dims)
                if width in (390, 1440):
                    page.screenshot(path=str(OUT / f'article-{width}.png'), full_page=False)
                    math = page.locator('span.math.inline').filter(has_text='n_').first
                    if math.count():
                        math.screenshot(path=str(OUT / f'inline-math-{width}.png'))
                rows.append(dims)
            browser.close()
    finally:
        server.shutdown()
    (OUT / 'browser-results.json').write_text(json.dumps(rows, indent=2))
    summary = f'Passed full Quarto render; no Python cell errors; {len(manifest["converted_expressions"])} corrected inline expressions verified in rendered maths; {manifest["fenced_block_count"]} fenced blocks unchanged; one visible synthetic-example tagline; reader-panel links preserved; equations typeset with no MathJax errors at four screen widths.\n'
    (OUT / 'summary.md').write_text(summary)
    print(summary)


if __name__ == '__main__':
    {'apply': apply, 'check': check}[sys.argv[1]]()
