from pathlib import Path
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from threading import Thread
from functools import partial
from playwright.sync_api import sync_playwright
import json,time
R=Path('/tmp/emc-release-test-c1baad3/_site');E=Path(__file__).parent;slug='posts/emc-part2-what-assays-can-tell-us/index.html'
class Quiet(SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
s=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(R)));Thread(target=s.serve_forever,daemon=True).start();base=f'http://127.0.0.1:{s.server_address[1]}'
out={'base':base,'viewports':[],'downloads':[],'browser_errors':[],'analytics_intercepted':0}
with sync_playwright() as p:
 browser=p.chromium.launch(executable_path='/home/matt/.cache/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-linux64/chrome-headless-shell',headless=True,args=['--no-sandbox','--disable-setuid-sandbox','--disable-dev-shm-usage'])
 out['chromium_version']=browser.version;out['playwright_version']='1.63.0';out['flags']=['--no-sandbox','--disable-setuid-sandbox','--disable-dev-shm-usage']
 for w in (320,390,768,1440):
  ctx=browser.new_context(viewport={'width':w,'height':900},device_scale_factor=1,accept_downloads=True)
  page=ctx.new_page()
  page.on('pageerror',lambda e:out['browser_errors'].append(str(e)))
  def route(r):
   if 'cloudflareinsights.com' in r.request.url or '/cdn-cgi/rum' in r.request.url:
    out['analytics_intercepted']+=1;r.abort()
   else:r.continue_()
  page.route('**/*',route)
  response=page.goto(base+'/'+slug,wait_until='domcontentloaded',timeout=30000)
  page.wait_for_timeout(1800)
  metrics=page.evaluate('''() => ({innerWidth:innerWidth,clientWidth:document.documentElement.clientWidth,scrollWidth:document.documentElement.scrollWidth,dpr:devicePixelRatio,mathSource:document.querySelectorAll('.math').length,mathRendered:document.querySelectorAll('mjx-container').length,bodyHeight:document.body.scrollHeight,tableCount:document.querySelectorAll('table').length,tocCount:document.querySelectorAll('#TOC, nav[role="doc-toc"]').length,outlineCount:document.querySelectorAll('main :is(h1,h2,h3)').length, favicon:document.querySelector('link[rel="icon"]')?.getAttribute('href')})''')
  metrics.update({'width':w,'http_status':response.status,'title':page.title()})
  page.screenshot(path=str(E/f'article-{w}.png'),full_page=True)
  if w==320:
   page.locator('.math').first.scroll_into_view_if_needed();page.screenshot(path=str(E/'math-320.png'))
   nav=page.locator('header').first
   if nav.count():nav.screenshot(path=str(E/'nav-320.png'))
  if w==1440:
   metrics['toc_visible']=page.locator('#TOC').first.is_visible() if page.locator('#TOC').count() else False
   metrics['table_horizontal_overflow']=page.evaluate('''() => [...document.querySelectorAll('table')].map(t=>({scroll:t.scrollWidth,client:t.clientWidth,parentScroll:t.parentElement?.scrollWidth,parentClient:t.parentElement?.clientWidth}))''')
   for label in ('evidence-notes','measurement-proposal'):
    link=page.locator(f'a[download][href*="{label}"]').first
    assert link.count()==1
    with page.expect_download() as d:link.click()
    item=d.value;out['downloads'].append({'label':label,'suggested_filename':item.suggested_filename,'failure':item.failure()})
  out['viewports'].append(metrics);ctx.close()
 browser.close()
(E/'browser-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert not out['browser_errors']
assert len(out['viewports'])==4 and len(out['downloads'])==2 and all(x['failure'] is None for x in out['downloads'])
for x in out['viewports']:
 assert x['http_status']==200 and x['innerWidth']==x['width'] and x['clientWidth']==x['width'] and x['scrollWidth']<=x['width']
 assert x['mathSource']==11 and x['mathRendered']==11 and x['outlineCount']>0 and x['tableCount']>0
assert out['viewports'][-1]['toc_visible']
s.shutdown()
