from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urljoin,urlparse,unquote
from urllib.request import urlopen
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from threading import Thread
from functools import partial
import json,hashlib,mimetypes,re,subprocess
R=Path('/tmp/emc-release-test-c1baad3/_site'); G=Path('/tmp/emc-release-guard-c1baad3/_site'); E=Path(__file__).parent
slug='emc-part2-what-assays-can-tell-us'
class P(HTMLParser):
 def __init__(self): super().__init__();self.refs=[];self.ids=[]
 def handle_starttag(self,tag,attrs):
  d=dict(attrs)
  if d.get('id'):self.ids.append(d['id'])
  for k in ('href','src','poster','data-src'):
   if d.get(k):self.refs.append((tag,k,d[k]))
  for item in d.get('srcset','').split(','):
   if item.strip():self.refs.append((tag,'srcset',item.strip().split()[0]))
class Quiet(SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(R)));Thread(target=server.serve_forever,daemon=True).start();port=server.server_address[1]
files=sorted(p for p in R.rglob('*') if p.is_file()); rels={str(p.relative_to(R)):p for p in files};parsed={}
for rel,p in rels.items():
 if rel.endswith('.html'):
  h=P();h.feed(p.read_text(errors='replace'));parsed[rel]=h
missing=[];bad_frag=[];local_refs=0;ext=[]
for rel,h in parsed.items():
 for tag,key,raw in h.refs:
  if raw.startswith(('mailto:','tel:','javascript:','data:')):continue
  u=urlparse(urljoin('https://local.invalid/writing/'+rel,raw))
  if u.netloc!='local.invalid':ext.append(raw);continue
  dest=unquote(u.path)
  if not dest.startswith('/writing/'):
   missing.append([rel,raw,'outside mount']);continue
  dest=dest[len('/writing/'):]
  if dest.endswith('/'):dest+='index.html'
  local_refs+=1
  if dest not in rels:missing.append([rel,raw,dest]);continue
  if u.fragment and dest in parsed and unquote(u.fragment) not in parsed[dest].ids:bad_frag.append([rel,raw,dest])
records=[]
for rel,p in rels.items():
 with urlopen(f'http://127.0.0.1:{port}/{rel}') as response:
  body=response.read();mime=response.headers.get_content_type();status=response.status
 expected=mimetypes.guess_type(rel)[0] or 'application/octet-stream'
 records.append({'path':rel,'bytes':len(body),'sha256':hashlib.sha256(body).hexdigest(),'status':status,'mime':mime,'expected_mime':expected,'matches_disk':body==p.read_bytes(),'mime_ok':mime==expected})
index_surfaces=['index.html','index.xml','search.json','sitemap.xml','listings.json']
def counts(root):
 return {x: (root/x).read_text(errors='replace').count(slug) if (root/x).exists() else None for x in index_surfaces}
gfiles=[str(p.relative_to(G)) for p in G.rglob('*') if p.is_file()]
guard={'guard_files':len(gfiles),'candidate_files':[x for x in gfiles if x.startswith('posts/'+slug+'/')], 'guard_stub_bytes':len((G/'posts'/slug/'index.html').read_bytes()) if (G/'posts'/slug/'index.html').exists() else None,'guard_index_counts':counts(G),'release_index_counts':counts(R),'guard_no_supplements':not any('/supplements/' in x and slug in x for x in gfiles),'no_global_resources_glob':not any(x.strip() in ('- posts/**', '- research/**', '- drafts/**') for x in Path('/tmp/emc-release-test-c1baad3/_quarto.yml').read_text().splitlines())}
source=Path('/tmp/emc-release-test-c1baad3/posts')/slug/'index.qmd';raw=source.read_text()
links=re.findall(r'https://github\.com/ClausonGeomet/writing/blob/([0-9a-f]{40})/([^\s)#]+)',raw)
checks=[]
for sha,path in links:
 ok=subprocess.run(['git','-C','/tmp/emc-publication-release','cat-file','-e',f'{sha}:{path}'],capture_output=True).returncode==0
 checks.append({'commit':sha,'path':path,'valid_local_git_object':ok})
local_research_links=re.findall(r'(?<!https://)[(\"\'](?:\.\./)*research/',raw)
summary={'file_count':len(files),'html_count':len(parsed),'local_refs':local_refs,'missing':missing,'bad_fragments':bad_frag,'external_refs':len(ext),'http_all_ok':all(x['status']==200 and x['matches_disk'] and x['mime_ok'] for x in records),'candidate_files':[x for x in rels if x.startswith('posts/'+slug+'/')],'forbidden_output_paths':[x for x in rels if re.search(r'(research/|dataset|archive|\.qmd$|\.py$|qa/|draft)',x,re.I)],'source_permalink_checks':checks,'local_research_relative_links':local_research_links,'guard':guard,'http_port':port}
(E/'resource-inventory.json').write_text(json.dumps({'summary':summary,'files':records},indent=2)+'\n')
print(json.dumps(summary,indent=2))
assert len(files)==34 and len(parsed)==6 and not missing and not bad_frag
assert summary['http_all_ok'] and not summary['forbidden_output_paths'] and not local_research_links
assert all(x['valid_local_git_object'] for x in checks)
assert guard['guard_no_supplements'] and guard['guard_stub_bytes']==90 and all(v==0 for v in guard['guard_index_counts'].values())
assert all(v>0 for v in guard['release_index_counts'].values()) and guard['no_global_resources_glob']
server.shutdown()
