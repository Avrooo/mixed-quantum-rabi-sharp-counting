import re,sys,os
REV='/home/claude/pra/rev'
def rd(f): return open(os.path.join(REV,f)).read()
def expand(s):
    return re.sub(r'\\input\{([^}]+)\}',lambda m: expand(rd(m.group(1)+('' if m.group(1).endswith('.tex') else '.tex'))),s)
def roman(n):
    r=''
    for v,c in [(1000,'M'),(900,'CM'),(500,'D'),(400,'CD'),(100,'C'),(90,'XC'),(50,'L'),(40,'XL'),(10,'X'),(9,'IX'),(5,'V'),(4,'IV'),(1,'I')]:
        while n>=v: r+=c; n-=v
    return r
def match_brace(s,i):  # s[i]=='{'
    d=0
    for j in range(i,len(s)):
        if s[j]=='{': d+=1
        elif s[j]=='}':
            d-=1
            if d==0: return j
    raise ValueError
def strip_cmd_arg(s,cmd):  # \cmd{X} -> X
    while True:
        i=s.find(cmd+'{')
        if i<0: return s
        j=match_brace(s,i+len(cmd)); s=s[:i]+s[i+len(cmd)+1:j]+s[j+1:]
def split_rows(b):
    rows=[];cur='';d=0;i=0
    while i<len(b):
        if b.startswith('\\begin{',i): d+=1
        if b.startswith('\\end{',i): d-=1
        if d==0 and b.startswith('\\\\',i): rows.append(cur);cur='';i+=2;continue
        cur+=b[i];i+=1
    rows.append(cur); return rows

body=expand(''.join(rd(f) for f in ['sec01_intro.tex','sec02_model.tex','sec03_termination.tex','sec04_main.tex','sec05_mechanism.tex','sec06_sharp.tex','sec07_n1.tex','sec08_negative.tex','sec09_validation.tex','sec10_discussion.tex','app_all.tex']))
abstract=rd('abstract_text.tex')
labels={}  # label -> display
cnt=dict(sec=0,sub=0,thm=0,prob=0,eq=0,fig=0,tab=0)
appendix=[False]
# --- 0. lists to paragraphs
def lists(s):
    def conv(m):
        kind=m.group(1); b=m.group(2); n=[0]; out=[]
        for it in re.split(r'\\item',b)[1:]:
            it=it.strip(); lab=None
            mm=re.match(r'\[((?:[^\[\]{}]|\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\})*)\]\s*(.*)',it,re.S)
            if mm: lab,it=mm.group(1),mm.group(2)
            else:
                n[0]+=1; lab={'itemize':'\\textbullet','enumerate':'(%d)'%n[0]}[kind]
            out.append('\n\n'+lab+'\\ '+it.strip()+'\n\n')
        return ''.join(out)
    return re.sub(r'\\begin\{(itemize|enumerate|description)\}(.*?)\\end\{\1\}',conv,s,flags=re.S)
body=lists(body)
body=re.sub(r'\\begin\{acknowledgments\}(.*?)\\end\{acknowledgments\}',lambda m:'\n\n\\section*{Acknowledgments}\n\n'+m.group(1).strip()+'\n\n',body,flags=re.S)
# --- 1. sections & labels
def sec(m):
    t=m.group(1); lab=m.group(2)
    if appendix[0]:
        cnt['sec']+=1; L=chr(64+cnt['sec']); disp=L; head='APPENDIX %s: %s'%(L,t)
    else:
        cnt['sec']+=1; disp=roman(cnt['sec']); head='%s. %s'%(disp,t)
    if lab: labels[lab]=disp
    cnt['sub']=0
    return '\n\n\\section*{%s}\n\n'%head
body=body.replace('\\appendix','@@APPENDIX@@')
parts=body.split('@@APPENDIX@@')
out=[]
SUBRE=r'\\subsection\{((?:[^{}]|\{[^{}]*\})*)\}'
SECRE=r'\\section\{((?:[^{}]|\{[^{}]*\})*)\}(?:\\label\{([^}]*)\})?'
def both(m):
    if m.group(0).startswith('\\subsection'):
        cnt['sub']+=1; L=chr(64+cnt['sub'])
        return '\n\n\\subsection*{%s. %s}\n\n'%(L,m.group(3))
    return sec(m)
for k,p in enumerate(parts):
    if k==1: appendix[0]=True; cnt['sec']=0
    if k==0:
        i0=p.index('\\section{Introduction}'); i1=p.index('\\section{Model')
        seg=re.sub(r'\\subsection\{[^}]*\}\s*','\n\n',p[i0:i1]); p=p[:i0]+seg+p[i1:]
    out.append(re.sub(SECRE+'|'+SUBRE,both,p))
body=''.join(out)
# --- 2. theorem-like envs
NAMES={'theorem':'Theorem','proposition':'Proposition','corollary':'Corollary','lemma':'Lemma','definition':'Definition','observation':'Observation','problem':'Open problem','remark':'Remark'}
def thm(m):
    env=m.group(1); opt=m.group(2); b=m.group(3)
    if env=='remark': head='\\textbf{Remark.}'
    else:
        key='prob' if env=='problem' else 'thm'
        cnt[key]+=1; n=cnt[key]
        lb=re.search(r'\\label\{([^}]*)\}',b)
        if lb: labels[lb.group(1)]=str(n); b=b.replace(lb.group(0),'')
        head='\\textbf{%s %d%s.}'%(NAMES[env],n,(' ('+opt+')') if opt else '')
    return '\n\n'+head+' '+b.strip()+'\n\n'
body=re.sub(r'\\begin\{(theorem|proposition|corollary|lemma|definition|observation|problem|remark)\}(?:\[((?:[^\[\]{}]|\{[^{}]*\})*)\])?(.*?)\\end\{\1\}',thm,body,flags=re.S)
body=re.sub(r'\\begin\{proof\}(.*?)\\end\{proof\}',lambda m:'\n\n\\textit{Proof.} '+m.group(1).strip()+'\n\n',body,flags=re.S)
body=body.replace('\\hfill$\\square$','\\quad$\\square$').replace('\\hfill\\ensuremath{\\square}','\\quad$\\square$')
# --- 3. equations
def clean_eq(b):
    b=re.sub(r'\\(label|tag)\{[^}]*\}','',b); b=re.sub(r'\\(nonumber|notag)','',b)
    b=strip_cmd_arg(b,'\\boxed').replace('\\begin{split}','\\begin{aligned}').replace('\\end{split}','\\end{aligned}')
    return b.strip().rstrip(',.;') if False else b.strip()
def numbered(b,force=None):
    lbs=re.findall(r'\\label\{([^}]*)\}',b); tg=re.search(r'\\tag\{([^}]*)\}',b)
    nonum=re.search(r'\\(nonumber|notag)',b) is not None
    if tg: disp=tg.group(1)
    elif nonum and not lbs: disp=None
    else:
        cnt['eq']+=1; disp=str(cnt['eq'])
    for l in lbs: labels[l]=disp
    if tg and not lbs: pass
    return disp
def eqblock(b):
    disp=numbered(b); cb=clean_eq(b)
    s='\n\n\\[%s\\]\n\n'%cb
    if disp: s+='@@EQ(%s)@@\n\n'%disp
    return s
def eq_env(m): return eqblock(m.group(1))
def align_env(m):
    rows=split_rows(m.group(1)); out=''
    for r in rows:
        if not r.strip(): continue
        out+=eqblock(r.replace('&',' '))
    return out
def single(m):
    env,b=m.group(1),m.group(2)
    if env=='equation': return eqblock(b)
    rows=split_rows(b); out=''
    for r in rows:
        if r.strip(): out+=eqblock(r.replace('&',' '))
    return out
body=re.sub(r'\\begin\{(equation|align|gather)\}(.*?)\\end\{\1\}',single,body,flags=re.S)
body=re.sub(r'\\begin\{(align|gather)\*\}(.*?)\\end\{\1\*\}',lambda m:'\n\n\\[\\begin{aligned}%s\\end{aligned}\\]\n\n'%m.group(2),body,flags=re.S)
# --- 4. figures & tables
def fig(m):
    cnt['fig']+=1; n=cnt['fig']; b=m.group(2); wide=m.group(1)=='*'
    g=re.search(r'\\includegraphics\[[^\]]*\]\{([^}]*)\}',b).group(1)
    ci=b.index('\\caption{'); cj=match_brace(b,ci+8); capt=b[ci+9:cj]
    lb=re.search(r'\\label\{([^}]*)\}',b)
    if lb: labels[lb.group(1)]=str(n)
    w='6.5in' if wide else '3.37in'
    r='\n\n@@IMG@@\\includegraphics[width=%s]{%s}\n\nFIG. %d. %s\n\n'%(w,g,n,capt)
    return '\n\n@@WS@@\n\n'+r+'\n\n@@WE@@\n\n' if wide else r
def tab(m):
    cnt['tab']+=1; n=roman(cnt['tab']); b=m.group(2)
    lb=re.search(r'\\label\{([^}]*)\}',b)
    if lb: labels[lb.group(1)]=n
    ci=b.index('\\caption{'); cj=match_brace(b,ci+8); capt=b[ci+9:cj]
    rest=b[:ci]+b[cj+1:]; rest=re.sub(r'\\label\{[^}]*\}','',rest).replace('\\centering','').replace('\\small','')
    r='\n\nTABLE %s. %s\n\n%s\n\n'%(n,capt,rest.strip())
    return '\n\n@@WS@@\n\n'+r+'\n\n@@WE@@\n\n' if m.group(1)=='*' else r
body=re.sub(r'\\begin\{figure(\*?)\}(?:\[[^\]]*\])?(.*?)\\end\{figure\1\}',fig,body,flags=re.S)
body=re.sub(r'\\begin\{table(\*?)\}(?:\[[^\]]*\])?(.*?)\\end\{table\1\}',tab,body,flags=re.S)
# --- 5. paragraph run-ins
body=re.sub(r'\\paragraph\*?\{((?:[^{}]|\{[^{}]*\})*)\}',r'\\textit{\1}',body)
# --- 6. refs / cites
def lab(k):
    if k not in labels: print('MISSING LABEL',k,file=sys.stderr); return '??'
    return labels[k]
body=re.sub(r'\\eqref\{([^}]*)\}',lambda m:'('+lab(m.group(1))+')',body)
body=re.sub(r'\\ref\{([^}]*)\}',lambda m:lab(m.group(1)),body)
# citations: order of first appearance
bib=rd('refs.tex'); items=dict(re.findall(r'\\bibitem\{([^}]*)\}\s*(.*?)(?=\\bibitem|\\end\{thebibliography\})',bib,re.S))
order=[]
def cit(m):
    ks=[k.strip() for k in m.group(1).split(',')]
    for k in ks:
        if k not in items: print('MISSING BIB',k,file=sys.stderr)
        if k not in order: order.append(k)
    return '@@CITE:%s@@'%','.join(ks)
body=re.sub(r'\\cite\{([^}]*)\}',cit,body)
def fmt(m):
    ns=sorted(order.index(k)+1 for k in m.group(1).split(','))
    parts=[];i=0
    while i<len(ns):
        j=i
        while j+1<len(ns) and ns[j+1]==ns[j]+1: j+=1
        parts.append(str(ns[i]) if j==i else ('%d,%d'%(ns[i],ns[j]) if j==i+1 else '%d--%d'%(ns[i],ns[j]))); i=j+1
    return '['+','.join(parts)+']'
body=re.sub(r'@@CITE:([^@]*)@@',fmt,body)
for k in items:
    if k not in order: print('UNCITED',k,file=sys.stderr); order.append(k)
refs=''.join('\n\n[%d] %s\n\n'%(i+1,items[k].strip()) for i,k in enumerate(order))
macros=''.join(l+'\n' for l in rd('macros.tex').split('\n') if l.startswith('\\newcommand'))
open('body.tex','w').write(macros+body+'\n\n\\section*{References}\n\n@@REFS@@\n\n'+refs)
open('abstract.tex','w').write(macros+'\n'+abstract)
print(cnt,len(order),file=sys.stderr)
