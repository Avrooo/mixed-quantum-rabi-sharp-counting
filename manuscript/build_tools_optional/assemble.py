import zipfile,re,copy,shutil,os,sys
from lxml import etree
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
M='http://schemas.openxmlformats.org/officeDocument/2006/math'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
A='http://schemas.openxmlformats.org/drawingml/2006/main'
def q(ns,t): return '{%s}%s'%(ns,t)
def w(t): return q(W,t)
TPL='/mnt/user-data/uploads/d2994624-PhysRev-template_zysvy0.docx'
def load(z,n): return etree.fromstring(zipfile.ZipFile(z).read(n))
body_doc=load('body_raw.docx','word/document.xml'); abs_doc=load('abs_raw.docx','word/document.xml')
brels=etree.fromstring(zipfile.ZipFile('body_raw.docx').read('word/_rels/document.xml.rels'))
imgmap={r.get('Id'):r.get('Target') for r in brels if r.get('Type').endswith('/image')}
TNR=lambda: etree.fromstring('<w:rFonts xmlns:w="%s" w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman" w:eastAsia="Times New Roman"/>'%W)
def sz(v): 
    e=etree.Element(w('sz')); e.set(w('val'),str(v)); return e
def szcs(v):
    e=etree.Element(w('szCs')); e.set(w('val'),str(v)); return e
RPR_ORDER=['rStyle','rFonts','b','bCs','i','iCs','caps','smallCaps','strike','color','sz','szCs','u','vertAlign']
def fix_rpr(rpr,size=20,mono=False):
    for st in rpr.findall(w('rStyle')): rpr.remove(st)
    for t in ('rFonts','sz','szCs'):
        for e in rpr.findall(w(t)): rpr.remove(e)
    f=TNR()
    if mono:
        for k in ('ascii','hAnsi','cs','eastAsia'): f.set(w(k),'Courier New')
    rpr.append(f); rpr.append(sz(size-(2 if mono else 0))); rpr.append(szcs(size-(2 if mono else 0)))
    kids=sorted(list(rpr),key=lambda e:RPR_ORDER.index(etree.QName(e).localname) if etree.QName(e).localname in RPR_ORDER else 99)
    for k in list(rpr): rpr.remove(k)
    for k in kids: rpr.append(k)
def style_runs(p,size=20):
    for r in p.iter(w('r')):
        par=r.getparent()
        if etree.QName(par).namespace==M: continue
        rpr=r.find(w('rPr'))
        mono=False
        if rpr is None:
            rpr=etree.Element(w('rPr')); r.insert(0,rpr)
        rs=rpr.find(w('rStyle'))
        if rs is not None and 'Verbatim' in rs.get(w('val')): mono=True
        fix_rpr(rpr,size,mono)
    for r in p.iter(q(M,'r')):
        rp=r.find(w('rPr'))
        if rp is None:
            rp=etree.Element(w('rPr')); mr=r.find(q(M,'rPr'))
            r.insert(1 if mr is not None else 0,rp)
        for t in ('sz','szCs','rFonts'):
            for e in rp.findall(w(t)): rp.remove(e)
        f=etree.Element(w('rFonts'))
        for k in ('ascii','hAnsi','cs'): f.set(w(k),'Cambria Math')
        rp.append(f); rp.append(sz(size)); rp.append(szcs(size))
def set_ppr(p,jc='both',before=0,after=100,keep_next=False,ind=None,style=None):
    ppr=p.find(w('pPr'))
    if ppr is None: ppr=etree.Element(w('pPr')); p.insert(0,ppr)
    for c in list(ppr): ppr.remove(c)
    if keep_next: ppr.append(etree.Element(w('keepNext')))
    sp=etree.Element(w('spacing')); sp.set(w('before'),str(before)); sp.set(w('after'),str(after)); sp.set(w('line'),'240'); sp.set(w('lineRule'),'auto'); ppr.append(sp)
    if ind:
        i=etree.Element(w('ind'))
        for k,v in ind.items(): i.set(w(k),str(v))
        ppr.append(i)
    j=etree.Element(w('jc')); j.set(w('val'),jc); ppr.append(j)
def text(p): return ''.join(t.text or '' for t in p.iter(w('t')))
def mk(xml): return etree.fromstring(xml.replace('<w:','<w:',1)) 
NS='xmlns:w="%s" xmlns:m="%s"'%(W,M)
def run(txt,b=False,i=False,sup=False,size=20):
    r=etree.Element(w('r')); rp=etree.SubElement(r,w('rPr'))
    if b: etree.SubElement(rp,w('b'))
    if i: etree.SubElement(rp,w('i'))
    if sup: va=etree.SubElement(rp,w('vertAlign')); va.set(w('val'),'superscript')
    fix_rpr(rp,size)
    t=etree.SubElement(r,w('t')); t.text=txt; t.set('{http://www.w3.org/XML/1998/namespace}space','preserve'); return r
def para(runs,jc='both',**kw):
    p=etree.Element(w('p')); set_ppr(p,jc,**kw)
    for r in runs: p.append(r)
    return p
def eq_table(mathp,num):
    tbl=etree.fromstring('<w:tbl %s><w:tblPr><w:tblW w:w="4320" w:type="dxa"/><w:tblLayout w:type="fixed"/><w:tblCellMar><w:left w:w="0" w:type="dxa"/><w:right w:w="0" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid><w:gridCol w:w="400"/><w:gridCol w:w="3520"/><w:gridCol w:w="400"/></w:tblGrid><w:tr><w:trPr><w:cantSplit/></w:trPr><w:tc><w:tcPr><w:tcW w:w="400" w:type="dxa"/></w:tcPr><w:p/></w:tc><w:tc><w:tcPr><w:tcW w:w="3520" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr></w:tc><w:tc><w:tcPr><w:tcW w:w="400" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr></w:tc></w:tr></w:tbl>'%NS)
    tcs=tbl.findall('.//'+w('tc'))
    p0=tcs[0].find(w('p')); set_ppr(p0,'left',0,0)
    mathp.find(w('pPr')) is not None and mathp.remove(mathp.find(w('pPr')))
    set_ppr(mathp,'center',0,0); tcs[1].append(mathp)
    tcs[2].append(para([run('('+num+')')],'right',before=0,after=0))
    return tbl
def sectp(cols):
    c='<w:cols w:space="360" w:num="2"/>' if cols==2 else '<w:cols w:space="720"/>'
    return etree.fromstring('<w:p %s><w:pPr><w:sectPr><w:type w:val="continuous"/><w:pgSz w:w="12240" w:h="15840" w:orient="portrait"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>%s</w:sectPr></w:pPr></w:p>'%(NS,c))
# ---------- body

def schema_fix(DOC):
    for rp in DOC.iter(q(M,'rPr')):
        if (rp.find(q(M,'scr')) is not None or rp.find(q(M,'nor')) is not None) and rp.find(q(M,'sty')) is not None: rp.remove(rp.find(q(M,'sty')))
    for e in list(DOC.iter(w('bookmarkStart')))+list(DOC.iter(w('bookmarkEnd'))): e.getparent().remove(e)
    for rp in DOC.iter(q(M,'rPr')):
        order=['lit','nor','scr','sty','brk','aln']
        ks=sorted(list(rp),key=lambda e:order.index(etree.QName(e).localname) if etree.QName(e).localname in order else 9)
        for k in list(rp): rp.remove(k)
        for k in ks: rp.append(k)
    for dp in DOC.iter(q(M,'dPr')):
        order=['begChr','sepChr','endChr','grow','shp','ctrlPr']
        ks=sorted(list(dp),key=lambda e:order.index(etree.QName(e).localname) if etree.QName(e).localname in order else 9)
        for k in list(dp): dp.remove(k)
        for k in ks: dp.append(k)
    
    for mc in DOC.iter(q(M,'mcPr')):
        ks=sorted(list(mc),key=lambda e:0 if etree.QName(e).localname=='count' else 1)
        for k in list(mc): mc.remove(k)
        for k in ks: mc.append(k)
    for e in list(DOC.iter(q(M,'sepChr'))):
        if e.get(q(M,'val'))=='': e.getparent().remove(e)
    
schema_fix(body_doc); schema_fix(abs_doc)
for e in list(abs_doc.iter(w('bookmarkStart')))+list(abs_doc.iter(w('bookmarkEnd'))): e.getparent().remove(e)
body=body_doc.find(w('body'))
kids=[k for k in body if etree.QName(k).localname!='sectPr']
out=[]; i=0
in_refs=False
while i<len(kids):
    k=kids[i]; i+=1
    if etree.QName(k).localname=='tbl':
        # table: style
        tblPr=k.find(w('tblPr'))
        for c in list(tblPr): tblPr.remove(c)
        tblPr.append(etree.fromstring('<w:tblW %s w:w="5000" w:type="pct"/>'%NS))
        tblPr.append(etree.fromstring('<w:jc %s w:val="center"/>'%NS))
        tblPr.append(etree.fromstring('<w:tblBorders %s><w:top w:val="single" w:sz="8" w:space="0" w:color="000000"/><w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/></w:tblBorders>'%NS))
        tblPr.append(etree.fromstring('<w:tblCellMar %s><w:left w:w="60" w:type="dxa"/><w:right w:w="60" w:type="dxa"/></w:tblCellMar>'%NS))
        tblPr.append(etree.fromstring('<w:tblLook %s w:val="0000"/>'%NS))
        rows=k.findall(w('tr'))
        for ri,tr in enumerate(rows):
            for tc in tr.findall(w('tc')):
                tcPr=tc.find(w('tcPr'))
                if tcPr is None: tcPr=etree.Element(w('tcPr')); tc.insert(0,tcPr)
                if ri==0:
                    tcPr.append(etree.fromstring('<w:tcBorders %s><w:bottom w:val="single" w:sz="4" w:space="0" w:color="000000"/></w:tcBorders>'%NS))
                for p in tc.findall(w('p')):
                    set_ppr(p,'left',0,0,keep_next=(ri<len(rows)-1)); style_runs(p,18)
            trPr=tr.find(w('trPr'))
            if trPr is None: trPr=etree.Element(w('trPr')); tr.insert(0,trPr)
            trPr.append(etree.Element(w('cantSplit')))
        tp=k.find(w('tr')); 
        out.append(k); out.append(para([],'left',after=60)); continue
    if etree.QName(k).localname=='p' and text(k).strip() in ('@@WS@@','@@WE@@'):
        out.append(sectp(2 if text(k).strip()=='@@WS@@' else 1)); continue
    p=k; t=text(p); sty=p.find(w('pPr')+'/'+w('pStyle')); sty=sty.get(w('val')) if sty is not None else ''
    if t.strip()=='@@REFS@@': in_refs=True; continue
    if sty=='Heading1':
        set_ppr(p,'center',240,120,keep_next=True); style_runs(p); 
        for r in p.iter(w('r')):
            rp=r.find(w('rPr')); 
            if rp.find(w('b')) is None: rp.insert(1,etree.Element(w('b')))
            rp.insert(2,etree.Element(w('caps')))
        out.append(p); continue
    if sty=='Heading2':
        set_ppr(p,'center',160,80,keep_next=True); style_runs(p)
        for r in p.iter(w('r')):
            rp=r.find(w('rPr'))
            if rp.find(w('b')) is None: rp.insert(1,etree.Element(w('b')))
        out.append(p); continue
    if p.find('.//'+q(M,'oMathPara')) is not None and not t.strip():
        nxt=kids[i] if i<len(kids) else None
        if nxt is not None and text(nxt).startswith('@@EQ('):
            num=re.match(r'@@EQ\((.*)\)@@',text(nxt)).group(1); i+=1
            style_runs(p); out.append(eq_table(p,num)); out.append(para([],'left',after=0,before=0)); continue
        set_ppr(p,'center',60,60); style_runs(p); out.append(p); continue
    if p.find('.//'+w('drawing')) is not None:
        for tt in p.iter(w('t')):
            if tt.text and '@@IMG@@' in tt.text: tt.text=tt.text.replace('@@IMG@@','')
        set_ppr(p,'center',120,60,keep_next=True); out.append(p); continue
    if t.startswith('FIG. ') or t.startswith('TABLE '):
        set_ppr(p,'both',60 if t.startswith('FIG') else 120,160 if t.startswith('FIG') else 80,keep_next=t.startswith('TABLE')); style_runs(p)
        _pp=p.find(w('pPr')); _kl=etree.Element(w('keepLines')); _kn=_pp.find(w('keepNext')); _pp.insert(1 if _kn is not None else 0,_kl)
        # bold label
        out.append(p); continue
    if in_refs:
        set_ppr(p,'left',0,60,ind={'left':400,'hanging':400}); style_runs(p,20); out.append(p); continue
    # normal body para
    set_ppr(p,'both',0,100,ind={'firstLine':0}); style_runs(p)
    out.append(p)
# ack heading: find paragraph starting "Acknowledg"
# ---------- abstract & front matter
ap=[x for x in abs_doc.find(w('body')) if etree.QName(x).localname=='p']
abp=ap[0]
for extra in ap[1:]:
    for r in list(extra): abp.append(r)
set_ppr(abp,'both',120,160); 
abp.insert(1,run('ABSTRACT',b=True)); abp.insert(2,run('. '))
style_runs(abp)
front=[]
tp=etree.Element(w('p'))
tpr=etree.SubElement(tp,w('pPr')); ps=etree.SubElement(tpr,w('pStyle')); ps.set(w('val'),'Title')
sp=etree.SubElement(tpr,w('spacing')); sp.set(w('before'),'0'); sp.set(w('after'),'120')
j=etree.SubElement(tpr,w('jc')); j.set(w('val'),'center')
tp.append(run('Exact polynomial eigenstates and sharp counting bounds in the biased mixed quantum Rabi model',b=True,size=24))
front.append(tp)
front.append(para([run('Sakib Shahrear Avro'),run('1,*',sup=True),run('†',sup=True)],'center',after=60))
front.append(para([run('1',sup=True),run('[Affiliation to be supplied by the author: department, institution, city, country]',i=True)],'center',after=40))
front.append(para([run('*',sup=True),run('Contact author: [e-mail address to be supplied by the author]')],'center',after=40))
front.append(para([run('†',sup=True),run('Prepared with the assistance of Claude (Anthropic); see the Acknowledgments.')],'center',after=160))
front.append(abp)
_fs=sectp(1); _s=_fs.find('.//'+w('sectPr')); _s.remove(_s.find(w('type')))
_s.insert(0,etree.fromstring('<w:footerReference xmlns:w="%s" xmlns:r="%s" w:type="default" r:id="rId9"/>'%(W,R)))
_s.find(w('cols')).addprevious(etree.fromstring('<w:pgNumType xmlns:w="%s" w:start="1"/>'%W))
front.append(_fs)
# ---------- write package
tmp='pkg'; shutil.rmtree(tmp,ignore_errors=True)
zipfile.ZipFile(TPL).extractall(tmp)
os.remove(tmp+'/word/media/image1.jpg') if os.path.exists(tmp+'/word/media/image1.jpg') else None
zb=zipfile.ZipFile('body_raw.docx')
rels=open(tmp+'/word/_rels/document.xml.rels').read()
rels=re.sub(r'<Relationship Id="rId(10|11|7|8)"[^>]*/>','',rels)
newrels=''
for n,(rid,tgt) in enumerate(sorted(imgmap.items())):
    name='fig%d_%s'%(n+1,os.path.basename(tgt)); open(tmp+'/word/media/'+name,'wb').write(zb.read('word/'+tgt))
    newrels+='<Relationship Id="rIdF%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/%s"/>'%(n+1,name)
    imgmap[rid]='rIdF%d'%(n+1)
rels=rels.replace('</Relationships>',newrels+'</Relationships>'); open(tmp+'/word/_rels/document.xml.rels','w').write(rels)
for e in body_doc.iter(q(A,'blip')):
    e.set(q(R,'embed'),imgmap[e.get(q(R,'embed'))])
# unique docPr ids
n=1
for e in body_doc.iter('{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}docPr'): e.set('id',str(n)); n+=1
tdoc=open(tmp+'/word/document.xml',encoding='utf8').read()
head=tdoc[:tdoc.index('<w:body>')+len('<w:body>')]
def ser(e): return etree.tostring(e,encoding='unicode')
sect='<w:sectPr><w:type w:val="continuous"/><w:pgSz w:w="12240" w:h="15840" w:orient="portrait"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/><w:cols w:space="360" w:num="2"/></w:sectPr>'
newbody=''.join(ser(e) for e in front+out)
open(tmp+'/word/document.xml','w',encoding='utf8').write(head+newbody+sect+'</w:body></w:document>')
ct=open(tmp+'/[Content_Types].xml').read()
if 'Extension="png"' not in ct: ct=ct.replace('<Default ','<Default Extension="png" ContentType="image/png"/><Default ',1)
ct=re.sub(r'<Override PartName="/word/footer2.xml"[^>]*/>','',ct)
open(tmp+'/[Content_Types].xml','w').write(ct)
os.remove(tmp+'/word/footer2.xml')
if os.path.exists(tmp+'/word/fonts'): pass
core=open(tmp+'/docProps/core.xml').read(); print(core[:600])
core=core.replace('<dc:creator>Jacilynn Otero</dc:creator>','<dc:title>Exact polynomial eigenstates and sharp counting bounds in the biased mixed quantum Rabi model</dc:title><dc:creator>Sakib Shahrear Avro</dc:creator>')
open(tmp+'/docProps/core.xml','w').write(core)
