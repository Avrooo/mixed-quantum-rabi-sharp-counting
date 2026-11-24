"""Post-run gate: every verification log must contain its expected pass marker and no failure marker.
Exit status 1 if anything is missing or failed. Pure standard library."""
import re,sys,os
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),"logs")
def rd(n): return open(os.path.join(D,n+".log")).read()
bad=[]
def need(n,pat,cnt=None,absent=None,exact=False):
    try: t=rd(n)
    except FileNotFoundError: bad.append(n+": log missing"); return
    if "Traceback" in t: bad.append(n+": Traceback")
    m=re.findall(pat,t,flags=re.M)
    if (cnt is None and not m) or (cnt is not None and ((exact and len(m)!=cnt) or (not exact and len(m)<cnt))): bad.append("%s: expected %s x '%s', found %d"%(n,cnt if cnt is not None else ">=1",pat,len(m)))
    for a in (absent or []):
        if re.search(a,t): bad.append("%s: failure marker '%s'"%(n,a))
need("test_exact",r"98/98 pass")
need("topw",r"constant ratio R1_top/\(Delta g\^par F1\): True, R0: True",21,[r"False"],exact=True)
need("degree_law",r"all pass: True",1)
need("verify_prop12",r"^ALL OK",1,[r"FAIL"])
need("checks1",r"PASS",5,[r"FAIL"])
need("checks3",r"status=certified",18,[r"FAIL"],exact=True)
chk=re.findall(r"B_n=(\d+) certified count=(\d+)",rd("checks3")) if os.path.exists(os.path.join(D,"checks3.log")) else []
if any(a!=b for a,b in chk) or len(chk)!=18: bad.append("checks3: certified count != B_n in some case")
need("checks4",r"PASS",1,[r"FAIL"])
need("symcheck",r"True$",4,[r"False"])
need("verify_krawtchouk",r"KRAWTCHOUK IDENTIFICATION VERIFIED",1)
need("audit_divisibility",r"PASSED.*True",1)
need("audit_tp",r"no \(u,v\)>0 solution: True",1)
need("audit_tp_certified",r"CERTIFIED: a1>0, b1>0, R<0 on the whole bracket: True",1)
need("audit_zero_spin",r"ALL PASS",1)
need("audit_limit_bezout",r"ALL PASS",1)
need("audit_parity_degree",r"RESULT: 120/120",1)
need("audit_appC",r"0 of 1140",1,[r": False"])
need("audit_n1_negative",r"0 of 2000",1)
need("exceptional",r"True$",5,[r"False"])
need("degdrop",r"G1 nonzero:True G0 nonzero:True",None,[r"nonzero:False",r"!=0:False"])
need("audit_independent",r"PASS$",7,[r"FAIL"],exact=True)
need("audit_symbolic_n1",r"True")
need("audit_n1_threshold",r"Q\(-1/2\)= -3")
need("audit_monotone",r"distinct real roots of numerator in \(-1/2,0\): 0")
need("audit_tstar",r"rel\|Res\|")
need("validation_table_run",r"1\.000000000000$",7,exact=True)
need("supp_table5",r"^\d,",7)
need("make_figures3",r"Wigner: int=1\.0000",3,exact=True)
for n in ("make_figures","make_figures2","make_workflow"): need(n,r".",None)
if bad:
    print("LOG CHECK FAILED:"); [print(" -",b) for b in bad]; sys.exit(1)
print("LOG CHECK PASSED (%d logs)"%len(os.listdir(D)))
