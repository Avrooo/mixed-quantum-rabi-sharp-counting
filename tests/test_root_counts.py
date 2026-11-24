import csv
from fractions import Fraction as F
def test_certified_counts():
    rows=list(csv.DictReader(open("data/CERTIFIED_ROOT_COUNTS.csv")))
    pos=[r for r in rows if F(r["t"])>0]; neg=[r for r in rows if F(r["t"])<0]
    assert len(pos)==187 and sum(int(r["count"])==int(r["bezout"]) for r in pos)==68
    assert all(int(r["count"])<=int(r["bezout"]) for r in pos)
    cert=[r for r in neg if r["status"].startswith("certified")]
    assert len(cert)==111 and all(r["count"]=="0" for r in cert)
    assert [r["status"] for r in neg if not r["status"].startswith("certified")]==["TIMEOUT_240s"]*2   # 2 uncertified rows are not claimed
