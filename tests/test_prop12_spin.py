import verify_prop12
def test_identity_and_sectors():
    ok,out=verify_prop12.run(12)
    assert ok
    for (l,npos,nzero,dims,_) in out:
        assert npos+npos+nzero==l+1 and nzero==(1 if l%2==0 else 0)
        assert dims[1]==npos           # sector without the zero eigenvalue has dimension ceil(l/2)
