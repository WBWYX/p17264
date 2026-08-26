import struct,sys
H=13
def load(p):
    d=open(p,'rb').read()
    cnt,t0=struct.unpack_from('<ii',d,0)
    off=8; sz=8*(4+4*H)
    out=[]
    for i in range(cnt):
        v=struct.unpack_from('<%dQ'%(4+4*H),d,off+i*sz)
        st={'x':v[0],'s':v[1],'w':v[2],'b':v[3],
            'aw':list(v[4:4+H]),'am':list(v[4+H:4+2*H]),
            'as':list(v[4+2*H:4+3*H]),'ab':list(v[4+3*H:4+4*H])}
        out.append(st)
    return t0,out
def check(st):
    sa=sum(st['as']); sw=sum(st['aw']); sb=sum(st['ab']); sm=sum(st['am'])
    K=st['w']+sw-4          # worker productions started
    B=st['b']+sb-1          # bases started
    num=st['s']+sa+st['w']+sw-10*st['b']-10*sb
    errs=[]
    if K<0: errs.append('K<0')
    if B<0: errs.append('B<0')
    if num<0 or num%8: errs.append('D not integer: %d'%num)
    D=num//8 if num>=0 and num%8==0 else None
    if D is not None:
        g=st['x']+sm-25+25*K+50*D+200*B
        if g<0 or g%4: errs.append('G not integer: %d'%g)
    return errs,K,B,D
for p in sys.argv[1:]:
    t0,sts=load(p)
    bad=0; ex=None
    for i,st in enumerate(sts):
        e,K,B,D=check(st)
        if e:
            bad+=1
            if ex is None: ex=(i,e,K,B,D)
    print('%s t0=%d n=%d  invalid=%d %s'%(p,t0,len(sts),bad,ex if ex else ''))
