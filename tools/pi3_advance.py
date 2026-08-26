#!/usr/bin/env python3
"""把束转储中的某个状态按 Π3（自校正补给的纯成长策略）推进到目标时刻，写成 .bin。"""
import struct,sys
H=13
def load(p):
    d=open(p,'rb').read(); cnt,t0=struct.unpack_from('<ii',d,0)
    sz=8*(4+4*H)
    return t0,[list(struct.unpack_from('<%dQ'%(4+4*H),d,8+i*sz)) for i in range(cnt)]
def growact(x,s,w,b,aw,am,asu,ab):
    v=2*ab[1]+ab[2]+ab[3]+2*b-s-asu[1]-asu[2]-asu[3]
    need=(v+7)//8 if v>0 else 0
    x1=x
    p=min(need,w,x1//50); x1-=50*p
    k=min(b,s,x1//25);    x1-=25*k
    q=min(w-p,x1//200)
    return k,p,q
def step(x,s,w,b,aw,am,asu,ab,k,p,q):
    aw=aw[:]; am=am[:]; asu=asu[:]; ab=ab[:]
    x-=25*k+50*p+200*q; s-=k; b-=k
    ab[2]+=k; aw[2]+=k
    aw[3]+=p; asu[3]+=8*p
    aw[12]+=q; asu[12]+=10*q; ab[12]+=q
    g=w-p-q; aw[1]+=g; am[1]+=4*g
    nx=x+am[1]; ns=s+asu[1]; nw=aw[1]; nb=b+ab[1]
    for i in range(1,H-1):
        aw[i]=aw[i+1]; am[i]=am[i+1]; asu[i]=asu[i+1]; ab[i]=ab[i+1]
    aw[H-1]=am[H-1]=asu[H-1]=ab[H-1]=0
    aw[0]=am[0]=asu[0]=ab[0]=0
    return nx,ns,nw,nb,aw,am,asu,ab
def main(src,idxs,tgt,out):
    t0,sts=load(src)
    res=[]
    for idx in idxs:
        v=sts[idx]
        x,s,w,b=v[:4]; aw=v[4:17]; am=v[17:30]; asu=v[30:43]; ab=v[43:56]
        for t in range(t0,tgt):
            k,p,q=growact(x,s,w,b,aw,am,asu,ab)
            assert k<=b and k<=s and p+q<=w and 25*k+50*p+200*q<=x
            x,s,w,b,aw,am,asu,ab=step(x,s,w,b,aw,am,asu,ab,k,p,q)
        res.append([x,s,w,b]+aw+am+asu+ab)
    d=struct.pack('<ii',len(res),tgt)
    for r in res: d+=struct.pack('<%dQ'%(4+4*H),*r)
    open(out,'wb').write(d)
    print("写出 %d 个 t=%d 状态 -> %s"%(len(res),tgt,out))
if __name__=='__main__':
    src=sys.argv[1]; tgt=int(sys.argv[2]); out=sys.argv[3]
    idxs=[int(z) for z in sys.argv[4].split(',')] if len(sys.argv)>4 else [0]
    main(src,idxs,tgt,out)
