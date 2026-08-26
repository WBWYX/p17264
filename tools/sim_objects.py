# 完全独立的模拟器：以秒为时间轴、显式到账队列，不折半、不复用 step() 结构。
import sys
from collections import Counter
def simulate(acts):
    idleW=4; idleB=1; minerals=50; supplies=6
    arrW=Counter(); arrB=Counter(); arrM=Counter(); arrS=Counter()
    for step,(k,p,q) in enumerate(acts):
        T=step*10
        minerals+=arrM.pop(T,0); supplies+=arrS.pop(T,0)
        idleW+=arrW.pop(T,0);    idleB+=arrB.pop(T,0)
        if k>idleB: raise AssertionError("t=%ds 空闲基地不足 %d>%d"%(T,k,idleB))
        if k>supplies: raise AssertionError("t=%ds 补给不足 %d>%d"%(T,k,supplies))
        if p+q>idleW: raise AssertionError("t=%ds 空闲工人不足 %d>%d"%(T,p+q,idleW))
        cost=50*k+100*p+400*q
        if cost>minerals: raise AssertionError("t=%ds 矿物不足 %d>%d"%(T,cost,minerals))
        minerals-=cost; supplies-=k; idleB-=k; idleW-=(p+q)
        arrB[T+20]+=k;  arrW[T+20]+=k
        arrW[T+30]+=p;  arrS[T+30]+=8*p
        arrW[T+120]+=q; arrS[T+120]+=10*q; arrB[T+120]+=q
        g=idleW; idleW=0
        arrW[T+10]+=g;  arrM[T+10]+=8*g
        if minerals<0 or supplies<0: raise AssertionError("负值 t=%ds"%T)
    T=len(acts)*10
    minerals+=arrM.pop(T,0); supplies+=arrS.pop(T,0)
    idleW+=arrW.pop(T,0);    idleB+=arrB.pop(T,0)
    return T,minerals,supplies,idleW,idleB
if __name__=='__main__':
    acts=[tuple(map(int,l.split())) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
    T,M,S,W,B=simulate(acts)
    print("步数 %d, 时间 %d 秒, 矿物 %d, 补给 %d, 空闲工人 %d, 空闲基地 %d"%(len(acts),T,M,S,W,B))
