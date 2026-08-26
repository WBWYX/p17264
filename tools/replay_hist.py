# 独立复算：用未折半的原始单位从 t=0 重放动作序列，逐状态与束转储比对。
import struct,sys
H=13
def loadbin(p):
    d=open(p,'rb').read(); cnt,t0=struct.unpack_from('<ii',d,0)
    sz=8*(4+4*H); out=[]
    for i in range(cnt):
        v=struct.unpack_from('<%dQ'%(4+4*H),d,8+i*sz)
        out.append(v)
    return t0,out
def replay(acts):
    # 原始单位
    M=50; S=6; W=4; B=1
    aw=[0]*H; am=[0]*H; asup=[0]*H; ab=[0]*H     # am 为原始矿物
    for (k,p,q) in acts:
        assert k<=B and k<=S and p+q<=W, "动作违反数量约束"
        cost=50*k+100*p+400*q
        assert cost<=M, "矿物不足"
        M-=cost; S-=k; B-=k
        ab[2]+=k; aw[2]+=k                 # 基地 20s 归位, 新工人 20s 后出现
        aw[3]+=p; asup[3]+=8*p             # 补给站 30s
        aw[12]+=q; asup[12]+=10*q; ab[12]+=q   # 基地 120s
        g=W-p-q; aw[1]+=g; am[1]+=8*g      # 采矿 10s, +8
        M+=am[1]; S+=asup[1]; W=aw[1]; B+=ab[1]
        for i in range(1,H-1):
            aw[i]=aw[i+1]; am[i]=am[i+1]; asup[i]=asup[i+1]; ab[i]=ab[i+1]
        aw[H-1]=am[H-1]=asup[H-1]=ab[H-1]=0
        aw[0]=am[0]=asup[0]=ab[0]=0
        assert M>=0 and S>=0
    return M,S,W,B,aw,am,asup,ab
def main(histp,binp,limit=None):
    t0,sts=loadbin(binp)
    f=open(histp); cnt,tt=map(int,f.readline().split())
    assert tt==t0, "时刻不符 %d vs %d"%(tt,t0)
    n=min(cnt,limit) if limit else cnt
    bad=0
    for i in range(cnt):
        parts=f.readline().split()
        if i>=n: continue
        ln=int(parts[0]); nums=list(map(int,parts[1:]))
        assert ln==tt, "历史长度 %d != %d"%(ln,tt)
        acts=[(nums[3*j],nums[3*j+1],nums[3*j+2]) for j in range(ln)]
        M,S,W,B,aw,am,asup,ab=replay(acts)
        v=sts[i]
        exp=(M//2,S,W,B)+tuple(aw)+tuple(x//2 for x in am)+tuple(asup)+tuple(ab)
        if M%2 or any(x%2 for x in am):
            print("状态 %d: 原始矿物为奇数!"%i); bad+=1; continue
        if exp!=v:
            bad+=1
            if bad<=3:
                print("状态 %d 不符:"%i)
                print("  重放 x=%d s=%d w=%d b=%d"%(M//2,S,W,B))
                print("  转储 x=%d s=%d w=%d b=%d"%(v[0],v[1],v[2],v[3]))
    print("%s: 检查 %d 个状态, 不符 %d"%(binp,n,bad))
if __name__=='__main__':
    main(sys.argv[1],sys.argv[2], int(sys.argv[3]) if len(sys.argv)>3 else None)
