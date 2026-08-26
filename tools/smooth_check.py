#!/usr/bin/env python3
"""归一化曲线平滑性检查：g[t]=2F[t]/rho^t 应当单调上升且增量平稳。
凹陷点即为兑现噪声（我们自己的次优），列出最差的若干个。"""
import re,sys
rho=1.053027864
p=sys.argv[1] if len(sys.argv)>1 else 'data/table_best.txt'
lo=int(sys.argv[2]) if len(sys.argv)>2 else 200
hi=int(sys.argv[3]) if len(sys.argv)>3 else 744
b=[int(x) for x in re.findall(r'(\d+)ULL', open(p).read())]
g={t: 2*b[t]/rho**t for t in range(60,745)}
bad=[]
for t in range(lo,hi+1):
    # 与两侧线性插值比较
    a=g[t-1]; c=g[t+1] if t+1<=744 else g[t-1]
    ref=max(g[u] for u in range(max(60,t-4),min(744,t+4)+1))
    bad.append(((g[t]-ref)/ref, t))
bad.sort()
print("最差 25 个凹陷 (相对局部最大):")
for r,t in bad[:25]:
    print("  t=%d  %.3e  2F=%d"%(t,r,2*b[t]))
inc=[(g[t]-g[t-1])/g[t-1] for t in range(lo,hi+1)]
print("每步归一化增量: 中位 %.3e, 最小 %.3e, 最大 %.3e"%(sorted(inc)[len(inc)//2],min(inc),max(inc)))
