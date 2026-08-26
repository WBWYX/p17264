#!/bin/bash
# 用 5 个被样例钉死的时刻检验一条候选轨迹是否与官方完全相容
# 用法: traj_check.sh <动作文件|-> [状态下标]
A="$1"; SI="${2:-0}"
cd /home/user/p17264
if [ "$A" = "-" ]; then EXTRA=""; else EXTRA="$A $SI"; fi
for seg in "331 333" "496 502"; do
  set -- $seg
  OMP_NUM_THREADS=4 ./build/tablelp data/beam105.bin 94 $1 $2 56 16 1 $EXTRA 2>/dev/null
done | python3 -c "
import sys,re
V={}
for l in sys.stdin:
    m=re.match(r'^(\d+) (\d+)\$',l.strip())
    if m: V[int(m.group(1))]=2*int(m.group(2))
S=[(331,'<=',9641625024),(333,'>=',10691286350),(496,'<=',48610229060554),
   (500,'>=',59770531908338),(502,'<=',66277611238090)]
allok=True
for t,op,b in S:
    if t not in V: continue
    ok=(V[t]<=b) if op=='<=' else (V[t]>=b)
    allok = allok and ok
    print('  t=%-4d 2F=%-18d %s %-18d %s  (%+.4e)'%(t,V[t],op,b,'ok' if ok else '不合格',V[t]/b-1))
print('  === 五点全部相容 ===' if allok else '  === 存在不合格 ===')
"
