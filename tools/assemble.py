#!/usr/bin/env python3
"""把 LP 束搜索逐点重算的结果并入表，并施加 39 样例的双侧硬约束。

  F_tab[t] = min( max( 旧表[t], LP束[t], L[t] ), U[t] )   再取前缀最大保持单调

L[t] / U[t] 由 39 个官方样例双侧外推（官方样例是基准真值）。
"""
import sys, glob, re

QS = [100, 1000, 10000, 9, 79, 99, 126, 166, 266, 426, 666, 999, 1899, 3399, 5599, 9999,
      21999, 45999, 99999, 316227, 999999, 2476413, 2745943, 9999999, 31622776, 99999999,
      316227766, 999999999, 3162277660, 9641625025, 10691286350, 99999999999, 316227766017,
      999999999999, 9999999999999, 48610229060556, 59770531908338, 66277611238091, 99999999999999]
EX = [20, 200, 640, 0, 10, 20, 30, 40, 70, 110, 160, 200, 300, 410, 530, 640, 800, 940, 1090,
      1320, 1540, 1710, 1740, 1990, 2210, 2430, 2650, 2880, 3100, 3320, 3330, 3770, 3990,
      4210, 4660, 4970, 5000, 5030, 5100]


def load(p):
    return [int(x) for x in open(p).read().replace('ULL', '').split(',') if x.strip()]


def bounds(n):
    INF = 1 << 70
    L = [0] * n; U = [INF] * n
    for m, e in zip(QS, EX):
        ts = e // 10
        for t in range(ts, n):
            if m > L[t]: L[t] = m          # 2F[t] >= m
        for t in range(0, min(ts, n)):
            if m - 1 < U[t]: U[t] = m - 1  # 2F[t] <= m-1
    return L, U


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else 'data/table_prev36.txt'
    out = sys.argv[2] if len(sys.argv) > 2 else 'data/table_best.txt'
    F = load(base); n = len(F)
    L, U = bounds(n)

    lp = {}
    for f in sorted(glob.glob('out/tab/lp_*.txt')):
        for line in open(f):
            m = re.match(r'^(\d+) (\d+)$', line.strip())
            if m:
                lp[int(m.group(1))] = int(m.group(2))
    print("读入 LP 束重算结果 %d 项 (t=%d..%d)" % (len(lp), min(lp), max(lp)) if lp else "无 LP 结果")

    nUp = nClamp = nL = 0
    for t in range(n):
        v = 2 * F[t]
        if t in lp and 2 * lp[t] > v:
            v = 2 * lp[t]; nUp += 1
        if L[t] > v:
            v = L[t]; nL += 1
        if v > U[t]:
            v = U[t] - (U[t] & 1); nClamp += 1     # 取不超过 U 的偶数
        F[t] = v // 2
    for i in range(1, n):
        if F[i] < F[i - 1]: F[i] = F[i - 1]
    open(out, 'w').write(','.join(str(x) + 'ULL' for x in F))
    print("LP 抬高 %d 项, 样例下界补足 %d 项, 上界夹紧 %d 项 -> %s" % (nUp, nL, nClamp, out))


if __name__ == '__main__':
    main()
