#!/usr/bin/env python3
"""P17264 打表校验：39 个官方样例 + 禁区上界余量 + 表不变量。

用法:  python3 src/check.py data/table_best.txt
退出码: 0 = 39/39 全对；否则 = 未通过的样例数
"""
import sys

# 官方两组样例（m -> 期望答案，单位秒）
QS = [100, 1000, 10000,
      9, 79, 99, 126, 166, 266, 426, 666, 999, 1899, 3399, 5599, 9999,
      21999, 45999, 99999, 316227, 999999, 2476413, 2745943, 9999999,
      31622776, 99999999, 316227766, 999999999, 3162277660, 9641625025,
      10691286350, 99999999999, 316227766017, 999999999999, 9999999999999,
      48610229060556, 59770531908338, 66277611238091, 99999999999999]
EX = [20, 200, 640,
      0, 10, 20, 30, 40, 70, 110, 160, 200, 300, 410, 530, 640,
      800, 940, 1090, 1320, 1540, 1710, 1740, 1990, 2210, 2430, 2650,
      2880, 3100, 3320, 3330, 3770, 3990, 4210, 4660, 4970, 5000, 5030, 5100]

RHO = 1.053027864          # 严格 von Neumann 上界，见 doc/solution.md 第 5 节


def load(path):
    txt = open(path).read().replace('ULL', '')
    return [int(x) for x in txt.split(',') if x.strip()]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'data/table_best.txt'
    F = load(path)
    n = len(F)

    print("表: %s  长度=%d" % (path, n))
    mono = all(F[i] >= F[i - 1] for i in range(1, n))
    print("单调非降: %s" % ("是" if mono else "*** 否，表已损坏 ***"))
    if not mono:
        return 99

    # ---- 39 样例 ----
    bad = []
    for m, e in zip(QS, EX):
        need = m // 2 + (m & 1)          # ceil(m/2)，矿物恒为偶数故折半存储
        t = 0
        while t < n and F[t] < need:
            t += 1
        got = t * 10
        if got != e:
            bad.append((m, got, e))

    print("\n=== 39 样例 ===")
    for m, got, e in bad:
        t = e // 10
        gap = m - F[t] * 2
        print("  MISS m=%-16d 输出 %-5d 期望 %-5d  2F[%d]=%d 差 %d (%.3e)"
              % (m, got, e, t, F[t] * 2, gap, gap / m))
    print("bad = %d / 39" % len(bad))

    # ---- 禁区：这些点的 2F[t] 必须严格小于阈值，否则会把原本正确的样例改坏 ----
    # 来源：样例 (m, 答案 10t) 蕴含 2F[t-1] < m
    FORBID = [(173, 2745943), (331, 9641625025),
              (496, 48610229060556), (502, 66277611238091)]
    print("\n=== 禁区上界（改进时不得顶破）===")
    print("%6s %20s %20s %14s" % ("t", "2F[t] 当前", "上界 m", "余量(相对)"))
    ok = True
    for t, m in FORBID:
        if t >= n:
            continue
        cur = F[t] * 2
        room = m - cur
        flag = "" if room > 0 else "   *** 已顶破 ***"
        if room <= 0:
            ok = False
        print("%6d %20d %20d %14.3e%s" % (t, cur, m, room / m, flag))

    # ---- 三个目标点还差多少 ----
    print("\n=== 目标点 ===")
    print("%6s %20s %20s %14s" % ("t", "2F[t] 当前", "需要 >= m", "缺口(相对)"))
    for t, m in [(171, 2476413), (333, 10691286350), (500, 59770531908338)]:
        if t >= n:
            continue
        cur = F[t] * 2
        gap = m - cur
        s = "已达标" if gap <= 0 else "%.3e" % (gap / m)
        print("%6d %20d %20d %14s" % (t, cur, m, s))

    if not ok:
        print("\n*** 警告：有禁区被顶破，该表不可用 ***")
        return 99
    return len(bad)


if __name__ == '__main__':
    sys.exit(main())
