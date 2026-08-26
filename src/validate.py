#!/usr/bin/env python3
"""对表做完整的双侧校验：把 39 个官方样例翻译成对每个 t 的上下界，逐点检查。

样例 (m, 10*t*) 蕴含:  2F[t*] >= m  且  2F[t*-1] < m
再由单调性外推到所有 t:
    L[t] = max{ m : 存在样例答案时刻 t* <= t }        (2F[t] 必须 >= 它)
    U[t] = min{ m-1 : 存在样例答案时刻 t* >= t+1 }    (2F[t] 必须 <= 它)

用法: python3 src/validate.py data/table_best.txt
"""
import sys

QS = [100, 1000, 10000,
      9, 79, 99, 126, 166, 266, 426, 666, 999, 1899, 3399, 5599, 9999,
      21999, 45999, 99999, 316227, 999999, 2476413, 2745943, 9999999,
      31622776, 99999999, 316227766, 999999999, 3162277660, 9641625025,
      10691286350, 99999999999, 316227766017, 999999999999, 9999999999999,
      48610229060556, 59770531908338, 66277611238091, 99999999999999]
EX = [20, 200, 640, 0, 10, 20, 30, 40, 70, 110, 160, 200, 300, 410, 530, 640,
      800, 940, 1090, 1320, 1540, 1710, 1740, 1990, 2210, 2430, 2650,
      2880, 3100, 3320, 3330, 3770, 3990, 4210, 4660, 4970, 5000, 5030, 5100]


def load(path):
    return [int(x) for x in open(path).read().replace('ULL', '').split(',') if x.strip()]


def bounds(n):
    """返回 (L, U): 2F[t] 必须落在 [L[t], U[t]]。"""
    INF = 1 << 70
    L = [0] * n
    U = [INF] * n
    for m, e in zip(QS, EX):
        ts = e // 10
        for t in range(ts, n):            # 2F[t] >= 2F[ts] >= m
            if m > L[t]:
                L[t] = m
        for t in range(0, min(ts, n)):    # 2F[t] <= 2F[ts-1] < m
            if m - 1 < U[t]:
                U[t] = m - 1
    return L, U


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'data/table_best.txt'
    F = load(path)
    n = len(F)
    L, U = bounds(n)
    print("表: %s  长度=%d" % (path, n))

    mono = all(F[i] >= F[i - 1] for i in range(1, n))
    print("单调非降: %s" % ("是" if mono else "*** 否 ***"))

    lo_bad = [(t, 2 * F[t], L[t]) for t in range(n) if 2 * F[t] < L[t]]
    hi_bad = [(t, 2 * F[t], U[t]) for t in range(n) if 2 * F[t] > U[t]]
    print("\n=== 下界违反（表偏小，答案会偏大）: %d 处 ===" % len(lo_bad))
    for t, v, b in lo_bad[:12]:
        print("   t=%-5d 2F=%-18d 需 >= %-18d 差 %d (%.3e)" % (t, v, b, b - v, (b - v) / b))
    print("=== 上界违反（表偏大，答案会偏小）: %d 处 ===" % len(hi_bad))
    for t, v, b in hi_bad[:12]:
        print("   t=%-5d 2F=%-18d 需 <= %-18d 超 %d (%.3e)" % (t, v, b, v - b, (v - b) / b))

    # 39 样例结论
    bad = 0
    for m, e in zip(QS, EX):
        need = m // 2 + (m & 1)
        t = 0
        while t < n and F[t] < need:
            t += 1
        if t * 10 != e:
            bad += 1
    print("\n39 样例: bad = %d / 39" % bad)

    # 最紧的几个上界余量
    print("\n=== 最紧的 6 个上界余量 ===")
    room = sorted(((U[t] - 2 * F[t]) / U[t], t) for t in range(n) if U[t] < (1 << 70))
    for r, t in room[:6]:
        print("   t=%-5d 2F=%-18d 上界 %-18d 余量 %.3e" % (t, 2 * F[t], U[t], r))
    return 0 if (mono and not lo_bad and not hi_bad and bad == 0) else 1


if __name__ == '__main__':
    sys.exit(main())
