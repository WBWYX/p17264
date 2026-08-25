#!/usr/bin/env python3
"""逐点取最大合并多张表，写出 merged.txt。

用法:  python3 src/merge.py out.txt in1.txt in2.txt [in3.txt ...]

合法性：表中每一项都是「时刻 t 可达的折半矿物量」的一个下界，恒 <= 真值。
逐点取最大仍是下界，故合并永远安全，只会让答案单调趋近正确。
合并后重做前缀最大以保持非降。
"""
import sys


def load(path):
    txt = open(path).read().replace('ULL', '')
    return [int(x) for x in txt.split(',') if x.strip()]


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    out, ins = sys.argv[1], sys.argv[2:]
    tabs = [load(p) for p in ins]
    n = max(len(t) for t in tabs)
    F = [0] * n
    for t in tabs:
        for i, v in enumerate(t):
            if v > F[i]:
                F[i] = v
    for i in range(1, n):
        if F[i] < F[i - 1]:
            F[i] = F[i - 1]
    open(out, 'w').write(','.join(str(x) + 'ULL' for x in F))
    print("合并 %d 张表 -> %s  (长度 %d)" % (len(tabs), out, n))
    for p, t in zip(ins, tabs):
        diff = sum(1 for i in range(min(len(t), n)) if t[i] < F[i])
        print("  %-40s 有 %d 项被其他表超过" % (p, diff))
    return 0


if __name__ == '__main__':
    sys.exit(main())
