#!/usr/bin/env python3
"""独立校验可达性证书：按题面原始语义（矿物不折半、显式事件表）从初始状态重放动作序列，
逐步断言四条约束，最后核对到达时刻与矿物量。

用法:  python3 src/verify_cert.py data/cert/t500.txt [...]
退出码: 0 = 全部通过
"""
import sys
from collections import defaultdict


class Sim:
    """题面语义：时间单位 10 秒。初始 1 空闲基地 / 4 空闲工人 / 6 补给 / 50 矿物。

    也可用 START 行从一个中间状态（折半口径的束状态）出发，此时把流水线数组
    展开成绝对时刻的事件表。"""

    def __init__(self, start=None):
        self.t = 0
        self.M, self.S, self.W, self.B = 50, 6, 4, 1
        self.ev = defaultdict(lambda: [0, 0, 0, 0])   # 时刻 -> [工人, 补给, 基地, 矿物]
        if start is not None:
            t0, v = start
            self.t = t0
            self.M, self.S, self.W, self.B = 2 * v[0], v[1], v[2], v[3]
            aw, am, as_, ab = v[4:17], v[17:30], v[30:43], v[43:56]
            for i in range(1, 13):
                if aw[i]:  self.ev[t0 + i][0] += aw[i]
                if as_[i]: self.ev[t0 + i][1] += as_[i]
                if ab[i]:  self.ev[t0 + i][2] += ab[i]
                if am[i]:  self.ev[t0 + i][3] += 2 * am[i]

    def act(self, k, p, q):
        # k 个空闲基地造工人 / p 个工人建补给站 / q 个工人建基地 / 其余工人采矿
        assert k >= 0 and p >= 0 and q >= 0, "t=%d 动作为负" % self.t
        assert k <= self.B, "t=%d: 造工人 %d 个但只有 %d 个空闲基地" % (self.t, k, self.B)
        assert k <= self.S, "t=%d: 造工人 %d 个但只有 %d 点补给" % (self.t, k, self.S)
        assert p + q <= self.W, "t=%d: 派出 %d 个工人但只有 %d 个空闲" % (self.t, p + q, self.W)
        cost = 50 * k + 100 * p + 400 * q
        assert cost <= self.M, "t=%d: 花费 %d 但只有 %d 矿物" % (self.t, cost, self.M)
        self.M -= cost
        self.S -= k
        self.B -= k
        g = self.W - p - q
        self.ev[self.t + 2][0] += k                                   # 造工人 20 秒
        self.ev[self.t + 2][2] += k
        self.ev[self.t + 3][0] += p                                   # 补给站 30 秒
        self.ev[self.t + 3][1] += 8 * p
        self.ev[self.t + 12][0] += q                                  # 基地 120 秒
        self.ev[self.t + 12][1] += 10 * q
        self.ev[self.t + 12][2] += q
        self.ev[self.t + 1][0] += g                                   # 采矿 10 秒 +8
        self.ev[self.t + 1][3] += 8 * g
        self.t += 1
        e = self.ev.pop(self.t, [0, 0, 0, 0])
        self.W = e[0]
        self.S += e[1]
        self.B += e[2]
        self.M += e[3]


def check(path):
    T = val = None
    start = None
    acts = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('T='):
            T = int(line.split()[0].split('=')[1])
            val = int(line.split()[1].split('=')[1])
            continue
        if line.startswith('START'):
            f = [int(x) for x in line.split()[1:]]
            start = (f[0], f[1:])
            continue
        acts.append(tuple(int(x) for x in line.split()))
    s = Sim(start)
    for a in acts:
        s.act(*a)
    ok = (s.t == T and s.M == val)
    print("%-22s 起点 t=%-4d 步数=%-4d 到达 t=%-4d (=%d 秒)  2F[t]=%-16d 声明 %-16d %s"
          % (path, start[0] if start else 0, len(acts), s.t, 10 * s.t, s.M, val, "通过" if ok else "*** 不符 ***"))
    return ok


if __name__ == '__main__':
    files = sys.argv[1:] or ['data/cert/t171.txt', 'data/cert/t333.txt', 'data/cert/t500.txt']
    sys.exit(0 if all([check(f) for f in files]) else 1)
