#!/usr/bin/env python3
"""从 target/improve2 的逐时刻束转储里反向匹配出一条完整动作链。

转储只存状态、不存父指针；但束在 t 时刻的每个状态一定由 t-1 时刻束里的某个状态
经候选动作 gen() 一步得到，逐层枚举匹配即可唯一还原动作序列。
"""
import struct, sys
H = 13

def loadbin(p):
    d = open(p, 'rb').read(); cnt, t0 = struct.unpack_from('<ii', d, 0)
    sz = 8 * (4 + 4 * H)
    return t0, [struct.unpack_from('<%dQ' % (4 + 4 * H), d, 8 + i * sz) for i in range(cnt)]

def unpack(v):
    return v[0], v[1], v[2], v[3], list(v[4:17]), list(v[17:30]), list(v[30:43]), list(v[43:56])

def step(st, k, p, q):
    x, s, w, b, aw, am, asu, ab = st
    aw = aw[:]; am = am[:]; asu = asu[:]; ab = ab[:]
    x -= 25 * k + 50 * p + 200 * q; s -= k; b -= k
    ab[2] += k; aw[2] += k
    aw[3] += p; asu[3] += 8 * p
    aw[12] += q; asu[12] += 10 * q; ab[12] += q
    g = w - p - q; aw[1] += g; am[1] += 4 * g
    nx = x + am[1]; ns = s + asu[1]; nw = aw[1]; nb = b + ab[1]
    for i in range(1, H - 1):
        aw[i] = aw[i + 1]; am[i] = am[i + 1]; asu[i] = asu[i + 1]; ab[i] = ab[i + 1]
    aw[H - 1] = am[H - 1] = asu[H - 1] = ab[H - 1] = 0
    aw[0] = am[0] = asu[0] = ab[0] = 0
    return (nx, ns, nw, nb, aw, am, asu, ab)

def pack(st):
    x, s, w, b, aw, am, asu, ab = st
    return tuple([x, s, w, b] + aw + am + asu + ab)

def needsup(st):
    x, s, w, b, aw, am, asu, ab = st
    v = 2 * ab[1] + ab[2] + ab[3] + 2 * b - s - asu[1] - asu[2] - asu[3]
    return (v + 7) // 8 if v > 0 else 0

def gen(st, DP, DK, DQ):
    x, s, w, b = st[0], st[1], st[2], st[3]
    out = []
    need = needsup(st); pmax = min(w, x // 50); p0 = min(need, pmax)
    for dp in range(-DP, DP + 1):
        pv = p0 + dp
        if pv < 0 or pv > pmax: continue
        x1 = x - 50 * pv
        k0 = min(b, s, x1 // 25)
        for dk in range(-DK, 1):
            kv = k0 + dk
            if kv < 0: continue
            x2 = x1 - 25 * kv
            q0 = min(w - pv, x2 // 200)
            for dq in range(-DQ, 1):
                qv = q0 + dq
                if qv < 0: continue
                if pv + qv > w: continue
                if 25 * kv + 50 * pv + 200 * qv > x: continue
                out.append((kv, pv, qv))
    return out

def main(od, t_lo, t_hi, idx_hi, DP=3, DK=2, DQ=3):
    """从 outdir 的 g<t_hi>.bin 第 idx_hi 个状态，反向还原到 g<t_lo>.bin"""
    tgt = None
    acts = []
    cur_target = None
    for t in range(t_hi, t_lo, -1):
        _, hi = loadbin('%s/g%d.bin' % (od, t))
        _, lo = loadbin('%s/g%d.bin' % (od, t - 1))
        goal = hi[idx_hi] if cur_target is None else cur_target
        found = None
        for i, v in enumerate(lo):
            st = unpack(v)
            for a in gen(st, DP, DK, DQ):
                if pack(step(st, *a)) == goal:
                    found = (i, a); break
            if found: break
        if not found:
            print("t=%d 无法匹配父状态" % t); return None
        acts.append(found[1])
        cur_target = lo[found[0]]
    acts.reverse()
    return acts

if __name__ == '__main__':
    od = sys.argv[1]; t_lo = int(sys.argv[2]); t_hi = int(sys.argv[3]); idx = int(sys.argv[4])
    a = main(od, t_lo, t_hi, idx)
    if a is None: sys.exit(1)
    out = sys.argv[5] if len(sys.argv) > 5 else '/dev/stdout'
    open(out, 'w').write('\n'.join('%d %d %d' % x for x in a) + '\n')
    print("还原 %d 步 (t=%d -> %d) -> %s" % (len(a), t_lo, t_hi, out))
