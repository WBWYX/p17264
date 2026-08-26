#!/usr/bin/env python3
"""把「束状态 -> Pi3 成长 -> rollout 改进 -> LP 束兑现」拼成一条完整动作序列证书。
证书以 START 行给出起点束状态（折半口径），后随每步 k p q。"""
import struct, sys
H = 13


def load(p):
    d = open(p, 'rb').read(); cnt, t0 = struct.unpack('<ii', d[:8])
    return t0, [list(struct.unpack('<56Q', d[8 + 448 * i:8 + 448 * (i + 1)])) for i in range(cnt)]


def unpack(v):
    return dict(x=v[0], s=v[1], w=v[2], b=v[3],
                aw=list(v[4:17]), am=list(v[17:30]), as_=list(v[30:43]), ab=list(v[43:56]))


def pack(e):
    return [e['x'], e['s'], e['w'], e['b']] + e['aw'] + e['am'] + e['as_'] + e['ab']


def step(st, k, p, q):
    c = {a: (list(st[a]) if isinstance(st[a], list) else st[a]) for a in st}
    c['x'] -= 25 * k + 50 * p + 200 * q; c['s'] -= k; c['b'] -= k
    c['ab'][2] += k; c['aw'][2] += k; c['aw'][3] += p; c['as_'][3] += 8 * p
    c['aw'][12] += q; c['as_'][12] += 10 * q; c['ab'][12] += q
    g = st['w'] - p - q; c['aw'][1] += g; c['am'][1] += 4 * g
    nx = c['x'] + c['am'][1]; ns = c['s'] + c['as_'][1]; nw = c['aw'][1]; nb = c['b'] + c['ab'][1]
    for a in ('aw', 'am', 'as_', 'ab'):
        for i in range(1, H - 1):
            c[a][i] = c[a][i + 1]
        c[a][H - 1] = 0; c[a][0] = 0
    c['x'] = nx; c['s'] = ns; c['w'] = nw; c['b'] = nb
    return c


def growAct(e):
    v = 2 * e['ab'][1] + e['ab'][2] + e['ab'][3] + 2 * e['b'] - e['s'] - e['as_'][1] - e['as_'][2] - e['as_'][3]
    need = (v + 7) // 8 if v > 0 else 0
    x1 = e['x']
    p = min(need, e['w'], x1 // 50); x1 -= 50 * p
    k = min(e['b'], e['s'], x1 // 25); x1 -= 25 * k
    q = min(e['w'] - p, x1 // 200)
    return k, p, q


def readacts(fn):
    tok = open(fn).read().split()
    return tok


def growacts(fn, si):
    tok = readacts(fn); tb = int(tok[0]); ns = int(tok[1]); nst = int(tok[2])
    b = 3 + 3 * nst * si
    return tb, [tuple(int(x) for x in tok[b + 3 * i:b + 3 * i + 3]) for i in range(nst)]


def cashacts(fn):
    tok = readacts(fn); t0 = int(tok[0]); si = int(tok[1]); n = int(tok[2]); val = int(tok[3])
    return t0, si, val, [tuple(int(x) for x in tok[4 + 3 * i:4 + 3 * i + 3]) for i in range(n)]


def build(out, beamfile, tgrow, growfile, cashfile, T):
    tb, beam = load(beamfile)
    t0c, si, val, CA = cashacts(cashfile)
    acts = []
    if growfile:
        tg, GA = growacts(growfile, si)
        # 找出哪个束状态按 Pi3 推进 (tg-tb) 步后可作为改进段起点
        src = None
        for idx, v in enumerate(beam):
            e = unpack(v); seq = []
            for t in range(tb, tg):
                a = growAct(e); seq.append(a); e = step(e, *a)
            e2 = dict(e); ok = True
            for a in GA:
                if a[0] > e2['b'] or a[0] > e2['s'] or a[1] + a[2] > e2['w'] or 25 * a[0] + 50 * a[1] + 200 * a[2] > e2['x']:
                    ok = False; break
                e2 = step(e2, *a)
            if ok and pack(e2) == pack(unpack(load('out/imp_probe.bin')[1][si])) if False else False:
                pass
            if ok:
                # 用兑现段验证：接着走 CA 应得到 val
                e3 = dict(e2); good = True
                for a in CA:
                    if a[0] > e3['b'] or a[0] > e3['s'] or a[1] + a[2] > e3['w'] or 25 * a[0] + 50 * a[1] + 200 * a[2] > e3['x']:
                        good = False; break
                    e3 = step(e3, *a)
                if good and e3['x'] == val:
                    src = idx; acts = seq + GA + CA; break
        if src is None:
            print("找不到匹配的束状态"); return None
    else:
        src = 0
        e = unpack(beam[0])
        for a in CA:
            e = step(e, *a)
        assert e['x'] == val, (e['x'], val)
        acts = CA
    st = beam[src]
    with open(out, 'w') as f:
        f.write("# P17264 动作序列证书。起点为 %s 中第 %d 个束状态（gen_table 从真实初始状态推进得到，折半口径）\n" % (beamfile, src))
        f.write("# START 行: t0 x s w b aw[0..12] am[0..12] as[0..12] ab[0..12]\n")
        f.write("# 其后每行 k p q = 本步 (基地造工人数, 建补给站数, 建基地数); 时间单位 10 秒\n")
        f.write("T=%d minerals=%d\n" % (T, 2 * val))
        f.write("START %d %s\n" % (tb, ' '.join(str(v) for v in st)))
        for k, p, q in acts:
            f.write("%d %d %d\n" % (k, p, q))
    print("%s: 起点 %s#%d t0=%d, %d 步 -> t=%d, 2F=%d" % (out, beamfile, src, tb, len(acts), tb + len(acts), 2 * val))
    return src


if __name__ == '__main__':
    build('data/cert/t171.txt', 'out/top171_beam105.bin', None, None, 'out/c171_beam105.txt', 171)
    build('data/cert/t333.txt', 'data/beam105.bin', 180, 'out/a333_grow.txt', 'out/a333_cash.txt', 333)
