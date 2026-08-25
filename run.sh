#!/usr/bin/env bash
# P17264 打表收尾流水线。用法: bash run.sh {build|selftest|check}
set -u
cd "$(dirname "$0")"
mkdir -p build out

CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:--O2 -std=c++17}

build() {
    echo "=== 编译 ==="
    for f in exact lp dumpst allcash finalcash gen_table; do
        printf '  %-12s' "$f"
        if $CXX $CXXFLAGS "src/$f.cpp" -o "build/$f" 2> "build/$f.log"; then
            echo "ok"
        else
            echo "失败，见 build/$f.log"; head -5 "build/$f.log"; return 1
        fi
    done
}

# 生成一个状态文件：t0 与四元组 (x,s,w,b)，流水线全零
mkstate() {  # mkstate <out.bin> <t0> <x> <s> <w> <b>
    python3 - "$@" <<'PY'
import struct, sys
out, t0, x, s, w, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
d = struct.pack('<ii', 1, t0) + b''.join(struct.pack('<Q', v) for v in [x, s, w, b] + [0]*52)
open(out, 'wb').write(d)
PY
}

selftest() {
    build || return 1
    local fail=0

    echo
    echo "=== 自校验 1/3：精确全枚举 vs 表（表在 t<=25 应逐项等于真实最优）==="
    ./build/exact 25 2000000 1 0 > out/exact25.txt 2>&1
    python3 - <<'PY' || fail=1
import re
tab = [int(v) for v in open('data/table_best.txt').read().replace('ULL','').split(',') if v.strip()]
got = {}
for line in open('out/exact25.txt', encoding='utf-8', errors='replace'):
    m = re.search(r't=\s*(\d+).*?maxX\*2=(\d+)', line)
    if m: got[int(m.group(1))] = int(m.group(2))
ok = True
# 全枚举 + 支配剪枝得到的是真值；表必须逐项相符
for t, want in [(10, 384), (20, 1014), (25, 1414)]:
    g = got.get(t)
    s = 'ok' if (g == want and tab[t]*2 == want) else '*** 不符 ***'
    print("  t=%-3d 全枚举 2F=%-8s 期望 %-8d 表中 %-8d %s" % (t, g, want, tab[t]*2, s))
    if g != want or tab[t]*2 != want: ok = False
# 逐项比对 t=1..25
bad = [t for t in got if got[t] != tab[t]*2]
print("  t=1..25 逐项比对：%s" % ("全部一致" if not bad else "*** 不一致于 %s ***" % bad))
if bad: ok = False
raise SystemExit(0 if ok else 1)
PY

    echo
    echo "=== 自校验 2/3：LP 必须是有效上界，且小窗口下取到真值 ==="
    mkstate out/init.bin 0 25 6 4 1
    python3 - <<'PY' || fail=1
import subprocess, re
def lpval(f, T):
    r = subprocess.run(['./build/lp', f, str(T), '1'], capture_output=True, text=True, errors='replace')
    m = re.search(r'\*2 = ([0-9.]+)', r.stdout)
    return float(m.group(1)) if m else -1
ok = True
# T=10/20：LP 恰好取到整数最优。T=30：LP=1968 而真值 1920，2.5% 的整数性缺口——
# 这不是错误，正是「小尺度下 LP 不可信、大尺度下 LP 才精确」这一核心判断的实证。
for T, exact_val, expect_eq in [(10, 384.0, True), (20, 1014.0, True), (30, 1920.0, False)]:
    v = lpval('out/init.bin', T)
    if v < exact_val - 1e-6:
        print("  *** T=%d LP=%.3f 低于真值 %.0f，LP 不再是上界，工具有误 ***" % (T, v, exact_val)); ok = False
    elif expect_eq and abs(v - exact_val) > 1e-6:
        print("  *** T=%d LP=%.3f 期望恰为 %.0f ***" % (T, v, exact_val)); ok = False
    else:
        gap = (v - exact_val) / exact_val
        print("  T=%-3d LP 2F=%-10.3f 真值 %-8.0f 整数性缺口 %.2e" % (T, v, exact_val, gap))
raise SystemExit(0 if ok else 1)
PY

    echo
    echo "=== 自校验 3/3：LP 对初始状态的一次齐次性 ==="
    mkstate out/init_s.bin 0 25000000 6000000 4000000 1000000
    python3 - <<'PY' || fail=1
import subprocess, re
def lpval(f, T):
    r = subprocess.run(['./build/lp', f, str(T), '1'], capture_output=True, text=True, errors='replace')
    m = re.search(r'x\((\d+)\) = ([0-9.]+)', r.stdout)
    return float(m.group(2)) if m else -1
ok = True
for T in (20, 30):
    a, b = lpval('out/init.bin', T), lpval('out/init_s.bin', T)
    ratio = b / a if a else 0
    s = 'ok' if abs(ratio - 1e6) / 1e6 < 1e-9 else '*** 不符 ***'
    print("  T=%-3d 缩放后/原 = %.6f  期望 1000000  %s" % (T, ratio, s))
    if abs(ratio - 1e6) / 1e6 >= 1e-9: ok = False
raise SystemExit(0 if ok else 1)
PY

    echo
    if [ "$fail" -eq 0 ]; then echo "=== 自校验全部通过 ==="; else echo "=== 自校验存在失败项，先修好再继续 ==="; fi
    return $fail
}

check() { python3 src/check.py "${1:-data/table_best.txt}"; }

case "${1:-check}" in
    build)    build ;;
    selftest) selftest ;;
    check)    check "${2:-data/table_best.txt}" ;;
    *)        echo "用法: bash run.sh {build|selftest|check [表文件]}"; exit 1 ;;
esac
