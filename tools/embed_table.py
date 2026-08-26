#!/usr/bin/env python3
"""把 data/table_best.txt 转成 src/table_inc.h（提交程序内嵌用）。"""
import sys
src = sys.argv[1] if len(sys.argv) > 1 else 'data/table_best.txt'
dst = sys.argv[2] if len(sys.argv) > 2 else 'src/table_inc.h'
F = [int(x) for x in open(src).read().replace('ULL', '').split(',') if x.strip()]
with open(dst, 'w') as f:
    for i in range(0, len(F), 8):
        f.write(','.join(str(v) + 'ULL' for v in F[i:i + 8]) + ',\n')
print("内嵌 %d 项 -> %s (%d 字节)" % (len(F), dst, len(open(dst).read())))
