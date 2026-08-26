// P17264 [ICPC 2017 Urumqi R] Gathering —— 提交程序
//
// 离线打表 + 在线二分。表中 F[t] 是「时刻 t（=10t 秒）可达的最大矿物量的一半」。
// 矿物量恒为偶数（初始 50，采矿 +8，成本 100/400/50），折半后 m<2^64 也不会溢出。
// 答案 = 10 * min{ t : F[t] >= ceil(m/2) }。
//
// 表由 data/table_best.txt 内嵌（见 tools/embed_table.py 生成）。无任何特判。
#include <cstdio>
#include <cstdint>

typedef unsigned long long u64;

static const u64 F[] = {
#include "table_inc.h"
};
static const int NT = sizeof(F) / sizeof(F[0]);

int main() {
    int T;
    if (scanf("%d", &T) != 1) return 0;
    while (T--) {
        u64 m;
        if (scanf("%llu", &m) != 1) break;
        u64 need = m / 2 + (m & 1);          // ceil(m/2)，不要写成 (m+1)/2（m=2^64-1 会溢出）
        int lo = 0, hi = NT - 1;
        while (lo < hi) {                     // F 非降，二分首个 F[t] >= need
            int mid = lo + (hi - lo) / 2;
            if (F[mid] >= need) hi = mid; else lo = mid + 1;
        }
        printf("%d\n", 10 * lo);
    }
    return 0;
}
