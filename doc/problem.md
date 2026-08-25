# P17264 [ICPC 2017 Urumqi R] Gathering

## 题目描述

In the game of Starcraft, Terran has $1$ idle base, $4$ idle workers, $6$ supplies and $50$ minerals at the beginning. When a worker is idle, the worker can choose one of the following $3$ things to do.

*   Gathering minerals: It will cost nothing, and the worker will be busy for $10$ seconds. After that, the minerals will increase by $8$ and the worker will be idle again.
*   Building a supply depot: It will cost $100$ minerals immediately, and the worker will be busy for $30$ seconds. After that, the supplies will increase by $8$ and the worker will be idle again.
*   Building a base: It will cost $400$ minerals immediately, and the worker will be busy for $120$ seconds. After that, a new idle base will be created, and the supplies will increase by $10$, and the worker will be idle again.

When a base is idle, it can cost $50$ minerals and $1$ supply immediately, and it will be busy for $20$ seconds. After that, an idle worker will be created and the base will be idle again.

The amount of minerals and supplies cannot be negative all the time. That is, we must have enough minerals and supplies to be cost before doing anything.

Your task is to reach $m$ minerals as soon as possible.

## 输入格式

The first line is the number of test cases which is up to $52$.

For each test case, there is a line containing only 1 integer $m (m < 2^{64})$.

## 输出格式

For each test case, output the minimum time to reach $m$ minerals.

## 输入输出样例 #1

### 输入 #1

```
3
100
1000
10000
```

### 输出 #1

```
20
200
640
```

## 输入输出样例 #2

### 输入 #2

```
36
9
79
99
126
166
266
426
666
999
1899
3399
5599
9999
21999
45999
99999
316227
999999
2476413
2745943
9999999
31622776
99999999
316227766
999999999
3162277660
9641625025
10691286350
99999999999
316227766017
999999999999
9999999999999
48610229060556
59770531908338
66277611238091
99999999999999
```

### 输出 #2

```
0
10
20
30
40
70
110
160
200
300
410
530
640
800
940
1090
1320
1540
1710
1740
1990
2210
2430
2650
2880
3100
3320
3330
3770
3990
4210
4660
4970
5000
5030
5100
```