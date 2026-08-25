// 用自校正成长规则跑到 T0，导出状态；同时给出 (a,b) 网格兑现到 T0+W 的最好值
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64;
static const int H=13;
struct St{u64 x,s,w,b;u64 aw[H],am[H],as_[H],ab[H];};
static St step(const St&t,u64 k,u64 p,u64 q){
    St c=t; c.x-=25*k+50*p+200*q; c.s-=k; c.b-=k;
    c.ab[2]+=k;c.aw[2]+=k; c.aw[3]+=p;c.as_[3]+=8*p;
    c.aw[12]+=q;c.as_[12]+=10*q;c.ab[12]+=q;
    u64 g=t.w-p-q; c.aw[1]+=g; c.am[1]+=4*g;
    u64 nx=c.x+c.am[1],ns=c.s+c.as_[1],nw=c.aw[1],nb=c.b+c.ab[1];
    for(int i=1;i+1<H;i++){c.aw[i]=c.aw[i+1];c.am[i]=c.am[i+1];c.as_[i]=c.as_[i+1];c.ab[i]=c.ab[i+1];}
    c.aw[H-1]=c.am[H-1]=c.as_[H-1]=c.ab[H-1]=0;
    c.aw[0]=c.am[0]=c.as_[0]=c.ab[0]=0;
    c.x=nx;c.s=ns;c.w=nw;c.b=nb; return c;
}
static void growAct(const St&e,u64&k,u64&p,u64&q){
    long long v = 2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                + 2LL*(long long)e.b-(long long)e.s
                -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
    u64 need=(v>0)?(u64)((v+7)/8):0;
    u64 x1=e.x;
    p=min(min(need,e.w),x1/50); x1-=50*p;
    k=min(min(e.b,e.s),x1/25); x1-=25*k;
    q=min(e.w-p,x1/200);
}
int GNC=98,KAC=64,CASHFF=0;
// pol 0=全采矿 1=只造工人 2=猛建补给+造工人 3=停基地但补给维持最小(自校正)
static void cashAct(const St&e,int pol,u64&k,u64&p,u64&q){
    k=p=q=0; if(pol==0) return; u64 x1=e.x;
    if(pol==3){ long long v=2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                          +2LL*(long long)e.b-(long long)e.s
                          -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
                u64 need=(v>0)?(u64)((v+7)/8):0;
                p=min(min(need,e.w),x1/50); x1-=50*p; }
    else if(pol==2){ u64 need;
        if(CASHFF){ // 推广自校正: 8p = beta*(ab1+ab3) + b + ab1 + (b+ab2) - s - (as1+as2+as3)
            long long B=GNC; // beta = GNC/20
            long long v = B*((long long)e.ab[1]+(long long)e.ab[3])/20
                        + (long long)e.b + (long long)e.ab[1] + (long long)e.b + (long long)e.ab[2]
                        - (long long)e.s
                        - (long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
            need=(v>0)?(u64)((v+7)/8):0;
        } else { u64 Dm=((u64)e.b*GNC+19)/20; need=(Dm>e.s)?(Dm-e.s+7)/8:0; }
        p=min(min(need,e.w),x1/50); x1-=50*p; }
    k=min(min(e.b,e.s),x1/25); if(pol!=3) k=k*KAC/64;
}
int main(int argc,char**argv){
    int T0=argc>1?atoi(argv[1]):200, WIN=argc>2?atoi(argv[2]):48; int CMAX=argc>5?atoi(argv[5]):0; int ASTEP=argc>7?atoi(argv[7]):1; int AMAX=argc>8?atoi(argv[8]):30; int BMAX=argc>9?atoi(argv[9]):24;
    if(argc>3)GNC=atoi(argv[3]); if(argc>4)KAC=atoi(argv[4]); if(argc>6)CASHFF=atoi(argv[6]);
    St c; memset(&c,0,sizeof(c)); c.x=25;c.s=6;c.w=4;c.b=1;
    for(int t=0;t<T0;t++){ u64 k,p,q; growAct(c,k,p,q); c=step(c,k,p,q); }
    FILE*f=fopen("st.bin","wb"); int one=1; fwrite(&one,4,1,f); fwrite(&T0,4,1,f); fwrite(&c,sizeof(St),1,f); fclose(f);
    // 与生成器一致：沿成长轨迹每一步都做 (a,b) 网格兑现，取到达 T0+WIN 的最大值
    u64 best=0; St g=c;
    for(int t=T0;t<=T0+WIN;t++){
        if(t==T0+WIN){ if(g.x>best)best=g.x; break; }
        for(int cc=0;cc<=CMAX;cc+=ASTEP) for(int a=0;a<=AMAX;a+=ASTEP) for(int bb=0;bb<=BMAX;bb+=ASTEP){
            St d=g;
            for(int u=1;u<=T0+WIN-t;u++){ int pol=(u<=cc)?3:((u<=cc+a)?2:((u<=cc+a+bb)?1:0));
                u64 k,p,q; cashAct(d,pol,k,p,q);
                if(25*k+50*p>d.x||p>d.w||k>d.b||k>d.s) break;
                d=step(d,k,p,q); if(t+u==T0+WIN&&d.x>best)best=d.x; }
        }
        u64 k,p,q; growAct(g,k,p,q); g=step(g,k,p,q);
    }
    printf("t0=%d 状态 x=%llu s=%llu b=%llu w=%llu\n",T0,c.x,c.s,c.b,c.w);
    printf("(a,b) 网格(步长1, GNC=%d KAC=%d) 兑现到 t=%d 的最好值 x*2 = %llu\n",GNC,KAC,T0+WIN,2*best);
    return 0;
}
