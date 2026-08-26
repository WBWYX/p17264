// 成长期 rollout 策略改进：以「按自校正成长策略精确推进到 T 的 u64 矿物量」为打分，
// 在 [tb,TG] 做束搜索。打分函数就是纯成长策略的真实值，故结果恒 >= 纯成长（rollout 改进）。
#include <bits/stdc++.h>
#include <omp.h>
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
static inline u64 needSup(const St&e){
    long long v = 2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                + 2LL*(long long)e.b-(long long)e.s
                -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
    return (v>0)?(u64)((v+7)/8):0;
}
static void growAct(const St&e,u64&k,u64&p,u64&q){
    u64 need=needSup(e), x1=e.x;
    p=min(min(need,e.w),x1/50); x1-=50*p;
    k=min(min(e.b,e.s),x1/25); x1-=25*k;
    q=min(e.w-p,x1/200);
}
// 纯成长打分：从 e（时刻 t）推进到 T
static u64 growTo(St c,int t,int T){
    for(;t<T;t++){ u64 k,p,q; growAct(c,k,p,q); c=step(c,k,p,q); }
    return c.x;
}
int DP=3,DK=2,DQ=3;
static void gen(const St&e,vector<array<u64,3> >&out){
    out.clear();
    u64 need=needSup(e);
    u64 pmax=min(e.w,e.x/50);
    u64 p0=min(need,pmax);
    for(int dp=-DP;dp<=DP;dp++){
        long long pv=(long long)p0+dp; if(pv<0||(u64)pv>pmax) continue;
        u64 p=(u64)pv, x1=e.x-50*p;
        u64 k0=min(min(e.b,e.s),x1/25);
        for(int dk=-DK;dk<=0;dk++){
            long long kv=(long long)k0+dk; if(kv<0) continue;
            u64 k=(u64)kv, x2=x1-25*k;
            u64 q0=min(e.w-p,x2/200);
            for(int dq=-DQ;dq<=0;dq++){
                long long qv=(long long)q0+dq; if(qv<0) continue;
                u64 q=(u64)qv;
                if(p+q>e.w) continue;
                if(25*k+50*p+200*q>e.x) continue;
                array<u64,3> a={k,p,q}; out.push_back(a);
            }
        }
    }
}

// 逐时刻转储版：在 [tb,TG] 做 rollout 改进束搜索，每步把整个束写到 outdir/g<t>.bin
int main(int argc,char**argv){
    if(argc<6){printf("用法: improve2 <in.bin> <TG> <T打分> <B> <outdir> [DP DK DQ] [K初始]\n");return 1;}
    const char*fn=argv[1]; int TG=atoi(argv[2]); int T=atoi(argv[3]);
    int B=atoi(argv[4]); const char*od=argv[5];
    if(argc>6)DP=atoi(argv[6]); if(argc>7)DK=atoi(argv[7]); if(argc>8)DQ=atoi(argv[8]);
    int K0=argc>9?atoi(argv[9]):B;
    FILE*f=fopen(fn,"rb"); if(!f){printf("no %s\n",fn);return 1;}
    int cnt,tb; if(fread(&cnt,4,1,f)!=1||fread(&tb,4,1,f)!=1)return 1;
    vector<St> cur(cnt);
    for(int i=0;i<cnt;i++) if(fread(&cur[i],sizeof(St),1,f)!=1){cnt=i;cur.resize(cnt);break;}
    fclose(f);
    u64 base=0; for(int i=0;i<cnt;i++){ u64 v=growTo(cur[i],tb,T); if(v>base)base=v; }
    fprintf(stderr,"%s: %d 状态 tb=%d -> TG=%d, 打分时域 T=%d, 束宽 %d, 候选(%d,%d,%d)\n",
            fn,cnt,tb,TG,T,B,DP,DK,DQ);
    fprintf(stderr,"基线(纯成长) x(%d)*2 = %llu\n",T,2*base);
    {   vector<pair<u64,int> > v(cnt);
        #pragma omp parallel for schedule(dynamic,16)
        for(int i=0;i<cnt;i++) v[i]=make_pair(growTo(cur[i],tb,T),i);
        sort(v.begin(),v.end(),greater<pair<u64,int> >());
        int kk=min(K0,cnt); vector<St> ns(kk);
        for(int i=0;i<kk;i++) ns[i]=cur[v[i].second];
        cur.swap(ns); }
    auto dump=[&](int t){ char p[512]; snprintf(p,sizeof(p),"%s/g%d.bin",od,t);
        FILE*fo=fopen(p,"wb"); if(!fo)return; int c2=cur.size();
        fwrite(&c2,4,1,fo); fwrite(&t,4,1,fo);
        for(int i=0;i<c2;i++) fwrite(&cur[i],sizeof(St),1,fo); fclose(fo); };
    dump(tb);
    double w0=omp_get_wtime();
    vector<St> nx; vector<array<u64,3> > cd;
    for(int t=tb;t<TG;t++){
        nx.clear();
        for(size_t i=0;i<cur.size();i++){ gen(cur[i],cd);
            for(size_t j=0;j<cd.size();j++) nx.push_back(step(cur[i],cd[j][0],cd[j][1],cd[j][2])); }
        sort(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))<0;});
        nx.erase(unique(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))==0;}),nx.end());
        int n=nx.size(); if(!n){fprintf(stderr,"t=%d 无后继\n",t);break;}
        vector<pair<u64,int> > sc(n);
        #pragma omp parallel for schedule(dynamic,16)
        for(int i=0;i<n;i++) sc[i]=make_pair(growTo(nx[i],t+1,T),i);
        int kk=min(B,n);
        partial_sort(sc.begin(),sc.begin()+kk,sc.end(),greater<pair<u64,int> >());
        vector<St> ns(kk); for(int i=0;i<kk;i++) ns[i]=nx[sc[i].second];
        cur.swap(ns);
        dump(t+1);
        if((t+1)%20==0) fprintf(stderr,"  t=%d 候选=%d 束=%d 成长值*2=%llu (%.0fs)\n",
            t+1,n,(int)cur.size(),2*sc[0].first,omp_get_wtime()-w0);
    }
    fprintf(stderr,"完成, 用时 %.0fs\n",omp_get_wtime()-w0);
    return 0;
}
