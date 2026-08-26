// 阶段 1：从束转储重建交接轨迹。
// 对每个束状态按自校正成长规则(pol=3)推进，按到达目标时刻 T 的成长矿物量排序，
// 导出前 N 个状态（已推进到 t0）供 LP 使用；同时报告纯成长值与 longVal 排序的重合度。
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
static inline bool okA(const St&e,u64 k,u64 p,u64 q){
    return p+q<=e.w && k<=e.b && k<=e.s && 25*k+50*p+200*q<=e.x;
}
// 自校正成长（= gen_table act(pol=3), SUPMODE=7）
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
// gen_table 的 longVal（double 松弛，SUPMODE>=4 分支）
int NLONG=300,NCASH=22;
static double longVal(const St&e){
    double x=e.x,ss=e.s,w=e.w,b=e.b;
    double AW[H],AM[H],AS[H],AB[H];
    for(int i=0;i<H;i++){AW[i]=e.aw[i];AM[i]=e.am[i];AS[i]=e.as_[i];AB[i]=e.ab[i];}
    int h=0; double sc=0;
    for(int j=0;j<NLONG+NCASH;j++){
        bool grow=j<NLONG;
        double wf=w,p=0,k=0,q=0;
        if(grow||j<NLONG+8){
            double a1=AB[(h+1)%H],a3=AB[(h+3)%H]; double need=(a1-9*a3)/8;
            if(need<0)need=0; if(ss<b&&(b-ss)/8>need) need=(b-ss)/8;
            p=min(min(need,wf),x/50); x-=50*p; wf-=p;
            int i3=(h+3)%H; AW[i3]+=p; AS[i3]+=8*p;
        }
        k=min(min(b,ss),x/25); x-=25*k; ss-=k; b-=k;
        { int i2=(h+2)%H; AB[i2]+=k; AW[i2]+=k; }
        if(grow){ q=min(wf,x/200); x-=200*q; wf-=q;
            int i12=(h+12)%H; AW[i12]+=q; AS[i12]+=10*q; AB[i12]+=q; }
        int i1=(h+1)%H; AW[i1]+=wf; AM[i1]+=4*wf;
        x+=AM[i1]; ss+=AS[i1]; w=AW[i1]; b+=AB[i1];
        AW[i1]=AM[i1]=AS[i1]=AB[i1]=0; h=i1;
        if(x>1e250||w>1e250){ double f=1e-250; x*=f;ss*=f;w*=f;b*=f;
            for(int i=0;i<H;i++){AW[i]*=f;AM[i]*=f;AS[i]*=f;AB[i]*=f;}
            sc+=log(1e250); }
    }
    return sc+log(x+1.0);
}
int main(int argc,char**argv){
    if(argc<5){ printf("用法: rank <beam.bin> <T> <t0> <N> [out.bin]\n"); return 1; }
    const char*fn=argv[1]; int T=atoi(argv[2]); int t0=atoi(argv[3]); int N=atoi(argv[4]);
    const char*of=argc>5?argv[5]:0;
    FILE*f=fopen(fn,"rb"); if(!f){printf("no %s\n",fn);return 1;}
    int cnt,tb; if(fread(&cnt,4,1,f)!=1||fread(&tb,4,1,f)!=1)return 1;
    vector<St> all(cnt);
    for(int i=0;i<cnt;i++) if(fread(&all[i],sizeof(St),1,f)!=1){cnt=i;all.resize(cnt);break;}
    fclose(f);
    if(t0<tb){ printf("t0=%d < 束时刻 %d\n",t0,tb); return 1; }
    vector<St> at0(cnt); vector<u64> gv(cnt,0); vector<double> lv(cnt);
    int alive=0;
    for(int i=0;i<cnt;i++){
        St c=all[i]; bool good=true;
        for(int t=tb;t<T;t++){ u64 k,p,q; growAct(c,k,p,q);
            if(!okA(c,k,p,q)){good=false;break;}
            if(t==t0) at0[i]=c;
            c=step(c,k,p,q); }
        if(t0==T) at0[i]=c;
        if(good){ gv[i]=c.x; alive++; }
        lv[i]=longVal(all[i]);
    }
    vector<int> idx(cnt); iota(idx.begin(),idx.end(),0);
    sort(idx.begin(),idx.end(),[&](int a,int b){return gv[a]>gv[b];});
    vector<int> idl(cnt); iota(idl.begin(),idl.end(),0);
    sort(idl.begin(),idl.end(),[&](int a,int b){return lv[a]>lv[b];});
    printf("%s: %d 状态 t0dump=%d -> 成长到 T=%d (存 t0=%d), 可行 %d\n",fn,cnt,tb,T,t0,alive);
    printf("  纯成长 x(%d)*2 最大 = %llu  中位 = %llu\n",T,2*gv[idx[0]],2*gv[idx[cnt/2]]);
    printf("  longVal 最大 = %.6f ; 其纯成长值*2 = %llu\n",lv[idl[0]],2*gv[idl[0]]);
    // 两种排序的前 N 交集大小
    {   set<int> s1(idx.begin(),idx.begin()+min(N,cnt));
        int inter=0; for(int i=0;i<min(N,cnt);i++) if(s1.count(idl[i]))inter++;
        printf("  前 %d 名交集 = %d\n",min(N,cnt),inter); }
    if(of){
        // 并集：成长值前 N + longVal 前 N
        vector<int> u;
        for(int i=0;i<min(N,cnt);i++) u.push_back(idx[i]);
        for(int i=0;i<min(N,cnt);i++) u.push_back(idl[i]);
        sort(u.begin(),u.end()); u.erase(unique(u.begin(),u.end()),u.end());
        // 只保留可行的
        vector<int> v; for(size_t i=0;i<u.size();i++) if(gv[u[i]]>0) v.push_back(u[i]);
        sort(v.begin(),v.end(),[&](int a,int b){return gv[a]>gv[b];});
        FILE*fo=fopen(of,"wb"); int c2=v.size(); fwrite(&c2,4,1,fo); fwrite(&t0,4,1,fo);
        for(int i=0;i<c2;i++) fwrite(&at0[v[i]],sizeof(St),1,fo); fclose(fo);
        printf("  写出 %d 个 t0=%d 的状态 -> %s\n",c2,t0,of);
    }
    return 0;
}
