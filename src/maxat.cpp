// 直接最大化指定时刻 T 的矿物量：全候选动作 + 大束搜索，打分用真实兑现 rollout。
// 用于在小 t 上把模型的真值逼到极限，与官方样例的上界对照。
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
// 打分：从状态 e（时刻 t）出发，用「先 a 步满额造工人(带按需补给)，再全采矿」的小网格兑现到 T
static u64 scoreTo(const St&e,int t,int T){
    u64 best=e.x; int hz=T-t; if(hz<=0) return e.x;
    for(int a=0;a<=hz;a++){
        St d=e; bool ok=true;
        for(int u=0;u<hz;u++){
            u64 k=0,p=0,q=0;
            if(u<a){ // 满额造工人；补给不足时先补一座补给站
                u64 x1=d.x;
                if(d.s<d.b && x1>=50 && d.w>=1){ p=1; x1-=50; }
                k=min(min(d.b,d.s),x1/25);
            }
            if(p+q>d.w||k>d.b||k>d.s||25*k+50*p+200*q>d.x){ok=false;break;}
            d=step(d,k,p,q);
            if(d.x>best) best=d.x;
        }
        (void)ok;
    }
    return best;
}
static void gen(const St&e,vector<array<u64,3> >&out,int FULL){
    out.clear();
    u64 kM=min(min(e.b,e.s),e.x/25);
    for(u64 k=0;k<=kM;k++){
        u64 x1=e.x-25*k, pM=min(e.w,x1/50);
        for(u64 p=0;p<=pM;p++){
            u64 x2=x1-50*p, qM=min(e.w-p,x2/200);
            for(u64 q=0;q<=qM;q++){ array<u64,3> a={k,p,q}; out.push_back(a); }
        }
    }
    (void)FULL;
}
int main(int argc,char**argv){
    int T=argc>1?atoi(argv[1]):63; int B=argc>2?atoi(argv[2]):200000;
    vector<St> cur(1); memset(&cur[0],0,sizeof(St));
    cur[0].x=25;cur[0].s=6;cur[0].w=4;cur[0].b=1;
    vector<St> nx; vector<array<u64,3> > cd;
    u64 BEST=25;
    for(int t=0;t<T;t++){
        nx.clear();
        for(size_t i=0;i<cur.size();i++){ gen(cur[i],cd,1);
            for(size_t j=0;j<cd.size();j++) nx.push_back(step(cur[i],cd[j][0],cd[j][1],cd[j][2])); }
        sort(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))<0;});
        nx.erase(unique(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))==0;}),nx.end());
        int n=nx.size();
        vector<u64> sc(n);
        #pragma omp parallel for schedule(dynamic,64)
        for(int i=0;i<n;i++) sc[i]=scoreTo(nx[i],t+1,T);
        for(int i=0;i<n;i++) if(sc[i]>BEST) BEST=sc[i];
        vector<int> v(n); iota(v.begin(),v.end(),0);
        int kk=min(B,n);
        partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return sc[a]>sc[b];});
        vector<St> ns(kk); for(int i=0;i<kk;i++) ns[i]=nx[v[i]];
        cur.swap(ns);
        if(t%10==9||t==T-1) { fprintf(stderr,"t=%d 候选=%d 束=%d 当前最好 2F[%d]>=%llu\n",t+1,n,(int)cur.size(),T,2*BEST); fflush(stderr);}
    }
    u64 fin=0; for(size_t i=0;i<cur.size();i++) if(cur[i].x>fin)fin=cur[i].x;
    if(fin>BEST)BEST=fin;
    printf("T=%d 束宽=%d  ->  2F[%d] >= %llu\n",T,B,T,2*BEST);
    return 0;
}
