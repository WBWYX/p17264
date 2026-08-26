// 用 LP 值束搜索逐点重算 F[t]：对每个目标时刻 t，取参考轨迹在 t-W 处的状态，
// 在 [t-W, t] 上做 LP 束搜索。各个 t 相互独立，用 OpenMP 并行。
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64; typedef long double D;
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
    long long v=2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
              +2LL*(long long)e.b-(long long)e.s
              -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
    u64 need=(v>0)?(u64)((v+7)/8):0, x1=e.x;
    p=min(min(need,e.w),x1/50); x1-=50*p;
    k=min(min(e.b,e.s),x1/25); x1-=25*k;
    q=min(e.w-p,x1/200);
}
// ---- LP ----
struct LPS {
    int Hz,N,M,W; vector<D> Tb;
    D X0,S0,B0,W0,AW[H],AM[H],AS[H],AB[H];
    inline D& at(int r,int c){ return Tb[(size_t)r*W+c]; }
    int Kx(int t){return t;} int Px(int t){return Hz+t;} int Qx(int t){return 2*Hz+t;}
    D G0(int u){ D r=W0; for(int i=1;i<H&&i<=u;i++) r+=AW[i]; return r; }
    D X0c(int u){ D r=X0; for(int i=1;i<H&&i<=u;i++) r+=AM[i]; return r; }
    D S0c(int u){ D r=S0; for(int i=1;i<H&&i<=u;i++) r+=AS[i]; return r; }
    D B0c(int u){ D r=B0; for(int i=1;i<H&&i<=u;i++) r+=AB[i]; return r; }
    vector<D> coef;
    D xExpr(int t){
        coef.assign(N,0); D c0=X0c(t);
        for(int u=0;u<Hz;u++){
            D ck=4.0L*max(0,t-u-2);
            int hp=min(t-1,u+2); D cp=(hp>=u)?(D)(hp-u+1):0;
            int hq=min(t-1,u+11); D cq=(hq>=u)?(D)(hq-u+1):0;
            coef[Kx(u)]+= ck-(u<t?25.0L:0.0L);
            coef[Px(u)]+= -4*cp-(u<t?50.0L:0.0L);
            coef[Qx(u)]+= -4*cq-(u<t?200.0L:0.0L);
        }
        for(int v=0;v<t;v++) c0+=4.0L*G0(v);
        return c0;
    }
    void setState(const St&e,D lam){
        X0=(D)e.x*lam;S0=(D)e.s*lam;W0=(D)e.w*lam;B0=(D)e.b*lam;
        for(int i=0;i<H;i++){AW[i]=(D)e.aw[i]*lam;AM[i]=(D)e.am[i]*lam;AS[i]=(D)e.as_[i]*lam;AB[i]=(D)e.ab[i]*lam;}
    }
    void build(int hz){
        Hz=hz;N=3*Hz;M=4*Hz;W=N+M+1; Tb.assign((size_t)(M+1)*W,0.0L);
        for(int t=0;t<Hz;t++){
            { int r=4*t; D c0=xExpr(t);
              for(int j=0;j<N;j++) at(r,j)=-coef[j];
              at(r,Kx(t))+=25; at(r,Px(t))+=50; at(r,Qx(t))+=200; at(r,N+M)=c0; }
            { int r=4*t+1; at(r,Kx(t))+=1;
              for(int u=0;u<=t-1;u++) at(r,Kx(u))+=1;
              for(int u=0;u<=t-3;u++) at(r,Px(u))-=8;
              for(int u=0;u<=t-12;u++) at(r,Qx(u))-=10;
              at(r,N+M)=S0c(t); }
            { int r=4*t+2; at(r,Kx(t))+=1; if(t>=1) at(r,Kx(t-1))+=1;
              for(int u=0;u<=t-12;u++) at(r,Qx(u))-=1;
              at(r,N+M)=B0c(t); }
            { int r=4*t+3;
              for(int d=0;d<=2;d++) if(t-d>=0) at(r,Px(t-d))+=1;
              for(int j=0;j<=11;j++) if(t-j>=0) at(r,Qx(t-j))+=1;
              for(int u=0;u<=t-2;u++) at(r,Kx(u))-=1;
              at(r,N+M)=G0(t); }
        }
        for(int r=0;r<M;r++) at(r,N+r)=1;
        { D c0=xExpr(Hz); for(int j=0;j<N;j++) at(M,j)=-coef[j]; at(M,N+M)=c0; }
    }
    D solve(vector<D>*sol){
        vector<int> basis(M); for(int r=0;r<M;r++) basis[r]=N+r;
        for(int it=0;it<400000;it++){
            int pc=-1; D bs=-1e-9L;
            for(int j=0;j<N+M;j++) if(at(M,j)<bs){bs=at(M,j);pc=j;}
            if(pc<0) break;
            int pr=-1; D bR=0;
            for(int r=0;r<M;r++){ D a=at(r,pc); if(a>1e-11L){ D rt=at(r,N+M)/a;
                if(pr<0||rt<bR-1e-12L||(fabsl(rt-bR)<=1e-12L&&basis[r]<basis[pr])){pr=r;bR=rt;} } }
            if(pr<0) break;
            D pv=at(pr,pc); D*P=&at(pr,0);
            for(int c=0;c<W;c++) P[c]/=pv;
            for(int r=0;r<=M;r++) if(r!=pr){ D fq=at(r,pc); if(fq!=0){ D*R=&at(r,0);
                for(int c=0;c<W;c++) R[c]-=fq*P[c]; } }
            basis[pr]=pc;
        }
        if(sol){ sol->assign(N,0); for(int r=0;r<M;r++) if(basis[r]<N) (*sol)[basis[r]]=at(r,N+M); }
        return at(M,N+M);
    }
};
static u64 lpBeamOne(const St&s0,int Hz,int B,int DW){
    vector<St> cur(1,s0); u64 BEST=s0.x;
    for(int t=0;t<Hz;t++){
        int hz=Hz-t, n=cur.size();
        vector<St> nx;
        for(int i=0;i<n;i++){
            const St&c=cur[i];
            vector<D> sol; D lam=1.0L; if(c.x>2000000ULL) lam=1000000.0L/(D)c.x;
            { LPS L; L.setState(c,lam); L.build(hz); L.solve(&sol); }
            long long bk=(long long)floorl(sol[0]/lam),bp=(long long)floorl(sol[hz]/lam),bq=(long long)floorl(sol[2*hz]/lam);
            if(bk<0)bk=0; if(bp<0)bp=0; if(bq<0)bq=0;
            for(int dk=-DW;dk<=DW+1;dk++) for(int dp=-DW;dp<=DW+1;dp++) for(int dq=-DW;dq<=DW+1;dq++){
                long long k=bk+dk,p=bp+dp,q=bq+dq;
                if(k<0||p<0||q<0) continue;
                if((u64)k>c.b||(u64)k>c.s||(u64)(p+q)>c.w) continue;
                if((unsigned long long)(25*k+50*p+200*q)>c.x) continue;
                nx.push_back(step(c,(u64)k,(u64)p,(u64)q));
            }
        }
        sort(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))<0;});
        nx.erase(unique(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))==0;}),nx.end());
        int m=nx.size(); if(!m) break;
        vector<D> sc(m);
        if(hz-1<=0){ for(int i=0;i<m;i++) sc[i]=(D)nx[i].x; }
        else for(int i=0;i<m;i++){ LPS L; D l=1.0L; if(nx[i].x>2000000ULL) l=1000000.0L/(D)nx[i].x;
            L.setState(nx[i],l); L.build(hz-1); sc[i]=L.solve(0)/l; }
        for(int i=0;i<m;i++) if(nx[i].x>BEST) BEST=nx[i].x;
        vector<int> v(m); iota(v.begin(),v.end(),0);
        int kk=min(B,m);
        partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return sc[a]>sc[b];});
        vector<St> ns(kk); for(int i=0;i<kk;i++) ns[i]=nx[v[i]];
        cur.swap(ns);
    }
    return BEST;
}
int main(int argc,char**argv){
    const char*bf=argv[1]; int idx=atoi(argv[2]);
    int T1=atoi(argv[3]), T2=atoi(argv[4]), W=atoi(argv[5]);
    int B=argc>6?atoi(argv[6]):16, DW=argc>7?atoi(argv[7]):1;
    const char*actf=(argc>8)?argv[8]:0; int asi=(argc>9)?atoi(argv[9]):0;
    FILE*f=fopen(bf,"rb"); int cnt,tb; if(fread(&cnt,4,1,f)!=1||fread(&tb,4,1,f)!=1)return 1;
    vector<St> all(cnt);
    for(int i=0;i<cnt;i++) if(fread(&all[i],sizeof(St),1,f)!=1){cnt=i;break;}
    fclose(f);
    // 参考轨迹：束状态 idx -> Pi3 推进到 tg -> （可选）改进段动作 -> 其后继续 Pi3
    int tg=-1; vector<array<u64,3> > ACT;
    if(actf){ FILE*fa=fopen(actf,"r"); if(!fa){fprintf(stderr,"打不开 %s\n",actf);return 1;}
        int ns=0,nst=0; if(fscanf(fa,"%d %d %d",&tg,&ns,&nst)!=3){fprintf(stderr,"动作文件头有误\n");return 1;}
        vector<array<u64,3> > a2((size_t)ns*nst);
        for(size_t i=0;i<a2.size();i++){ unsigned long long u1,u2,u3;
            if(fscanf(fa,"%llu %llu %llu",&u1,&u2,&u3)!=3){fprintf(stderr,"动作文件长度不足\n");return 1;}
            a2[i]={u1,u2,u3}; }
        fclose(fa);
        for(int i=0;i<nst;i++) ACT.push_back(a2[(size_t)asi*nst+i]); }
    int TMAX=T2+2; vector<St> traj(TMAX+1);
    { St c=all[idx];
      for(int t=tb;t<=TMAX;t++){ traj[t]=c; if(t>=TMAX) break;
        u64 k,p,q;
        if(tg>=0&&t>=tg&&t<tg+(int)ACT.size()){ k=ACT[t-tg][0];p=ACT[t-tg][1];q=ACT[t-tg][2];
            if(p+q>c.w||k>c.b||k>c.s||25*k+50*p+200*q>c.x){fprintf(stderr,"改进段动作 t=%d 不可行\n",t);return 1;} }
        else growAct(c,k,p,q);
        c=step(c,k,p,q); } }
    printf("# 参考轨迹: %s#%d, t0=%d%s, 窗口 W=%d, 束宽=%d, DELW=%d\n",bf,idx,tb,
           actf?"(+改进段)":"",W,B,DW);
    fflush(stdout);
    int n=T2-T1+1; vector<u64> res(n,0);
    double w0=omp_get_wtime();
    #pragma omp parallel for schedule(dynamic,1)
    for(int i=0;i<n;i++){
        int T=T1+i, t0=max(tb,T-W);
        res[i]=lpBeamOne(traj[t0],T-t0,B,DW);
    }
    for(int i=0;i<n;i++) printf("%d %llu\n",T1+i,res[i]);
    fprintf(stderr,"完成 %d 个时刻, 用时 %.0fs\n",n,omp_get_wtime()-w0);
    return 0;
}
