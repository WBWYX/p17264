// 阶段 2-3：逐状态解 LP（缩放）+ 带进位取整 + 坐标下降修复 + u64 精确重放
// 用法: drive <statefile.bin> <T> <K> <MODE> [CW] [GNC]
//   MODE 位掩码: 1=LP 上界排序  2=固定策略族(cc,a,b)网格基线  4=LP 驱动取整
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64;
typedef long double D;
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
int GNC=98;
// pol: 0=全采矿 1=只造工人 2=猛补给(GNC)+造工人 3=完全成长 4=自校正补给+造工人(停基地)
static void act(const St&e,int pol,u64&k,u64&p,u64&q){
    k=p=q=0; if(pol==0) return;
    u64 x1=e.x, need=0;
    if(pol==2){ u64 Dm=((u64)e.b*(u64)GNC+19)/20; need=(Dm>e.s)?(Dm-e.s+7)/8:0; }
    else if(pol>=3){ long long v=2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                              +2LL*(long long)e.b-(long long)e.s
                              -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
                     need=(v>0)?(u64)((v+7)/8):0; }
    if(pol>=2){ p=min(min(need,e.w),x1/50); x1-=50*p; }
    k=min(min(e.b,e.s),x1/25); x1-=25*k;
    if(pol==3) q=min(e.w-p,x1/200);
}

// ---------------- LP ----------------
struct LPS {
    int Hz,N,M,W; vector<D> Tb;
    D X0,S0,B0,W0,AW[H],AM[H],AS[H],AB[H];
    inline D& at(int r,int c){ return Tb[(size_t)r*W+c]; }
    int Kx(int t){return t;} int Px(int t){return Hz+t;} int Qx(int t){return 2*Hz+t;}
    D G0(int tau){ D r=W0; for(int i=1;i<H&&i<=tau;i++) r+=AW[i]; return r; }
    D X0c(int tau){ D r=X0; for(int i=1;i<H&&i<=tau;i++) r+=AM[i]; return r; }
    D S0c(int tau){ D r=S0; for(int i=1;i<H&&i<=tau;i++) r+=AS[i]; return r; }
    D B0c(int tau){ D r=B0; for(int i=1;i<H&&i<=tau;i++) r+=AB[i]; return r; }
    vector<D> coef;
    D xExpr(int t){
        coef.assign(N,0); D c0=X0c(t);
        for(int u=0;u<Hz;u++){
            D ck=4.0L*max(0,t-u-2);
            int hp=min(t-1,u+2); D cp=(hp>=u)?(D)(hp-u+1):0;
            int hq=min(t-1,u+11); D cq=(hq>=u)?(D)(hq-u+1):0;
            coef[Kx(u)]+= ck - (u<t?25.0L:0.0L);
            coef[Px(u)]+= -4*cp - (u<t?50.0L:0.0L);
            coef[Qx(u)]+= -4*cq - (u<t?200.0L:0.0L);
        }
        for(int v=0;v<t;v++) c0+=4.0L*G0(v);
        return c0;
    }
    void setState(const St&e,D lam){
        X0=(D)e.x*lam;S0=(D)e.s*lam;W0=(D)e.w*lam;B0=(D)e.b*lam;
        for(int i=0;i<H;i++){AW[i]=(D)e.aw[i]*lam;AM[i]=(D)e.am[i]*lam;AS[i]=(D)e.as_[i]*lam;AB[i]=(D)e.ab[i]*lam;}
    }
    void build(int hz){
        Hz=hz; N=3*Hz; M=4*Hz; W=N+M+1;
        Tb.assign((size_t)(M+1)*W,0.0L);
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
    long long iters=0;
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
            basis[pr]=pc; iters++;
        }
        if(sol){ sol->assign(N,0); for(int r=0;r<M;r++) if(basis[r]<N) (*sol)[basis[r]]=at(r,N+M); }
        return at(M,N+M);
    }
};

static D lpBound(const St&e,int hz,vector<D>*sol,D*lamOut=0){
    LPS L; D lam=1.0L; if(e.x>2000000ULL) lam=1000000.0L/(D)e.x;
    L.setState(e,lam); L.build(hz); D v=L.solve(sol);
    if(lamOut)*lamOut=lam;
    return v/lam;   // 回到原尺度
}

// ---------- 带进位取整 + 精确 u64 重放 ----------
struct RunRes{ u64 val; vector<array<u64,3> > A; };
static RunRes carryRound(const St&e0,const vector<D>&sol,int Hz,D lam){
    St c=e0; D rk=0,rp=0,rq=0; RunRes R; R.A.resize(Hz);
    for(int t=0;t<Hz;t++){
        D dk=sol[t]/lam+rk, dp=sol[Hz+t]/lam+rp, dq=sol[2*Hz+t]/lam+rq;
        long long k=(long long)floorl(dk+0.5L), p=(long long)floorl(dp+0.5L), q=(long long)floorl(dq+0.5L);
        if(k<0)k=0; if(p<0)p=0; if(q<0)q=0;
        if((u64)k>c.b)k=(long long)c.b; if((u64)k>c.s)k=(long long)c.s;
        if((u64)(p+q)>c.w){ if((u64)p>c.w){p=(long long)c.w;q=0;} else q=(long long)c.w-p; }
        while((unsigned long long)(25*k+50*p+200*q)>c.x){ if(q>0)q--; else if(p>0)p--; else if(k>0)k--; else break; }
        rk=dk-(D)k; rp=dp-(D)p; rq=dq-(D)q;
        if(rk>4)rk=4; if(rk<-4)rk=-4; if(rp>4)rp=4; if(rp<-4)rp=-4; if(rq>4)rq=4; if(rq<-4)rq=-4;
        R.A[t]={(u64)k,(u64)p,(u64)q};
        c=step(c,(u64)k,(u64)p,(u64)q);
    }
    R.val=c.x; return R;
}
// 坐标下降修复：对每一步做小扰动（含矿物中性代换），尾部按记录动作夹紧重放
static u64 repair(const St&e0,vector<array<u64,3> >&A,int Hz,int passes){
    vector<St> TR(Hz+1);
    auto roll=[&](int from)->u64{
        for(int i=from;i<Hz;i++){ u64 k=A[i][0],p=A[i][1],q=A[i][2]; const St&c=TR[i];
            if(k>c.b)k=c.b; if(k>c.s)k=c.s;
            if(p+q>c.w){ if(p>c.w){p=c.w;q=0;} else q=c.w-p; }
            while(25*k+50*p+200*q>c.x){ if(q)q--; else if(p)p--; else if(k)k--; else break; }
            A[i][0]=k;A[i][1]=p;A[i][2]=q; TR[i+1]=step(c,k,p,q); }
        return TR[Hz].x; };
    TR[0]=e0; u64 cur=roll(0);
    for(int pass=0;pass<passes;pass++){ bool imp=false;
        for(int i=Hz-1;i>=0;i--){
            const St c=TR[i]; long long K=A[i][0],P=A[i][1],Q=A[i][2];
            static const long long DD[10]={1,-1,2,-2,4,-4,8,-8,16,-16};
            vector<array<long long,3> > tr;
            for(int a=0;a<10;a++){ long long d=DD[a];
                tr.push_back({K+d,P,Q}); tr.push_back({K,P+d,Q}); tr.push_back({K,P,Q+d});
                tr.push_back({K-8*d,P,Q+d});    // 矿物中性
                tr.push_back({K,P-4*d,Q+d});
                tr.push_back({K-2*d,P+d,Q}); }
            tr.push_back({(long long)min(min(c.b,c.s),c.x/25),P,Q});
            for(size_t z=0;z<tr.size();z++){
                long long k=tr[z][0],p=tr[z][1],q=tr[z][2];
                if(k<0||p<0||q<0) continue;
                if((u64)k>c.b||(u64)k>c.s||(u64)(p+q)>c.w) continue;
                if((unsigned long long)(25*k+50*p+200*q)>c.x) continue;
                vector<array<u64,3> > bak(A.begin()+i,A.end());
                A[i]={(u64)k,(u64)p,(u64)q};
                u64 v2=roll(i);
                if(v2>cur){ cur=v2; imp=true; }
                else { copy(bak.begin(),bak.end(),A.begin()+i); roll(i); }
            }
        }
        if(!imp) break; }
    return cur;
}
// 滚动时域：每步用当前整数状态重解 LP，取首步动作带进位取整
static RunRes recedeDrive(const St&e0,int Hz,int RESOLVE){
    St c=e0; RunRes R; R.A.resize(Hz); D rk=0,rp=0,rq=0;
    vector<D> sol; D lam=1.0L; int since=1<<30;
    for(int t=0;t<Hz;t++){
        int hz=Hz-t;
        if(since>=RESOLVE){
            LPS L; lam=1.0L; if(c.x>2000000ULL) lam=1000000.0L/(D)c.x;
            L.setState(c,lam); L.build(hz); L.solve(&sol); since=0;
            rk=rp=rq=0;
        }
        // sol 的下标基于重解时刻的窗口长度 hzr；本步是该窗口内的第 since 步
        int o=since, hzr=hz+since;
        D dk=sol[o]/lam+rk, dp=sol[hzr+o]/lam+rp, dq=sol[2*hzr+o]/lam+rq;
        long long k=(long long)floorl(dk+0.5L),p=(long long)floorl(dp+0.5L),q=(long long)floorl(dq+0.5L);
        if(k<0)k=0; if(p<0)p=0; if(q<0)q=0;
        if((u64)k>c.b)k=(long long)c.b; if((u64)k>c.s)k=(long long)c.s;
        if((u64)(p+q)>c.w){ if((u64)p>c.w){p=(long long)c.w;q=0;} else q=(long long)c.w-p; }
        while((unsigned long long)(25*k+50*p+200*q)>c.x){ if(q>0)q--; else if(p>0)p--; else if(k>0)k--; else break; }
        rk=dk-(D)k; rp=dp-(D)p; rq=dq-(D)q;
        if(rk>4)rk=4; if(rk<-4)rk=-4; if(rp>4)rp=4; if(rp<-4)rp=-4; if(rq>4)rq=4; if(rq<-4)rq=-4;
        R.A[t]={(u64)k,(u64)p,(u64)q};
        c=step(c,(u64)k,(u64)p,(u64)q); since++;
    }
    R.val=c.x; return R;
}
// LP 前瞻取整：每步解 LP 得分数首步动作，对 floor/ceil 组合逐个用「后继状态的 LP 值」打分取最优。
// 这样每一步的取整决策都由精确的价值函数（LP）裁决，而不是盲目四舍五入。
static u64 lpLookahead(const St&e0,int Hz,vector<array<u64,3> >*out){
    St c=e0; vector<array<u64,3> > A(Hz);
    for(int t=0;t<Hz;t++){
        int hz=Hz-t;
        vector<D> sol; D lam=1.0L; if(c.x>2000000ULL) lam=1000000.0L/(D)c.x;
        { LPS L; L.setState(c,lam); L.build(hz); L.solve(&sol); }
        D fk=sol[0]/lam, fp=sol[hz]/lam, fq=sol[2*hz]/lam;
        long long bk=(long long)floorl(fk), bp=(long long)floorl(fp), bq=(long long)floorl(fq);
        if(bk<0)bk=0; if(bp<0)bp=0; if(bq<0)bq=0;
        long long bestK=-1,bestP=0,bestQ=0; D bestV=-1;
        for(int dk=0;dk<=1;dk++) for(int dp=0;dp<=1;dp++) for(int dq=0;dq<=1;dq++){
            long long k=bk+dk,p=bp+dp,q=bq+dq;
            if((u64)k>c.b||(u64)k>c.s) continue;
            if((u64)(p+q)>c.w) continue;
            if((unsigned long long)(25*k+50*p+200*q)>c.x) continue;
            St d=step(c,(u64)k,(u64)p,(u64)q);
            D v;
            if(hz==1) v=(D)d.x;
            else { LPS L2; D l2=1.0L; if(d.x>2000000ULL) l2=1000000.0L/(D)d.x;
                   L2.setState(d,l2); L2.build(hz-1); v=L2.solve(0)/l2; }
            if(v>bestV){ bestV=v; bestK=k; bestP=p; bestQ=q; }
        }
        if(bestK<0){ bestK=bestP=bestQ=0; }   // 兜底：全采矿
        A[t]={(u64)bestK,(u64)bestP,(u64)bestQ};
        c=step(c,(u64)bestK,(u64)bestP,(u64)bestQ);
    }
    if(out)*out=A;
    return c.x;
}
// LP 束搜索：每步对每个束状态解 LP 取分数首步动作，在其邻域枚举整数候选，
// 用「后继状态到 T 的 LP 上界」给每个后继打分，保留前 B 个。
// LP 值是这一尺度上最准的价值函数（整数性缺口 ~1e-7），比任何固定策略族强得多。
struct BeamEnt{ St s; vector<array<u64,3> > A; };
static u64 lpBeam(const vector<St>&init,int Hz,int B,int DW,vector<array<u64,3> >*bestA,int*bestI,int verbose){
    vector<BeamEnt> cur(init.size());
    for(size_t i=0;i<init.size();i++){ cur[i].s=init[i]; }
    vector<int> origin(init.size()); for(size_t i=0;i<init.size();i++) origin[i]=i;
    u64 BEST=0; vector<array<u64,3> > BESTA; int BESTI=-1;
    for(int t=0;t<Hz;t++){
        int hz=Hz-t;
        int n=cur.size();
        vector<BeamEnt> nx; vector<int> nxo;
        vector<vector<array<u64,3> > > cand(n);
        #pragma omp parallel for schedule(dynamic,1)
        for(int i=0;i<n;i++){
            const St&c=cur[i].s;
            vector<D> sol; D lam=1.0L; if(c.x>2000000ULL) lam=1000000.0L/(D)c.x;
            { LPS L; L.setState(c,lam); L.build(hz); L.solve(&sol); }
            long long bk=(long long)floorl(sol[0]/lam), bp=(long long)floorl(sol[hz]/lam), bq=(long long)floorl(sol[2*hz]/lam);
            if(bk<0)bk=0; if(bp<0)bp=0; if(bq<0)bq=0;
            for(int dk=-DW;dk<=DW+1;dk++) for(int dp=-DW;dp<=DW+1;dp++) for(int dq=-DW;dq<=DW+1;dq++){
                long long k=bk+dk,p=bp+dp,q=bq+dq;
                if(k<0||p<0||q<0) continue;
                if((u64)k>c.b||(u64)k>c.s||(u64)(p+q)>c.w) continue;
                if((unsigned long long)(25*k+50*p+200*q)>c.x) continue;
                cand[i].push_back({(u64)k,(u64)p,(u64)q});
            }
        }
        for(int i=0;i<n;i++) for(size_t j=0;j<cand[i].size();j++){
            BeamEnt e; e.s=step(cur[i].s,cand[i][j][0],cand[i][j][1],cand[i][j][2]);
            e.A=cur[i].A; e.A.push_back(cand[i][j]);
            nx.push_back(e); nxo.push_back(origin[i]);
        }
        // 去重
        { vector<int> ord(nx.size()); iota(ord.begin(),ord.end(),0);
          sort(ord.begin(),ord.end(),[&](int a,int b){return memcmp(&nx[a].s,&nx[b].s,sizeof(St))<0;});
          vector<BeamEnt> n2; vector<int> o2;
          for(size_t z=0;z<ord.size();z++){
            if(!n2.empty()&&memcmp(&n2.back().s,&nx[ord[z]].s,sizeof(St))==0) continue;
            n2.push_back(nx[ord[z]]); o2.push_back(nxo[ord[z]]); }
          nx.swap(n2); nxo.swap(o2); }
        int m=nx.size();
        vector<D> sc(m);
        if(hz-1<=0){ for(int i=0;i<m;i++) sc[i]=(D)nx[i].s.x; }
        else {
            #pragma omp parallel for schedule(dynamic,1)
            for(int i=0;i<m;i++){ LPS L; D l=1.0L; if(nx[i].s.x>2000000ULL) l=1000000.0L/(D)nx[i].s.x;
                L.setState(nx[i].s,l); L.build(hz-1); sc[i]=L.solve(0)/l; }
        }
        for(int i=0;i<m;i++) if(nx[i].s.x>BEST){BEST=nx[i].s.x;BESTA=nx[i].A;BESTI=nxo[i];}
        vector<int> v(m); iota(v.begin(),v.end(),0);
        int kk=min(B,m);
        partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return sc[a]>sc[b];});
        vector<BeamEnt> ns(kk); vector<int> no(kk);
        for(int i=0;i<kk;i++){ ns[i]=nx[v[i]]; no[i]=nxo[v[i]]; }
        cur.swap(ns); origin.swap(no);
        if(verbose&&(t%10==9||t==Hz-1)) fprintf(stderr,"  LP束 t=+%d 候选=%d 束=%d 当前最好*2=%llu\n",t+1,m,(int)cur.size(),2*BEST);
    }
    if(bestA)*bestA=BESTA; if(bestI)*bestI=BESTI;
    return BEST;
}
// 三段固定策略族网格（从状态 e 出发，hz 步到 T）
static u64 gridCash(const St&e,int hz,int CCM,int CAM,int CBM){
    u64 best=e.x;
    for(int cc=0;cc<=CCM&&cc<=hz;cc++) for(int a=0;a+cc<=hz&&a<=CAM;a++) for(int b2=0;b2+a+cc<=hz&&b2<=CBM;b2++){
        St d=e; bool good=true;
        for(int u=0;u<hz;u++){ int pol=(u<cc)?4:((u<cc+a)?2:((u<cc+a+b2)?1:0));
            u64 k,p,q; act(d,pol,k,p,q); if(!okA(d,k,p,q)){good=false;break;}
            d=step(d,k,p,q); }
        if(good&&d.x>best)best=d.x;
    }
    return best;
}
// 精确 u64 重放并逐步断言四条约束
static bool verify(const St&e0,const vector<array<u64,3> >&A,u64 expect,int&badstep){
    St c=e0;
    for(size_t i=0;i<A.size();i++){
        u64 k=A[i][0],p=A[i][1],q=A[i][2];
        if(!(p+q<=c.w)){badstep=i;return false;}
        if(!(k<=c.b)){badstep=i;return false;}
        if(!(k<=c.s)){badstep=i;return false;}
        if(!(25*k+50*p+200*q<=c.x)){badstep=i;return false;}
        c=step(c,k,p,q);
    }
    badstep=-1; return c.x==expect;
}

// 成长期滚动 LP 值束搜索：每步用「到 +HZ 的 LP 上界」给后继打分，保留前 B 个。
// 逐时刻把整个束写出，供随后的兑现段使用。
int main(int argc,char**argv){
    if(argc<7){printf("用法: growlp <beam.bin> <TEND> <HZ> <B> <K> <outdir> [DELW] [DUMPEVERY]\n");return 1;}
    const char*fn=argv[1]; int TEND=atoi(argv[2]), HZ=atoi(argv[3]);
    int B=atoi(argv[4]), K=atoi(argv[5]); const char*od=argv[6];
    int DW=argc>7?atoi(argv[7]):1; int DE=argc>8?atoi(argv[8]):1;
    FILE*f=fopen(fn,"rb"); if(!f){printf("no %s\n",fn);return 1;}
    int cnt,t0; if(fread(&cnt,4,1,f)!=1||fread(&t0,4,1,f)!=1)return 1;
    vector<St> all(cnt);
    for(int i=0;i<cnt;i++) if(fread(&all[i],sizeof(St),1,f)!=1){cnt=i;all.resize(cnt);break;}
    fclose(f);
    fprintf(stderr,"%s: %d 状态 t0=%d -> TEND=%d, 视野 HZ=%d, 束 %d\n",fn,cnt,t0,TEND,HZ,B);
    vector<pair<D,int> > rk(cnt);
    #pragma omp parallel for schedule(dynamic,4)
    for(int i=0;i<cnt;i++){ D v=lpBound(all[i],HZ,0); rk[i]=make_pair(v,i); }
    sort(rk.begin(),rk.end(),[](const pair<D,int>&a,const pair<D,int>&b){return a.first>b.first;});
    int lim=min(K,cnt);
    vector<St> cur; for(int i=0;i<lim;i++) cur.push_back(all[rk[i].second]);
    fprintf(stderr,"起始束 %d 个, 最好 LP*2=%.1Lf\n",lim,2*rk[0].first);
    double w0=omp_get_wtime();
    auto dump=[&](int t){ char p[512]; snprintf(p,sizeof(p),"%s/g%d.bin",od,t);
        FILE*fo=fopen(p,"wb"); if(!fo){fprintf(stderr,"无法写 %s\n",p);return;}
        int c2=cur.size(); fwrite(&c2,4,1,fo); fwrite(&t,4,1,fo);
        for(int i=0;i<c2;i++) fwrite(&cur[i],sizeof(St),1,fo); fclose(fo); };
    if(t0%DE==0) dump(t0);
    for(int t=t0;t<TEND;t++){
        int n=cur.size();
        vector<vector<array<u64,3> > > cand(n);
        #pragma omp parallel for schedule(dynamic,1)
        for(int i=0;i<n;i++){
            const St&c=cur[i];
            vector<D> sol; D lam=1.0L; if(c.x>2000000ULL) lam=1000000.0L/(D)c.x;
            { LPS L; L.setState(c,lam); L.build(HZ); L.solve(&sol); }
            long long bk=(long long)floorl(sol[0]/lam), bp=(long long)floorl(sol[HZ]/lam), bq=(long long)floorl(sol[2*HZ]/lam);
            if(bk<0)bk=0; if(bp<0)bp=0; if(bq<0)bq=0;
            for(int dk=-DW;dk<=DW+1;dk++) for(int dp=-DW;dp<=DW+1;dp++) for(int dq=-DW;dq<=DW+1;dq++){
                long long k=bk+dk,p=bp+dp,q=bq+dq;
                if(k<0||p<0||q<0) continue;
                if((u64)k>c.b||(u64)k>c.s||(u64)(p+q)>c.w) continue;
                if((unsigned long long)(25*k+50*p+200*q)>c.x) continue;
                cand[i].push_back({(u64)k,(u64)p,(u64)q});
            }
        }
        vector<St> nx;
        for(int i=0;i<n;i++) for(size_t j=0;j<cand[i].size();j++)
            nx.push_back(step(cur[i],cand[i][j][0],cand[i][j][1],cand[i][j][2]));
        sort(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))<0;});
        nx.erase(unique(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))==0;}),nx.end());
        int m=nx.size(); if(m==0){fprintf(stderr,"t=%d 无后继\n",t);break;}
        vector<D> sc(m);
        #pragma omp parallel for schedule(dynamic,1)
        for(int i=0;i<m;i++) sc[i]=lpBound(nx[i],HZ,0);
        vector<int> v(m); iota(v.begin(),v.end(),0);
        int kk=min(B,m);
        partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return sc[a]>sc[b];});
        vector<St> ns(kk); for(int i=0;i<kk;i++) ns[i]=nx[v[i]];
        cur.swap(ns);
        if((t+1)%DE==0||t+1==TEND) dump(t+1);
        if((t+1)%10==0) fprintf(stderr,"t=%d 候选=%d 束=%d 最好x*2=%llu 用时 %.0fs\n",
            t+1,m,(int)cur.size(),2*cur[0].x,omp_get_wtime()-w0);
    }
    return 0;
}
