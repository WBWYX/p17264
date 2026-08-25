// 从任意实际状态出发的精确 LP：给定 t0 的状态与到达流水线，max x(T)
#include <bits/stdc++.h>
using namespace std;
typedef long double D; typedef unsigned long long u64;
static const int H=13;
struct St{u64 x,s,w,b;u64 aw[H],am[H],as_[H],ab[H];};
static St stepSim(const St&t,u64 k,u64 p,u64 q){
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
int Hz,N,M;                       // Hz = 窗口步数
vector<D> Tb;
inline D& at(int r,int c){ return Tb[(size_t)r*(N+M+1)+c]; }
int Kx(int t){return t;} int Px(int t){return Hz+t;} int Qx(int t){return 2*Hz+t;}
D X0,S0,B0,W0; D AW[H],AM[H],AS[H],AB[H];
static D G0(int tau){ D r=W0; for(int i=1;i<H&&i<=tau;i++) r+=AW[i]; return r; }
static D X0c(int tau){ D r=X0; for(int i=1;i<H&&i<=tau;i++) r+=AM[i]; return r; }
static D S0c(int tau){ D r=S0; for(int i=1;i<H&&i<=tau;i++) r+=AS[i]; return r; }
static D B0c(int tau){ D r=B0; for(int i=1;i<H&&i<=tau;i++) r+=AB[i]; return r; }
// x(tau) 线性表达
static D xExpr(int t, vector<D>&coef){
    coef.assign(N,0); D c0=X0c(t);
    for(int u=0;u<Hz;u++){
        // 4*sum_{v<=t-1} g(v)， g(v)=G0(v)+sum_{u<=v-2}k - (p(v)+p(v-1)+p(v-2)) - sum_{j=0..11}q(v-j)
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
int main(int argc,char**argv){
    const char*fn=argc>1?argv[1]:"beam.bin";
    int T=argc>2?atoi(argv[2]):171;
    int TOPN=argc>3?atoi(argv[3]):40;
    FILE*f=fopen(fn,"rb"); if(!f){printf("no %s\n",fn);return 1;}
    int cnt,t0; if(fread(&cnt,4,1,f)!=1||fread(&t0,4,1,f)!=1)return 1;
    vector<St> all(cnt);
    for(int i=0;i<cnt;i++) if(fread(&all[i],sizeof(St),1,f)!=1){cnt=i;all.resize(cnt);break;}
    fclose(f);
    Hz=T-t0; if(Hz<=0){printf("T<=t0\n");return 1;}
    N=3*Hz; M=4*Hz;
    // 先按当前矿物量排序取前 TOPN（也把矿物最多的和"总工人最多的"都覆盖）
    vector<int> idx(cnt); iota(idx.begin(),idx.end(),0);
    auto tot=[&](const St&e){ D r=e.w; for(int i=1;i<H;i++)r+=e.aw[i]; return r; };
    sort(idx.begin(),idx.end(),[&](int a,int b){ return tot(all[a])>tot(all[b]); });
    int lim=min(cnt,TOPN);
    printf("状态数=%d  t0=%d  T=%d  窗口=%d  评估前 %d 个(按总工人数)\n",cnt,t0,T,Hz,lim);
    D best=0; int bi=-1; vector<pair<double,int> > LPV;
    Tb.assign((size_t)(M+1)*(N+M+1),0.0L);
    vector<D> coef;
    for(int ii=0;ii<lim;ii++){
        const St&e=all[idx[ii]];
        X0=e.x;S0=e.s;W0=e.w;B0=e.b;
        for(int i=0;i<H;i++){AW[i]=e.aw[i];AM[i]=e.am[i];AS[i]=e.as_[i];AB[i]=e.ab[i];}
        fill(Tb.begin(),Tb.end(),0.0L);
        for(int t=0;t<Hz;t++){
            { int r=4*t; D c0=xExpr(t,coef);
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
        { D c0=xExpr(Hz,coef); for(int j=0;j<N;j++) at(M,j)=-coef[j]; at(M,N+M)=c0; }
        vector<int> basis(M); for(int r=0;r<M;r++) basis[r]=N+r;
        int W=N+M+1;
        for(int it=0;it<200000;it++){
            int pc=-1; D bs=-1e-9L;
            for(int j=0;j<N+M;j++) if(at(M,j)<bs){bs=at(M,j);pc=j;}
            if(pc<0) break;
            int pr=-1; D bR=0;
            for(int r=0;r<M;r++){ D a=at(r,pc); if(a>1e-11L){ D rt=at(r,N+M)/a;
                if(pr<0||rt<bR-1e-12L||(fabsl(rt-bR)<=1e-12L&&basis[r]<basis[pr])){pr=r;bR=rt;} } }
            if(pr<0) break;
            D pv=at(pr,pc); for(int c=0;c<W;c++) at(pr,c)/=pv;
            for(int r=0;r<=M;r++) if(r!=pr){ D fq=at(r,pc); if(fq!=0){ D*R=&at(r,0),*P=&at(pr,0);
                for(int c=0;c<W;c++) R[c]-=fq*P[c]; } }
            basis[pr]=pc;
        }
        D v=at(M,N+M); LPV.push_back(make_pair((double)v,idx[ii]));
        if(v>best){best=v;bi=idx[ii];
            if(getenv("LPDUMP")){ vector<D> sol(N,0);
                for(int r=0;r<M;r++) if(basis[r]<N) sol[basis[r]]=at(r,N+M);
                printf("  tau      k         p         q\n");
                for(int t=0;t<Hz;t++) printf("  %3d  %9.2Lf %9.4Lf %9.4Lf\n",t,sol[Kx(t)],sol[Px(t)],sol[Qx(t)]);
            } }
    }
    printf("最优 LP 上界 x(%d) = %.3Lf   *2 = %.3Lf   (来自状态 #%d)\n",T,best,2*best,bi);
    if(getenv("LPROUND")&&bi>=0){
        // 重解最优状态的 LP，取原始解并取整模拟
        const St&e=all[bi];
        X0=e.x;S0=e.s;W0=e.w;B0=e.b;
        for(int i=0;i<H;i++){AW[i]=e.aw[i];AM[i]=e.am[i];AS[i]=e.as_[i];AB[i]=e.ab[i];}
        fill(Tb.begin(),Tb.end(),0.0L);
        vector<D> cf;
        for(int t=0;t<Hz;t++){
            { int r=4*t; D c0=xExpr(t,cf);
              for(int j=0;j<N;j++) at(r,j)=-cf[j];
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
        { D c0=xExpr(Hz,cf); for(int j=0;j<N;j++) at(M,j)=-cf[j]; at(M,N+M)=c0; }
        vector<int> bs2(M); for(int r=0;r<M;r++) bs2[r]=N+r;
        int W2=N+M+1;
        for(int it=0;it<200000;it++){
            int pc=-1; D b2=-1e-9L;
            for(int j=0;j<N+M;j++) if(at(M,j)<b2){b2=at(M,j);pc=j;}
            if(pc<0) break;
            int pr=-1; D bR=0;
            for(int r=0;r<M;r++){ D a=at(r,pc); if(a>1e-11L){ D rt=at(r,N+M)/a;
                if(pr<0||rt<bR-1e-12L||(fabsl(rt-bR)<=1e-12L&&bs2[r]<bs2[pr])){pr=r;bR=rt;} } }
            if(pr<0) break;
            D pv=at(pr,pc); for(int c=0;c<W2;c++) at(pr,c)/=pv;
            for(int r=0;r<=M;r++) if(r!=pr){ D fq=at(r,pc); if(fq!=0){ D*R=&at(r,0),*P=&at(pr,0);
                for(int c=0;c<W2;c++) R[c]-=fq*P[c]; } }
            bs2[pr]=pc;
        }
        vector<D> sol(N,0);
        for(int r=0;r<M;r++) if(bs2[r]<N) sol[bs2[r]]=at(r,N+M);
        // 取整模拟：floor / round 两种
        for(int mode=0;mode<2;mode++){
            St d=e; bool ok=true;
            for(int t=0;t<Hz;t++){
                long long k=(long long)(mode? floorl(sol[Kx(t)]+0.5L) : floorl(sol[Kx(t)]));
                long long p=(long long)(mode? floorl(sol[Px(t)]+0.5L) : floorl(sol[Px(t)]));
                long long q=(long long)(mode? floorl(sol[Qx(t)]+0.5L) : floorl(sol[Qx(t)]));
                if(k<0)k=0; if(p<0)p=0; if(q<0)q=0;
                if((u64)k>d.b)k=d.b; if((u64)k>d.s)k=d.s;
                if((u64)(p+q)>d.w){ q=(long long)d.w-p; if(q<0){q=0;p=(long long)d.w;} }
                while(25*k+50*p+200*q>(long long)d.x){ if(q>0)q--; else if(p>0)p--; else if(k>0)k--; else {ok=false;break;} }
                if(!ok)break;
                d=stepSim(d,(u64)k,(u64)p,(u64)q);
            }
            if(ok){
                printf("LP解取整(%s) 模拟得 x(%d)*2 = %llu\n",mode?"round":"floor",T,2*d.x);
                // ---- 坐标下降修复 ----
                vector<array<u64,3> > A(Hz); vector<St> TR(Hz+1);
                { St c=e;
                  for(int t=0;t<Hz;t++){
                    long long k=(long long)(mode? floorl(sol[Kx(t)]+0.5L) : floorl(sol[Kx(t)]));
                    long long p=(long long)(mode? floorl(sol[Px(t)]+0.5L) : floorl(sol[Px(t)]));
                    long long q=(long long)(mode? floorl(sol[Qx(t)]+0.5L) : floorl(sol[Qx(t)]));
                    if(k<0)k=0; if(p<0)p=0; if(q<0)q=0;
                    if((u64)k>c.b)k=c.b; if((u64)k>c.s)k=c.s;
                    if((u64)(p+q)>c.w){ q=(long long)c.w-p; if(q<0){q=0;p=(long long)c.w;} }
                    while(25*k+50*p+200*q>(long long)c.x){ if(q>0)q--; else if(p>0)p--; else if(k>0)k--; else break; }
                    A[t]={(u64)k,(u64)p,(u64)q}; TR[t]=c; c=stepSim(c,(u64)k,(u64)p,(u64)q); }
                  TR[Hz]=c; }
                auto roll=[&](int from)->u64{
                    for(int i=from;i<Hz;i++){ u64 k=A[i][0],p=A[i][1],q=A[i][2]; const St&c=TR[i];
                        if(k>c.b)k=c.b; if(k>c.s)k=c.s;
                        if(p+q>c.w){ if(p>c.w){p=c.w;q=0;} else q=c.w-p; }
                        while(25*k+50*p+200*q>c.x){ if(q)q--; else if(p)p--; else if(k)k--; else break; }
                        A[i][0]=k;A[i][1]=p;A[i][2]=q; TR[i+1]=stepSim(c,k,p,q); }
                    return TR[Hz].x; };
                u64 cur=roll(0);
                for(int pass=0;pass<40;pass++){ bool imp=false;
                    for(int i=Hz-1;i>=0;i--){
                        const St&c=TR[i]; long long K=A[i][0],P=A[i][1],Q=A[i][2];
                        long long D1[9]={1,-1,2,-2,4,-4,8,-8,16};
                        vector<array<long long,3> > tr;
                        for(int a=0;a<9;a++){ tr.push_back({K+D1[a],P,Q}); tr.push_back({K,P,Q+D1[a]});
                                              tr.push_back({K,P+D1[a],Q}); tr.push_back({K+D1[a],P,Q-D1[a]}); }
                        tr.push_back({(long long)min(min(c.b,c.s),c.x/25),P,Q});
                        for(size_t z=0;z<tr.size();z++){
                            long long k=tr[z][0],p=tr[z][1],q=tr[z][2];
                            if(k<0||p<0||q<0) continue;
                            if((u64)k>c.b||(u64)k>c.s||(u64)(p+q)>c.w) continue;
                            if(25*k+50*p+200*q>(long long)c.x) continue;
                            vector<array<u64,3> > bak(A.begin()+i,A.end());
                            A[i]={(u64)k,(u64)p,(u64)q};
                            u64 v2=roll(i);
                            if(v2>cur){ cur=v2; imp=true; }
                            else { copy(bak.begin(),bak.end(),A.begin()+i); roll(i); }
                        }
                    }
                    if(!imp) break; }
                printf("  -> 坐标下降后 x(%d)*2 = %llu\n",T,2*cur);
            }
        }
    }
    if(getenv("LPTOP")){ int K=atoi(getenv("LPTOP"));
        sort(LPV.begin(),LPV.end(),[](const pair<double,int>&a,const pair<double,int>&b){return a.first>b.first;});
        int kk=min((int)LPV.size(),K);
        FILE*fo=fopen("top.bin","wb"); fwrite(&kk,4,1,fo); fwrite(&t0,4,1,fo);
        for(int i=0;i<kk;i++) fwrite(&all[LPV[i].second],sizeof(St),1,fo); fclose(fo);
        printf("按 LP 值写出前 %d 个状态到 top.bin (最好 %.1f, 第%d好 %.1f)\n",kk,2*LPV[0].first,kk,2*LPV[kk-1].first); }
    return 0;
}
