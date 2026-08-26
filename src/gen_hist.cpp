// P17264 打表器：束搜索优化开局 + ORDER=2 纯策略接管成长期
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64;
typedef long double ld;
static const u64 CAP=(u64)1<<63;
static const int H=13;

struct St{u64 x,s,w,b;u64 aw[H],am[H],as_[H],ab[H];};
static inline u64 addc(u64 a,u64 b){u64 c=a+b;return (c<a||c>CAP)?CAP:c;}

static St step(const St&t,u64 k,u64 p,u64 q){
    St c=t; c.x-=25*k+50*p+200*q; c.s-=k; c.b-=k;
    c.ab[2]+=k;c.aw[2]+=k; c.aw[3]+=p;c.as_[3]+=8*p;
    c.aw[12]+=q;c.as_[12]+=10*q;c.ab[12]+=q;
    u64 g=t.w-p-q; c.aw[1]+=g; c.am[1]=addc(c.am[1],4*g);
    u64 nx=addc(c.x,c.am[1]),ns=c.s+c.as_[1],nw=c.aw[1],nb=c.b+c.ab[1];
    for(int i=1;i+1<H;i++){c.aw[i]=c.aw[i+1];c.am[i]=c.am[i+1];c.as_[i]=c.as_[i+1];c.ab[i]=c.ab[i+1];}
    c.aw[H-1]=c.am[H-1]=c.as_[H-1]=c.ab[H-1]=0;
    c.aw[0]=c.am[0]=c.as_[0]=c.ab[0]=0;
    c.x=nx;c.s=ns;c.w=nw;c.b=nb; return c;
}

// ---------- 参考策略 ORDER=2：补给站 -> 工人 -> 基地 ----------
u64 GNA[2]={20,21}; u64 GD=20;      // 成长期：周期 2 的补给目标倍数
u64 GNC=60, KAC=64, QAC=64;
u64 HGN[2]={20,21}; u64 HGNC=98, HKAC=64; int INHAND=0;   // 交接阶段专用的成长期补给参数         // 兑现期：补给目标倍数 / 造工人比例 / 建基地比例(/64)
int TSPLIT=200;                     // 该时刻之前兑现仍用成长期参数
// pol 0=全采矿 1=只造工人 2=补给+工人 3=完全成长 ; ph=时刻奇偶
// pol<=2 视为兑现期，用兑现期参数
int SUPMODE=2; int TDUMP=-1; int CSTEP=2,CHOR=48,CAMAX=30,CBMAX=18,CCMAX=0,CASHTOP=0;   // 0=旧 GNA 反馈  1=精确前馈  2=精确前馈+兜底
static void act(const St&e,int pol,u64&k,u64&p,u64&q,int ph){
    k=p=q=0;
    if(pol==0) return;
    bool cash = (pol<=2) && (ph>=TSPLIT);
    u64 gn = cash?(INHAND?HGNC:GNC):(INHAND?HGN[ph&1]:GNA[ph&1]);
    u64 x1=e.x;
    if(pol>=2){ u64 need;
        // SUPMODE: 0=旧规则 5=自校正(仅交接成长) 6=自校正(所有 pol==3) 7=自校正(pol>=2 且非兑现期)
        bool useff=false;
        if(SUPMODE==5) useff=(pol==3)&&INHAND;
        else if(SUPMODE==6) useff=(pol==3);
        else if(SUPMODE==7) useff=!cash;
        if(SUPMODE==8||(SUPMODE==9&&(!cash||INHAND))){   // 统一自校正： beta=1 成长 / beta=GNC/20 兑现
            long long B = cash ? (long long)(INHAND?HGNC:GNC) : 20LL;
            long long v = B*((long long)e.ab[1]+(long long)e.ab[3])/20
                        + 2LL*(long long)e.b + (long long)e.ab[1] + (long long)e.ab[2]
                        - (long long)e.s
                        - (long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
            need=(v>0)?(u64)((v+7)/8):0;
        }
        else if(!useff){ u64 D=(e.b*gn+GD-1)/GD; need=(D>e.s)?(D-e.s+7)/8:0; }   // 旧反馈规则
        else{   // 自校正： 8p >= 2ab1+ab2+ab3+2b-s-(as1+as2+as3)
            long long v = 2LL*(long long)e.ab[1] + (long long)e.ab[2] + (long long)e.ab[3]
                        + 2LL*(long long)e.b - (long long)e.s
                        - (long long)e.as_[1] - (long long)e.as_[2] - (long long)e.as_[3];
            need = (v>0)? (u64)((v+7)/8) : 0; }
        p=min(min(need,e.w),x1/50); x1-=50*p; }
    k=min(min(e.b,e.s),x1/25);
    if(cash) k=k*(INHAND?HKAC:KAC)/64;
    x1-=25*k;
    if(pol<=2||pol==4) return;
    q=min(e.w-p,x1/200);
}

// ---------- 线性估值（仅用于最粗一级筛选） ----------
ld RHO,eW,eB,eS; static const int HH=61; ld vW[HH],vB[HH];
static void calcE(){
    ld lo=1.0000001L,hi=1.2L;
    for(int i=0;i<300;i++){ld r=(lo+hi)/2,wW=4/(r-1),i2=powl(r,-2),i12=powl(r,-12);
        ld wB=(i2*wW-25)/(1-i2); if(i12*wB-200-wW*(1-i12)>0)lo=r;else hi=r;}
    RHO=(lo+hi)/2; eW=4/(RHO-1); eB=(powl(RHO,-2)*eW-25)/(1-powl(RHO,-2));
    eS=powl(RHO,3)*(eW*(1-powl(RHO,-3))+50)/8;
    for(int h=0;h<HH;h++){ vW[h]=vB[h]=0;
        if(h>=1)vW[h]=vW[h-1]+4;
        if(h>=12)vW[h]=max(vW[h],vW[h-12]+vB[h-12]-200);
        if(h>=1)vB[h]=vB[h-1];
        if(h>=2)vB[h]=max(vB[h],vB[h-2]+vW[h-2]-25); }
}
static ld sH(const St&t,int h,ld sw){
    ld r=(ld)t.x+vW[h]*t.w+vB[h]*t.b+sw*t.s;
    for(int i=1;i<H&&i<=h;i++) r+=(ld)t.am[i]+vW[h-i]*t.aw[i]+vB[h-i]*t.ab[i]+sw*t.as_[i];
    return r;
}
static ld sInf(const St&t){
    ld r=(ld)t.x+eW*t.w+eB*t.b+eS*t.s,d=1;
    for(int i=1;i<H;i++){d/=RHO; r+=d*((ld)t.am[i]+eW*t.aw[i]+eB*t.ab[i]+eS*t.as_[i]);}
    return r;
}

extern int SUPMODE;
// ---------- 长时域评分：double 松弛 + 环形缓冲，ORDER=2 策略 ----------
int NLONG=200, NCASH=22;
static double longVal(const St&e,int ph0){
    double x=e.x,ss=e.s,w=e.w,b=e.b;
    double AW[H],AM[H],AS[H],AB[H];
    for(int i=0;i<H;i++){AW[i]=e.aw[i];AM[i]=e.am[i];AS[i]=e.as_[i];AB[i]=e.ab[i];}
    int h=0; double sc=0;
    for(int j=0;j<NLONG+NCASH;j++){
        bool grow=j<NLONG;
        double wf=w,p=0,k=0,q=0;
        if(grow||j<NLONG+8){
            double need;
            if(SUPMODE>=4){ double a1=AB[(h+1)%H], a3=AB[(h+3)%H]; need=(a1-9*a3)/8;
                            if(need<0)need=0; if(ss<b&&(b-ss)/8>need) need=(b-ss)/8; }
            else { double D=b*((double)GNA[(ph0+j)&1]/GD); need=(D>ss)?(D-ss)/8:0; }
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

// ---------- 候选动作 ----------
static void gen(const St&e, vector<array<u64,3> >&out){
    static vector<u64> ks,qs,ps; static vector<array<u64,3> > tmp;
    u64 kM=min(min(e.b,e.s),e.x/25);
    ks.clear();
    if(kM<=20){ for(u64 k=0;k<=kM;k++) ks.push_back(k); }
    else { for(u64 j=0;j<=8;j++) ks.push_back(kM*j/8); for(u64 d=1;d<=3;d++) ks.push_back(kM-d); }
    sort(ks.begin(),ks.end()); ks.erase(unique(ks.begin(),ks.end()),ks.end());
    tmp.clear();
    u64 muls[7]={1,2,3,4,6,8,12};
    for(size_t ki=0;ki<ks.size();ki++){
        u64 k=ks[ki], x1=e.x-25*k;
        ps.clear(); ps.push_back(0);
        for(int mi=0;mi<7;mi++){
            u64 D=e.b*muls[mi], need=(D>e.s)?(D-e.s+7)/8:0;
            u64 p=min(min(need,e.w),x1/50); if(p) ps.push_back(p);
        }
        // 关键补充：自校正补给量（维持 s=b 的精确前馈值），旧的 D=mul*b 族无法表达
        { long long v = 2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                      + 2LL*(long long)e.b-(long long)e.s
                      -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
          u64 need=(v>0)?(u64)((v+7)/8):0;
          u64 p0=min(min(need,e.w),x1/50);
          for(int d=-2;d<=2;d++){ long long pv=(long long)p0+d;
              if(pv>0&&pv<=(long long)min(e.w,x1/50)) ps.push_back((u64)pv); }
          // 兑现口径（beta=GNC/20）的补给量也加入
          long long v2 = (long long)GNC*((long long)e.ab[1]+(long long)e.ab[3])/20
                       + 2LL*(long long)e.b + (long long)e.ab[1] + (long long)e.ab[2]
                       - (long long)e.s
                       -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
          u64 n2=(v2>0)?(u64)((v2+7)/8):0;
          u64 p2=min(min(n2,e.w),x1/50); if(p2) ps.push_back(p2);
          if(p2>1) ps.push_back(p2/2); }
        sort(ps.begin(),ps.end()); ps.erase(unique(ps.begin(),ps.end()),ps.end());
        for(size_t pi=0;pi<ps.size();pi++){
            u64 p=ps[pi], x2=x1-50*p, qM=min(e.w-p,x2/200);
            qs.clear();
            if(qM<=20){ for(u64 q=0;q<=qM;q++) qs.push_back(q); }
            else { for(u64 j=0;j<=8;j++) qs.push_back(qM*j/8); for(u64 d=1;d<=3;d++) qs.push_back(qM-d);
                   qs.push_back(1); qs.push_back(2); }
            sort(qs.begin(),qs.end()); qs.erase(unique(qs.begin(),qs.end()),qs.end());
            for(size_t qi=0;qi<qs.size();qi++){ array<u64,3> a={k,p,qs[qi]}; tmp.push_back(a); }
        }
    }
    out.clear();
    for(size_t i=0;i<tmp.size();i++){ array<u64,3> c=tmp[i];
        if(c[1]+c[2]>e.w)continue; if(25*c[0]+50*c[1]+200*c[2]>e.x)continue;
        if(c[0]>e.b||c[0]>e.s)continue; out.push_back(c); }
}

int TMAX=830, TBEAM=260, RH=64, P0=800, P1=200, BS=20, BL=400, LMAX=3000, NKEEP=6, TKEEP=260;
vector<u64> F;
vector<St> KEEP; vector<int> KEEPT; vector<St> HD;
static const int NBR=6; static int BRJ[NBR]={12,20,32,44,56,64};
static int HSEL[15]={0,2,4,6,9,12,16,20,25,30,36,42,48,56,64};
struct RV{ u64 v[4][65]; u64 br[NBR]; };

static void rollout(const St&e,int t,RV&r){
    for(int pol=0;pol<4;pol++){
        St c=e;
        for(int j=0;j<=RH;j++){
            r.v[pol][j]=c.x;
            if(t+j<=TMAX && c.x>F[t+j]) F[t+j]=c.x;
            if(j==RH) break;
            u64 k,p,q; act(c,pol,k,p,q,t+j); c=step(c,k,p,q);
        }
    }
    St c=e; int bi=0;
    for(int j=0;j<=RH&&bi<NBR;j++){
        while(bi<NBR&&BRJ[bi]==j){
            u64 best=c.x;
            // 三段兑现的代表形状： cc 步停基地(补给最小) -> aa 步猛补给+工人 -> bb 步只工人 -> 采矿
            static const int SH[4][3]={{0,8,6},{4,10,4},{8,12,3},{12,14,2}};
            for(int si=0;si<4;si++){
                int cc=SH[si][0],aa=SH[si][1],b2=SH[si][2];
                St d=c;
                for(int u=0;u<34;u++){
                    int pol=(u<cc)?4:((u<cc+aa)?2:((u<cc+aa+b2)?1:0));
                    u64 k,p,q; act(d,pol,k,p,q,t+j+u); d=step(d,k,p,q);
                    if(d.x>best)best=d.x;
                    if(t+j+u+1<=TMAX && d.x>F[t+j+u+1]) F[t+j+u+1]=d.x; }
            }
            r.br[bi]=best; bi++;
        }
        if(j<RH){ u64 k,p,q; act(c,3,k,p,q,t+j); c=step(c,k,p,q); }
    }
    while(bi<NBR){ r.br[bi]=c.x; bi++; }
}

int main(int argc,char**argv){
    if(argc>1)TMAX=atoi(argv[1]);
    if(argc>2)P1=atoi(argv[2]);
    if(argc>3)BS=atoi(argv[3]);
    if(argc>4)P0=atoi(argv[4]);
    if(argc>5)BL=atoi(argv[5]);
    if(argc>6)NLONG=atoi(argv[6]);
    if(argc>7)GNA[1]=atoi(argv[7]);
    if(argc>8)TKEEP=atoi(argv[8]); if(argc>9)TBEAM=atoi(argv[9]); if(argc>10)NKEEP=atoi(argv[10]); if(argc>11)GNC=atoi(argv[11]); if(argc>12)KAC=atoi(argv[12]); if(argc>13)TSPLIT=atoi(argv[13]); if(argc>15)GNA[0]=atoi(argv[14]),GNA[1]=atoi(argv[15]); if(argc>17)HGN[0]=atoi(argv[16]),HGN[1]=atoi(argv[17]); if(argc>18)HGNC=atoi(argv[18]); if(argc>19)HKAC=atoi(argv[19]); if(argc>20)SUPMODE=atoi(argv[20]); if(argc>21)TDUMP=atoi(argv[21]); if(argc>22)CSTEP=atoi(argv[22]); if(argc>23)CHOR=atoi(argv[23]); if(argc>24)CAMAX=atoi(argv[24]); if(argc>25)CBMAX=atoi(argv[25]); if(argc>26)CCMAX=atoi(argv[26]); if(argc>27)CASHTOP=atoi(argv[27]);
    calcE(); F.assign(TMAX+RH+40,0);
    vector<St> cur(1); memset(&cur[0],0,sizeof(St));
    cur[0].x=25;cur[0].s=6;cur[0].w=4;cur[0].b=1;
    { RV r0; rollout(cur[0],0,r0); KEEP.push_back(cur[0]); KEEPT.push_back(0); }
    // ---- 动作历史（父指针池）----
    vector<int> POOLpar; vector<array<u64,3> > POOLact;
    vector<int> curNode(1,-1);
    vector<array<u64,3> > cd; vector<St> nx; vector<int> nxpar; vector<array<u64,3> > nxact;
    for(int t=0;t<min(TMAX,TBEAM);t++){
        nx.clear(); nxpar.clear(); nxact.clear();
        for(size_t ei=0;ei<cur.size();ei++){ gen(cur[ei],cd);
            for(size_t ci=0;ci<cd.size();ci++){ nx.push_back(step(cur[ei],cd[ci][0],cd[ci][1],cd[ci][2]));
                nxpar.push_back((int)ei); nxact.push_back(cd[ci]); } }
        { vector<int> ord(nx.size()); iota(ord.begin(),ord.end(),0);
          sort(ord.begin(),ord.end(),[&](int a,int b){return memcmp(&nx[a],&nx[b],sizeof(St))<0;});
          vector<St> n2; vector<int> p2; vector<array<u64,3> > a2;
          n2.reserve(nx.size()); p2.reserve(nx.size()); a2.reserve(nx.size());
          for(size_t z=0;z<ord.size();z++){
            if(!n2.empty()&&memcmp(&n2.back(),&nx[ord[z]],sizeof(St))==0) continue;
            n2.push_back(nx[ord[z]]); p2.push_back(nxpar[ord[z]]); a2.push_back(nxact[ord[z]]); }
          nx.swap(n2); nxpar.swap(p2); nxact.swap(a2); }
        int n=nx.size();
        // 一级：线性粗筛
        vector<int> pre;
        {
            vector<ld> sc(n); vector<int> v(n); int kk=min(P0,n);
            auto take=[&](void){ iota(v.begin(),v.end(),0);
                partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return sc[a]>sc[b];});
                for(int i=0;i<kk;i++)pre.push_back(v[i]); };
            int HL[12]={0,3,6,9,12,16,20,25,30,40,50,60};
            ld sws[2]={(ld)0,eS};
            for(int si=0;si<2;si++) for(int hi=0;hi<12;hi++){
                for(int i=0;i<n;i++)sc[i]=sH(nx[i],HL[hi],sws[si]); take(); }
            for(int i=0;i<n;i++)sc[i]=sInf(nx[i]); take();
            for(int i=0;i<n;i++)sc[i]=(ld)nx[i].s; take();
            sort(pre.begin(),pre.end()); pre.erase(unique(pre.begin(),pre.end()),pre.end());
        }
        // 二级：短 rollout 细筛
        {
            int n0=pre.size(); static const int NSC=7; int GJ[NSC]={0,2,4,7,10,14,18};
            vector<array<u64,NSC> > sh(n0);
            for(int i=0;i<n0;i++){
                St c=nx[pre[i]]; int gi=0; u64 kk,pp,qq;
                for(int j=0;j<=18&&gi<NSC;j++){
                    while(gi<NSC&&GJ[gi]==j){
                        St d=c; u64 bst=d.x;
                        for(int u=0;u<14;u++){ act(d,u<6?1:0,kk,pp,qq,t+1+j+u); d=step(d,kk,pp,qq); if(d.x>bst)bst=d.x; }
                        sh[i][gi]=bst; gi++;
                    }
                    if(j<18){ act(c,3,kk,pp,qq,t+1+j); c=step(c,kk,pp,qq); }
                }
            }
            vector<int> keep,v(n0); int kk2=min(P1,n0);
            for(int c2=0;c2<NSC;c2++){
                iota(v.begin(),v.end(),0);
                partial_sort(v.begin(),v.begin()+kk2,v.end(),[&](int a,int b){return sh[a][c2]>sh[b][c2];});
                for(int i=0;i<kk2;i++) keep.push_back(v[i]);
            }
            // 长时域准则也参与细筛前的保留
            sort(keep.begin(),keep.end()); keep.erase(unique(keep.begin(),keep.end()),keep.end());
            vector<int> np; for(size_t i=0;i<keep.size();i++) np.push_back(pre[keep[i]]);
            pre.swap(np);
        }
        int m=pre.size();
        vector<RV> rv(m);
        for(int i=0;i<m;i++) rollout(nx[pre[i]],t+1,rv[i]);
        // 长时域评分
        int mm=min(m,LMAX);
        vector<double> lv(mm);
        for(int i=0;i<mm;i++) lv[i]=longVal(nx[pre[i]],t+1);
        vector<int> sel;
        {
            vector<int> v(m); int kk=min(BS,m);
            for(int pol=0;pol<4;pol++) for(int hi=0;hi<15;hi++){
                int h=HSEL[hi]; if(h>RH)continue;
                iota(v.begin(),v.end(),0);
                partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return rv[a].v[pol][h]>rv[b].v[pol][h];});
                for(int i=0;i<kk;i++) sel.push_back(v[i]);
            }
            for(int bi=0;bi<NBR;bi++){
                iota(v.begin(),v.end(),0);
                partial_sort(v.begin(),v.begin()+kk,v.end(),[&](int a,int b){return rv[a].br[bi]>rv[b].br[bi];});
                for(int i=0;i<kk;i++) sel.push_back(v[i]);
            }
            vector<int> u(mm); iota(u.begin(),u.end(),0);
            int kl=min(BL,mm);
            partial_sort(u.begin(),u.begin()+kl,u.end(),[&](int a,int b){return lv[a]>lv[b];});
            for(int i=0;i<kl;i++) sel.push_back(u[i]);
            sort(sel.begin(),sel.end()); sel.erase(unique(sel.begin(),sel.end()),sel.end());
        }
        // 记录交接种子 + 对束内长时域最优的前 CASHTOP 个状态做 3 段兑现网格
        if(mm>0){
            vector<int> u(mm); iota(u.begin(),u.end(),0);
            int kn=min(max(NKEEP,CASHTOP),mm);
            partial_sort(u.begin(),u.begin()+kn,u.end(),[&](int a,int b){return lv[a]>lv[b];});
            if(t+1<=TKEEP) for(int i=0;i<min(NKEEP,mm);i++){ KEEP.push_back(nx[pre[u[i]]]); KEEPT.push_back(t+1); }
            // 兑现网格的作用对象：多准则并集（longVal / 矿物量 / 各分岔值 / 成长策略远期值）
            vector<int> cand;
            { int nc=min((int)u.size(),CASHTOP);
              for(int z=0;z<nc;z++) cand.push_back(u[z]);
              vector<int> w2(m); int kk2=min(CASHTOP,m);
              auto take2=[&](auto cmp){ iota(w2.begin(),w2.end(),0);
                  partial_sort(w2.begin(),w2.begin()+kk2,w2.end(),cmp);
                  for(int i=0;i<kk2;i++) cand.push_back(w2[i]); };
              take2([&](int a,int b){return nx[pre[a]].x>nx[pre[b]].x;});
              for(int bi2=0;bi2<NBR;bi2++) take2([&,bi2](int a,int b){return rv[a].br[bi2]>rv[b].br[bi2];});
              take2([&](int a,int b){return rv[a].v[3][RH]>rv[b].v[3][RH];});
              take2([&](int a,int b){return rv[a].v[2][RH]>rv[b].v[2][RH];});
              sort(cand.begin(),cand.end()); cand.erase(unique(cand.begin(),cand.end()),cand.end()); }
            int t1=t+1;
            for(size_t z=0;z<cand.size();z++){
                const St&e0=nx[pre[cand[z]]];
                for(int cc=0;cc<=CCMAX;cc+=CSTEP) for(int a=0;a<=CAMAX;a+=CSTEP) for(int bb=0;bb<=CBMAX;bb+=CSTEP){
                    St d=e0;
                    for(int uu=1;uu<=CHOR&&t1+uu<=TMAX;uu++){
                        int pol=(uu<=cc)?4:((uu<=cc+a)?2:((uu<=cc+a+bb)?1:0));
                        u64 k,p,q; act(d,pol,k,p,q,t1+uu-1); d=step(d,k,p,q);
                        if(d.x>F[t1+uu]) F[t1+uu]=d.x;
                    }
                }
            }
        }
        vector<St> ns; vector<int> nsNode;
        for(size_t i=0;i<sel.size();i++){ int gi=pre[sel[i]];
            ns.push_back(nx[gi]);
            POOLpar.push_back(curNode[nxpar[gi]]); POOLact.push_back(nxact[gi]);
            nsNode.push_back((int)POOLpar.size()-1); }
        cur.swap(ns); curNode.swap(nsNode);
        if(TDUMP>0&&(t+1==TDUMP||t+1==TDUMP-10||t+1==TDUMP-20||t+1==TDUMP-30)){ char fn[64]; sprintf(fn,"beam%d.bin",t+1); FILE*fd=fopen(fn,"wb"); int cnt=cur.size(); int tt=t+1;
            fwrite(&cnt,4,1,fd); fwrite(&tt,4,1,fd);
            for(int i=0;i<cnt;i++) fwrite(&cur[i],sizeof(St),1,fd); fclose(fd);
            fprintf(stderr,"dumped %d states at t=%d\n",cnt,tt);
            char fh[64]; sprintf(fh,"hist%d.txt",tt); FILE*fp=fopen(fh,"w");
            fprintf(fp,"%d %d\n",cnt,tt);
            for(int i=0;i<cnt;i++){ vector<array<u64,3> > h; int nd=curNode[i];
                while(nd>=0){ h.push_back(POOLact[nd]); nd=POOLpar[nd]; }
                reverse(h.begin(),h.end());
                fprintf(fp,"%d",(int)h.size());
                for(size_t z=0;z<h.size();z++) fprintf(fp," %llu %llu %llu",h[z][0],h[z][1],h[z][2]);
                fprintf(fp,"\n"); }
            fclose(fp); fprintf(stderr,"histories -> %s\n",fh); }
        if(t%50==0) fprintf(stderr,"t=%d n=%d pre=%d beam=%d F=%llu\n",t,n,m,(int)cur.size(),F[t]);
    }
    INHAND=1;
    // ---------- 交接：ORDER=2 策略跑到底，沿途多切换点兑现 ----------
    {
        vector<u64> FB=F; vector<u64> G(TMAX+1,0);
        for(size_t i=0;i<KEEP.size();i++){
            St c=KEEP[i]; int t0=KEEPT[i];
            for(int t=t0;t<=TMAX;t++){
                if(c.x>G[t])G[t]=c.x;
                if(t==TDUMP) HD.push_back(c);
                for(int cc=0;cc<=CCMAX;cc+=CSTEP) for(int a=0;a<=CAMAX;a+=CSTEP) for(int bb=0;bb<=CBMAX;bb+=CSTEP){
                    St d=c;
                    for(int u=1;u<=CHOR&&t+u<=TMAX;u++){
                        int pol = (u<=cc)?4:((u<=cc+a)?2:((u<=cc+a+bb)?1:0));
                        u64 k,p,q; act(d,pol,k,p,q,t+u); d=step(d,k,p,q);
                        if(d.x>G[t+u])G[t+u]=d.x;
                    }
                }
                if(t==TMAX)break;
                u64 k,p,q; act(c,3,k,p,q,t); c=step(c,k,p,q);
            }
        }
        for(int t=1;t<=TMAX;t++) if(G[t]<G[t-1]) G[t]=G[t-1];
        int nb=0; for(int t=0;t<=TMAX;t++){ if(G[t]>F[t]){nb++; F[t]=G[t];} }
        if(TDUMP>0&&!HD.empty()){ FILE*fd=fopen("hand.bin","wb"); int cnt=HD.size();
            fwrite(&cnt,4,1,fd); fwrite(&TDUMP,4,1,fd);
            for(int i=0;i<cnt;i++) fwrite(&HD[i],sizeof(St),1,fd); fclose(fd);
            fprintf(stderr,"handoff dumped %d states at t=%d\n",cnt,TDUMP); }
        fprintf(stderr,"handoff seeds=%d improved=%d  ",(int)KEEP.size(),nb);
        for(int t=200;t<=TMAX;t+=100) if(FB[t]>0)
            fprintf(stderr,"[t=%d beam=%.5g hand=%.5g] ",t,2.0*FB[t],2.0*G[t]);
        fprintf(stderr,"\n");
    }
    INHAND=0;
    for(int t=1;t<=TMAX;t++) if(F[t]<F[t-1]) F[t]=F[t-1];
    u64 QQ[39]={100,1000,10000,9,79,99,126,166,266,426,666,999,1899,3399,5599,9999,21999,45999,99999,316227,999999,2476413,2745943,9999999,31622776,99999999,316227766,999999999,3162277660ULL,9641625025ULL,10691286350ULL,99999999999ULL,316227766017ULL,999999999999ULL,9999999999999ULL,48610229060556ULL,59770531908338ULL,66277611238091ULL,99999999999999ULL};
    u64 EE[39]={20,200,640,0,10,20,30,40,70,110,160,200,300,410,530,640,800,940,1090,1320,1540,1710,1740,1990,2210,2430,2650,2880,3100,3320,3330,3770,3990,4210,4660,4970,5000,5030,5100};
    int bad=0;
    for(int i=0;i<39;i++){ u64 tg=QQ[i]/2+(QQ[i]&1); int t=0; while(t<=TMAX&&F[t]<tg)t++;
        u64 got=(u64)t*10;
        if(got!=EE[i]){bad++;printf("MISS m=%llu got=%llu exp=%llu diff=%+lld\n",QQ[i],got,EE[i],(long long)got-(long long)EE[i]);}}
    printf("bad=%d F150x2=%llu F510x2=%llu\n",bad,2*F[min(150,TMAX)],2*F[min(510,TMAX)]);
    FILE*fo=fopen("table68.txt","w");
    for(int t=0;t<=TMAX;t++) fprintf(fo,"%lluULL,%s",F[t],(t%8==7)?"\n":" ");
    fclose(fo);
    return 0;
}
