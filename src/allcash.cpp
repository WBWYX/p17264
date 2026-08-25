// 对文件中全部状态施加步长1的三段兑现网格（不含成长前缀），求 max x(T)
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
static inline bool ok(const St&e,u64 k,u64 p,u64 q){
    return p+q<=e.w && k<=e.b && k<=e.s && 25*k+50*p+200*q<=e.x;
}
int GNC=98;
static void act(const St&e,int pol,u64&k,u64&p,u64&q){
    k=p=q=0; if(pol==0) return;
    u64 x1=e.x, need=0;
    if(pol==2){ u64 Dm=((u64)e.b*GNC+19)/20; need=(Dm>e.s)?(Dm-e.s+7)/8:0; }
    else if(pol>=3){ long long v=2LL*(long long)e.ab[1]+(long long)e.ab[2]+(long long)e.ab[3]
                              +2LL*(long long)e.b-(long long)e.s
                              -(long long)e.as_[1]-(long long)e.as_[2]-(long long)e.as_[3];
                     need=(v>0)?(u64)((v+7)/8):0; }
    if(pol>=2){ p=min(min(need,e.w),x1/50); x1-=50*p; }
    k=min(min(e.b,e.s),x1/25); x1-=25*k;
    if(pol==3) q=min(e.w-p,x1/200);
}
int main(int argc,char**argv){
    const char*fn=argc>1?argv[1]:"beam135.bin";
    int T=argc>2?atoi(argv[2]):171;
    int GMAX=argc>3?atoi(argv[3]):0;    // 允许的成长前缀步数上限（0=不成长）
    int CCM=argc>4?atoi(argv[4]):20, CAM=argc>5?atoi(argv[5]):28, CBM=argc>6?atoi(argv[6]):14;
    if(argc>7)GNC=atoi(argv[7]);
    FILE*f=fopen(fn,"rb"); if(!f){printf("no file\n");return 1;}
    int cnt,t0; if(fread(&cnt,4,1,f)!=1||fread(&t0,4,1,f)!=1)return 1;
    vector<St> all(cnt);
    for(int i=0;i<cnt;i++) if(fread(&all[i],sizeof(St),1,f)!=1){cnt=i;all.resize(cnt);break;}
    fclose(f);
    int Hz=T-t0;
    printf("%s: %d 状态 t0=%d T=%d 窗口=%d  成长前缀<=%d  网格步长1 cc<=%d a<=%d b<=%d GNC=%d\n",
           fn,cnt,t0,T,Hz,GMAX,CCM,CAM,CBM,GNC);
    fflush(stdout);
    u64 BEST=0; int bi=-1,bg=0,bcc=0,ba=0,bb=0;
    for(int i=0;i<cnt;i++){
        St g0=all[i];
        for(int g=0;g<=GMAX&&g<=Hz;g++){
            if(g>0){ u64 k,p,q; act(g0,3,k,p,q); if(!ok(g0,k,p,q)) break; g0=step(g0,k,p,q); }
            int hz=Hz-g;
            for(int cc=0;cc<=CCM&&cc<=hz;cc++)
            for(int a=0;a+cc<=hz&&a<=CAM;a++)
            for(int b2=0;b2+a+cc<=hz&&b2<=CBM;b2++){
                St d=g0; bool good=true;
                for(int u=0;u<hz;u++){
                    int pol=(u<cc)?4:((u<cc+a)?2:((u<cc+a+b2)?1:0));
                    u64 k,p,q; act(d,pol,k,p,q); if(!ok(d,k,p,q)){good=false;break;}
                    d=step(d,k,p,q);
                }
                if(good&&d.x>BEST){BEST=d.x;bi=i;bg=g;bcc=cc;ba=a;bb=b2;}
            }
        }
        if((i+1)%500==0){ fprintf(stderr,"%d/%d  best2=%llu\n",i+1,cnt,2*BEST); fflush(stderr); }
    }
    printf("最好 x(%d)*2 = %llu  (状态#%d, 成长%d + cc=%d a=%d b=%d + 采矿%d)\n",
           T,2*BEST,bi,bg,bcc,ba,bb,Hz-bg-bcc-ba-bb);
    return 0;
}
