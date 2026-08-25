// 开局精确全枚举：每步枚举全部可行 (k,p,q)，仅去重（可选支配剪枝），统计状态数
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
// A 支配 B：所有分量都 >=
static inline bool dom(const St&a,const St&b){
    if(a.x<b.x||a.s<b.s||a.w<b.w||a.b<b.b) return false;
    for(int i=1;i<H;i++){
        if(a.aw[i]<b.aw[i]||a.am[i]<b.am[i]||a.as_[i]<b.as_[i]||a.ab[i]<b.ab[i]) return false;
    }
    return true;
}
int main(int argc,char**argv){
    int TMAX=argc>1?atoi(argv[1]):40;
    size_t LIM=argc>2?(size_t)atoll(argv[2]):3000000;
    int DOMON=argc>3?atoi(argv[3]):1; int DUMPT=argc>4?atoi(argv[4]):0;
    vector<St> cur(1); memset(&cur[0],0,sizeof(St));
    cur[0].x=25;cur[0].s=6;cur[0].w=4;cur[0].b=1;
    vector<St> nx;
    for(int t=0;t<TMAX;t++){
        nx.clear();
        for(size_t i=0;i<cur.size();i++){
            const St&e=cur[i];
            u64 kM=min(min(e.b,e.s),e.x/25);
            for(u64 k=0;k<=kM;k++){
                u64 x1=e.x-25*k;
                u64 pM=min(e.w,x1/50);
                for(u64 p=0;p<=pM;p++){
                    u64 x2=x1-50*p, qM=min(e.w-p,x2/200);
                    for(u64 q=0;q<=qM;q++) nx.push_back(step(e,k,p,q));
                }
            }
        }
        sort(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))<0;});
        nx.erase(unique(nx.begin(),nx.end(),[](const St&a,const St&b){return memcmp(&a,&b,sizeof(St))==0;}),nx.end());
        size_t raw=nx.size();
        if(DOMON&&raw<=40000){
            // O(n^2) 支配剪枝（只在规模允许时做）
            vector<char> dead(raw,0);
            for(size_t i=0;i<raw;i++){ if(dead[i])continue;
                for(size_t j=0;j<raw;j++){ if(i==j||dead[j])continue;
                    if(dom(nx[i],nx[j])) dead[j]=1; } }
            vector<St> keep; for(size_t i=0;i<raw;i++) if(!dead[i]) keep.push_back(nx[i]);
            nx.swap(keep);
        }
        u64 bx=0; for(size_t i=0;i<nx.size();i++) if(nx[i].x>bx)bx=nx[i].x;
        printf("t=%2d  去重后=%zu  支配剪枝后=%zu  maxX*2=%llu\n",t+1,raw,nx.size(),2*bx);
        fflush(stdout);
        if(DUMPT>0&&t+1==DUMPT){ FILE*fo=fopen("exact.bin","wb"); int c2=nx.size(),tt=t+1;
            fwrite(&c2,4,1,fo); fwrite(&tt,4,1,fo);
            for(int i=0;i<c2;i++) fwrite(&nx[i],sizeof(St),1,fo); fclose(fo);
            printf("已导出 %d 个精确前沿状态到 exact.bin (t=%d)\n",c2,tt); }
        if(nx.size()>LIM){ printf("超过上限 %zu，停止\n",LIM); break; }
        cur.swap(nx);
    }
    return 0;
}
