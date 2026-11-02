#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>
#include <cassert>
#include <cfloat>
#include <cfenv>
#if defined(__FAST_MATH__) || defined(NDEBUG)
#error Strict arithmetic and enabled assertions are required.
#endif
static_assert(FLT_RADIX==2 && DBL_MANT_DIG>=53, "binary64 or better required");
#ifndef COEFFICIENT_FILE
#define COEFFICIENT_FILE "coefficients_N640_Q65547_validated.h"
#endif
#include COEFFICIENT_FILE
using namespace std;
using i64=int64_t;using i128=__int128_t;
constexpr i64 U=i64(1)<<40,L=i64(1)<<30;
constexpr double eps=0x1p-32;
int m, bound_mask=15; i64 Cscale; vector<i64> table,nexttable; vector<int> hints,nexthints; vector<i64> grid; vector<int> lookup;
inline int indexof(int i,int j){return i*(i+1)/2+j;}
inline i64 floordiv(i128 a,i64 b){return a>=0?i64(a/b):-i64((-a+b-1)/b);}
inline i64 interpolate(double y,double x){
 i64 Y=min(Cscale,max(i64(0),i64(ceil((y+eps)*Cscale))));
 i64 X=max(i64(0),i64(floor((x-eps)*Cscale)));
 if(X>=Y)return 0;
 int i=lookup[min(65536,int(Y>>14))];while(i>0 && grid[i]>=Y)--i;while(i<m-1 && grid[i+1]<Y)++i;
 int j=lookup[min(65536,int(X>>14))];while(j<m-1 && grid[j+1]<=X)++j;j=min(i,j);
 i64 dy=grid[i+1]-grid[i],dx=grid[j+1]-grid[j],a=Y-grid[i],b=X-grid[j];
 i64 denominator=dx*dy; i128 sum;
 if(b*dy<=a*dx)sum=i128((dy-a)*dx)*table[indexof(i,j)]+i128(a*dx-b*dy)*table[indexof(i+1,j)]+i128(b*dy)*table[indexof(i+1,j+1)];
 else sum=i128((dx-b)*dy)*table[indexof(i,j)]+i128(b*dy-a*dx)*table[indexof(i,j+1)]+i128(a*dx)*table[indexof(i+1,j+1)];
 return floordiv(sum,denominator);
}

// Keep every original rate-tree interval and split each finite original
// leaf once. The final unbounded interval is not split.
inline int partition_mid(int lo,int hi){
 assert(hi-lo>1);
 if(hi-lo==2)return lo+1;
 assert(lo%2==0);
 assert(hi%2==0 || hi==Q);
 int oldlo=lo/2,oldhi=hi==Q?(Q+1)/2:hi/2;
 int mid=2*((oldlo+oldhi)/2);
 assert(lo<mid && mid<hi);
 return mid;
}


constexpr int potential_bits=16;
constexpr double potential_units=0x1p24;
int slack_bits=28;
i64 slack_units=i64(1)<<(40-28);
struct Node{
 double x,y,t,s,bmin,p,z,negative,positive,global_floor;
 mutable unsigned calls=0,leaves=0,coverage=0;
 inline double ph(int r)const{return 1-y+y*(t+ii[r]+ee[r]*(1-t-dt));}
 inline i64 A(int lo){
  double phi=ph(lo), endpoint_numerator=cap*x+y*s*dt*vv[lo];
  i64 continuation=0;
  if(endpoint_numerator<cap*y+eps*y){double xp=ee[lo]*endpoint_numerator/(cap*phi),yp=y*ee[lo]/phi;continuation=interpolate(yp,xp);}
  i64 multiplier=i64(ceil((phi+eps)*U))+1;
  return floordiv(i128(multiplier)*continuation,U);
 }
 inline double log_lower(double ratio) const {
  assert(ratio>=1);
  int exponent;double mant=frexp(ratio,&exponent)*2;--exponent;
  double w=(mant-1)/(mant+1),w2=w*w,term=w,sum=w;
  for(int d=3;d<=23;d+=2){term*=w2;sum+=term/d;}
  constexpr double ln2_lower=0x1.62e42fefa39eep-1;
  return 2*sum+exponent*ln2_lower-eps;
 }
 inline i64 B(int hi,bool shifted=true){
  double running=min(p,0.0)*dt-eps;
  if(hi<Q){
   double alpha=(y*bmin*ee[hi]-min(p,0.0)*ii[hi]/dt)*dt*dt;
   double beta=bmin*(1-y+t*y)*dt+y*bmin*ii[hi]*dt;
   double ze=beta*er[hi], ell=ze-alpha/2;
   double base=min(p,0.0)*dt+max(p,0.0)*ii[hi];
   if(alpha>eps && beta>eps){
    double freebound;
    if(beta<=alpha){double ratio=beta/alpha;freebound=beta*(1-ratio/4);}
    else freebound=alpha*(.75+.5*log_lower(beta/alpha));
    running=max(running,base+freebound-eps);
   }
   if(shifted && rr[hi]<=64 && ze>alpha+eps){
    double sharp=base+ze+alpha*rr[hi]/2-alpha*alpha/(24*ell);
    running=max(running,sharp-eps);
   }
  }
  return i64(floor((running-eps)*U))-1;
 }
 inline i64 rectA(int lo,int hi){
  double ratio=x/y+s*dt*vv[lo]/cap;
  if(ratio>=1+eps)return 0;
  ratio=min(1.,ratio);
  double phi=ph(lo), ya=y*ee[lo]/phi, xlo=0;
  if(hi<Q && rr[hi]*dt<=4-eps){double yb=y*ee[hi]/ph(hi);xlo=yb*ratio;}
  i64 v=interpolate(ya,xlo);
  i64 mult=i64(ceil((phi+eps)*U))+1;
  return floordiv(i128(mult)*v,U);
 }
 inline i64 Bnode(int lo,int hi){
  i64 free=B(hi,false);
  if(hi==Q || rr[hi]>64)return free;
  double negp=min(p,0.0),base=negp*dt+max(p,0.0)*ii[hi];
  double alpha_a=(y*bmin*ee[lo]-negp*ii[lo]/dt)*dt*dt;
  double alpha_b=(y*bmin*ee[hi]-negp*ii[hi]/dt)*dt*dt;
  double beta_b=bmin*(1-y+t*y)*dt+y*bmin*ii[hi]*dt;
  double z_b=beta_b*er[hi];
  constexpr double eta=0x1p-48;
  double aup=alpha_a+16*eta,alow=max(0.,alpha_b-16*eta),zlow=max(0.,z_b-16*eta);
  if(zlow>aup+eps){
   double rlow=max(0.,rr[lo]-16*eta),ell=zlow-aup/2;
   double sharp=base+zlow+alow*rlow/2-aup*aup/(24*ell);
   free=max(free,i64(floor((sharp-2*eps)*U))-1);
  }
  return free;
 }
 inline i64 lower(int lo,int hi){return rectA(lo,hi)+Bnode(lo,hi);}
 inline i64 allowance(int index){
  // Candidate initialization only: these nonnegative integer subtractions do
  // not enter the fixed original FRLP rows or any node pruning bound.
  double w=y-x;
  double future=ee[index]*max(0.,w-y*s*dt*vv[index]/cap);
  double drop=max(0.,min(w,w-future));
  return i64(ceil(drop*potential_units));
 }
 inline i64 leaf(int index){leaves++;return A(index)+B(index+1)-allowance(index);}
 void search(int lo,int hi,i64 &best,int &which){
  calls++;
  if(lower(lo,hi)>=best){coverage+=hi-lo;return;}
  if(hi-lo==1){coverage++;i64 v=leaf(lo)-slack_units;if(v<best){best=v;which=lo;}return;}
  int mid=partition_mid(lo,hi);
  if(which<mid){search(lo,mid,best,which);search(mid,hi,best,which);}
  else{search(mid,hi,best,which);search(lo,mid,best,which);}
 }
 i64 solve(int &which){
  i64 best=min(i64(0),leaf(which)-slack_units);
  if(which!=Q-1){i64 v=leaf(Q-1)-slack_units;if(v<best){best=v;which=Q-1;}}
  if(lower(0,Q)>=best){coverage=Q;assert(coverage==Q);return best;}
  int ll=0,hh=Q-1;
  for(int it=0;it<18 && hh-ll>2;it++){
   int aa=ll+(hh-ll)/3,bb=hh-(hh-ll)/3;
   i64 va=leaf(aa)-slack_units,vb=leaf(bb)-slack_units;
   if(va<best){best=va;which=aa;}if(vb<best){best=vb;which=bb;}
   if(va<=vb)hh=bb-1;else ll=aa+1;
  }
  for(int a:{ll,(ll+hh)/2,hh}){i64 v=leaf(a)-slack_units;if(v<best){best=v;which=a;}}
  search(0,Q,best,which);assert(coverage==Q);return best;
 }
};
void lower_monotone(){
 for(int i=1;i<=m;i++)for(int j=i-1;j>=0;j--)table[indexof(i,j)]=min(table[indexof(i,j)],table[indexof(i,j+1)]);
 for(int i=1;i<=m;i++)for(int j=0;j<i;j++)table[indexof(i,j)]=min(table[indexof(i,j)],table[indexof(i-1,j)]);
}
void check_monotone(){
 for(int i=1;i<=m;i++)for(int j=0;j<i;j++){
  assert(table[indexof(i,j)]<=table[indexof(i,j+1)]);
  assert(table[indexof(i,j)]<=table[indexof(i-1,j)]);
 }
}
void save_checkpoint(const string& dir,int k){
 string prefix=dir+"/layer_"+to_string(k);
 {ofstream f(prefix+".table.tmp",ios::binary);f.write((const char*)table.data(),table.size()*sizeof(i64));assert(f);}
 {ofstream f(prefix+".next.tmp",ios::binary);f.write((const char*)nexttable.data(),nexttable.size()*sizeof(i64));assert(f);}
 {ofstream f(prefix+".hints.tmp",ios::binary);f.write((const char*)hints.data(),hints.size()*sizeof(int));assert(f);}
 for(string suffix:{".table",".next",".hints"})assert(rename((prefix+suffix+".tmp").c_str(),(prefix+suffix).c_str())==0);
 {ofstream f(dir+"/latest.k.tmp");f<<k<<"\n";assert(f);}
 assert(rename((dir+"/latest.k.tmp").c_str(),(dir+"/latest.k").c_str())==0);
}

int main(int argc,char**argv){
 assert(std::fesetround(FE_TONEAREST)==0);
 assert(argc>=4);int threads=atoi(argv[1]);string knotsfile=argv[2],dir=argv[3];
 if(argc>4)slack_bits=atoi(argv[4]);assert(slack_bits>=20&&slack_bits<=40);slack_units=i64(1)<<(40-slack_bits);
 bool resume=argc>5 && string(argv[5])=="resume";
 int stop_k=argc>6?atoi(argv[6]):0;
 ifstream knots(knotsfile);assert(knots);knots>>m;grid.resize(m+1);for(auto &v:grid)knots>>v;assert(knots);
 omp_set_num_threads(threads);Cscale=L;assert(grid[0]==0&&grid[m]==L);for(int i=1;i<=m;i++)assert(grid[i]>grid[i-1]);lookup.resize(65537);int gi=0;for(int li=0;li<=65536;li++){while(gi<m && grid[gi+1]<=(i64(li)<<14))++gi;lookup[li]=gi;}
 int size=(m+1)*(m+2)/2;table.assign(size,0);nexttable=table;hints.assign(size,0);nexthints=hints;
 int start_k=N-1;
 if(resume){int saved;ifstream meta(dir+"/latest.k");meta>>saved;assert(meta);start_k=saved-1;string prefix=dir+"/layer_"+to_string(saved);ifstream tf(prefix+".table",ios::binary);tf.read((char*)table.data(),size*sizeof(i64));assert(tf);ifstream hf(prefix+".hints",ios::binary);hf.read((char*)hints.data(),size*sizeof(int));assert(hf);check_monotone();}
 double start=omp_get_wtime();cerr<<setprecision(15)<<"literal leaf solver m="<<m<<" N="<<N<<" Q="<<Q<<" slack_bits="<<slack_bits<<" potential_bits="<<potential_bits<<" start_k="<<start_k<<"\n";
 for(int k=start_k;k>=stop_k;k--){double t=k*dt,s=bs[k],bmin=bl[k],D=DD[k];unsigned long long calls=0,leaves=0;
  #pragma omp parallel for schedule(dynamic,2) reduction(+:calls,leaves)
  for(int i=1;i<=m;i++){assert(std::fesetround(FE_TONEAREST)==0);assert(std::fegetround()==FE_TONEAREST);for(int j=0;j<i;j++){
   double y=double(grid[i])/L,x=double(grid[j])/L,p=cap*x+y*(D-bmin);
   Node node={x,y,t,s,bmin,p,0,0,0,-INFINITY};int which=hints[indexof(i,j)];
   nexttable[indexof(i,j)]=node.solve(which);nexthints[indexof(i,j)]=which;calls+=node.calls;leaves+=node.leaves;
  }}
  table.swap(nexttable);hints.swap(nexthints);
  vector<i64> raw=table;lower_monotone();check_monotone();
  i64 lowered=0,maxlower=0;for(int a=0;a<size;a++)if(table[a]<raw[a]){lowered++;maxlower=max(maxlower,raw[a]-table[a]);}
  cerr<<k<<" "<<double(HU+table[indexof(m,0)])/U<<" calls="<<calls<<" leaves="<<leaves<<" lowered="<<lowered<<" maxlower="<<maxlower<<" seconds="<<omp_get_wtime()-start<<"\n";
  if(k%8==0||k==635||k==start_k||k==stop_k)save_checkpoint(dir,k);
 }
 i64 numerator=HU+table[indexof(m,0)];
 cout<<"{\"m\":"<<m<<",\"N\":"<<N<<",\"Q\":"<<Q<<",\"slack_bits\":"<<slack_bits<<",\"potential_bits\":"<<potential_bits<<",\"last_k\":"<<stop_k<<",\"numerator\":"<<numerator<<",\"denominator\":"<<U<<",\"lower_bound\":"<<setprecision(17)<<double(numerator)/U<<",\"seconds\":"<<omp_get_wtime()-start<<"}\n";
}
