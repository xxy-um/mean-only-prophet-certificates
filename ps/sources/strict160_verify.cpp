#include "interval.hpp"
#include "log_plane.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <map>
#include <list>
#include <vector>
using std::vector;
constexpr int64_t PDEN=10000000000LL,MDEN=4294967296LL,SDEN=4294967296LL,RDEN=1048576;
struct Plane {int owner,branch;int64_t mix,p,q,H,rho;std::string kind;};
struct Proposal {
 int coarse,K;std::array<int64_t,640> params;vector<int64_t> tail;vector<Plane> planes;
 explicit Proposal(const std::string& path){
  std::ifstream f(path);std::string magic;
  if(!(f>>magic>>coarse>>K)||magic!="MEAN_ONLY_AFFINE_V1"||K!=160||coarse%160||coarse<=0)
   throw std::runtime_error("proposal header");
  for(auto& v:params)if(!(f>>v)||v<0)throw std::runtime_error("proposal parameters");
  for(int j=0;j<K;++j)if(params[j]<=0||params[K+j]>params[j])throw std::runtime_error("coefficient feasibility");
  for(int j=0;j<160;++j){if(params[480+j])throw std::runtime_error("eta must be zero");if(j&&params[j]<params[j-1])throw std::runtime_error("nonmonotone b");}
  i128 totalb=0;for(int j=0;j<160;++j)totalb+=params[j];if(totalb>160*i128(PDEN))throw std::runtime_error("threshold too high");
  tail.resize(coarse);for(auto& v:tail)if(!(f>>v)||v<0||v>RDEN)throw std::runtime_error("tail proposal");
  int count;if(!(f>>count)||count<0||count>100000)throw std::runtime_error("plane count");
  vector<i128> totals(3*coarse);
  for(int j=0;j<count;++j){Plane p{};
   if(!(f>>p.owner>>p.mix>>p.kind)||p.owner<0||p.owner>=3*coarse||p.mix<0||p.mix>MDEN)throw std::runtime_error("plane owner or mixture");
   totals[p.owner]+=p.mix;
   if(p.kind=="log") {if(!(f>>p.branch>>p.p>>p.q)||p.owner>=2*coarse||p.branch!=1)throw std::runtime_error("log proposal");}
   else if(p.kind=="C") {if(!(f>>p.rho)||p.rho<0||p.rho>PDEN||p.owner<2*coarse)throw std::runtime_error("C proposal");}
   else throw std::runtime_error("plane kind");
   planes.push_back(p);
  }
  for(i128 x:totals)if(x>MDEN)throw std::runtime_error("mixtures exceed one");
  if(f>>magic)throw std::runtime_error("trailing proposal input");
 }
};
static std::pair<Rat,Rat> route(int j){
 j/=4;
 if(j==16)return {Rat(3,4),Rat(1)};
 if(j==17)return {Rat(7,10),Rat(9,10)};
 if(j>=18&&j<28)return {Rat(1,2),Rat(1,2)};
 if(j==32)return {Rat(1,4),Rat(1,2)};
 if(j>=33&&j<38)return {Rat(1,4),Rat(1,4)};
 return {Rat(),Rat(1)};
}

static I div_constant(I a,Rat r){
 if(r.n>0)return a/r.box();
 if(r.n<0)return -(a/Rat(-r.n,r.d).box());
 throw std::runtime_error("zero constant denominator");
}
static I primitive(Rat A,Rat C,Rat gamma,Rat s){
 if(s.n==0){if(C.n&&A.n>0)return div_constant(A.box()*(log_rat(A)-I::integer(1)),C);return I();}
 Rat z=A+C*s;if(z.n<=0)throw std::runtime_error("kernel argument");
 I out=C.n?div_constant(z.box()*(log_rat(z)-I::integer(1)),C):s.box()*log_rat(A);
 return out-s.box()*(log_rat(gamma*s)-I::integer(1));
}
static I integrate(Rat beta,Rat gamma,Rat v,Rat l,Rat u){
 if(!(l<u)||v.n==v.d)return I();
 Rat one(1),cut=(v-beta)/(one-beta),middle=rmax(l,rmin(cut,u));I ans;
 if(l<middle)ans+=primitive(one-v,gamma,gamma,middle)-primitive(one-v,gamma,gamma,l);
 if(middle<u)ans+=primitive(one-beta,gamma+beta-one,gamma,u)-primitive(one-beta,gamma+beta-one,gamma,middle);
 return ans;
}
struct Common {
 int N;const std::array<int64_t,640> params;vector<I> future;std::map<int,I> negative;
 Common(int n,const std::array<int64_t,640>& p):N(n),params(p),future(n){
  for(int j=0;j<N;++j)future[j]=kernel(Rat(j+1,N),Rat(j,N),true)*I::rat(1,N);
 }
 I kernel(Rat v,Rat a=Rat(),bool restrict=false)const{
  I result;
  for(int c=0;c<160;++c)if(params[160+c]){
   auto bg=route(c);Rat l(c,160),u(c+1,160);
   if(restrict)u=rmin(u,a/bg.second);
   result+=I::rat(params[160+c],PDEN)*integrate(bg.first,bg.second,v,l,u);
  }
  return result;
 }
 I neglog(int t){auto it=negative.find(t);if(it!=negative.end())return it->second;return negative.emplace(t,kernel(Rat(t,N))).first->second;}
 I baseline(int t)const{
  i128 numerator=i128(t)*PDEN;
  for(int c=0;c<160;++c)numerator+=i128(params[c])*std::max((c+1)*(N/160)-std::max(c*(N/160),t),0);
  return I::rat(numerator,i128(N)*PDEN);
 }
};
static Common& common(int N,const Proposal& p){
 static std::map<std::pair<int,std::array<int64_t,640>>,std::unique_ptr<Common>> cache;
 auto key=std::make_pair(N,p.params);auto it=cache.find(key);if(it!=cache.end())return *it->second;
 auto ptr=std::make_unique<Common>(N,p.params);auto& out=*ptr;cache.emplace(key,std::move(ptr));return out;
}

// One fixed affine proof on a complete rectangle. No clock interpolation.
struct Bound {int64_t lower,constLo,constHi,suffix,gradientWidth;int minIndex;};
static Bound evaluate(const Proposal& p,int N,int lo,int width){
 if(N<=0||N%p.coarse||N%160||lo<0||width<=0||lo+width>N||lo/(N/160)!=(lo+width-1)/(N/160))throw std::runtime_error("mesh or rectangle");
 int hi=lo+width,refine=N/p.coarse;Common& co=common(N,p);
 Rat ell(lo,N),upper(hi,N),one(1),half(1,2),h(1,N);I H=h.box();
 vector<I> g=co.future;vector<i128> diffLo(N+1),diffHi(N+1);
 I constant=meet_min(co.baseline(lo),co.baseline(hi))-co.neglog(lo);
 i128 negativeK=0;
 for(int c=0;c<160;++c)negativeK+=i128(p.params[320+c])*std::clamp(hi-c*(N/160),0,N/160);
 constant-=I::rat(negativeK,i128(N)*PDEN);
 for(int j=N-1;j>=0;--j){
  int c=j/(N/160);auto bg=route(c);Rat cc(j+1,N),cut=bg.first+(one-bg.first)*cc;
  g[j]+=-I::rat(p.params[c],PDEN)*H*rmin(cc,upper).box();
  g[j]+=I::rat(p.params[160+c],PDEN)*H*(one-rmax(upper,cut)).box();
  if(j<hi)g[j]+=I::rat(p.params[320+c],i128(N)*PDEN);
 }
 auto addr=[&](int l,int u,I value){if(l<u){diffLo[l]+=value.lo;diffLo[u]-=value.lo;diffHi[l]+=value.hi;diffHi[u]-=value.hi;}};
 for(const Plane& ap:p.planes){
  int typ=ap.owner/p.coarse,cj=ap.owner%p.coarse,l=cj*refine,u=(cj+1)*refine,c=l/(N/160);
  int64_t parameter=typ==0?p.params[c]-p.params[160+c]:typ==1?p.params[160+c]:p.params[320+c];
  if(parameter==0||ap.mix==0)continue;
  I weight=I::rat(parameter,i128(N)*PDEN)*I::rat(ap.mix,MDEN);
  for(int j=l;j<u;++j){
   Rat aa(j,N),cc(j+1,N);
   if(ap.kind=="C"){
    Rat rho(ap.rho,PDEN);I power=pow_unit(rho,cc),A=power*(one+cc-cc*rho).box(),B=power*rho.box();
    constant+=weight*(A-B*H*I::rat(1,8));g[j]-=weight*(A+B*H*half.box());addr(0,j,-weight*B*H);continue;
   }
   auto bg=route(c);Rat d,du;
   if(typ==0){d=rmax(ell-cc,Rat());du=rmax(upper-aa,Rat());}
   else {d=rmax(rmin(ell,bg.first+(one-bg.first)*aa)-cc,Rat());du=rmax(rmin(upper,bg.first+(one-bg.first)*cc)-aa,Rat());}
   if(d.n==0)continue;
   Rat delta=rmin(aa,(one-du)*half),where=(cc-delta)*Rat(N);
   int start=std::clamp(int(ceildiv(where.n,where.d)),0,j);
   Rat w(j-start,N),slopeP(ap.p,SDEN),slopeQ(ap.q,SDEN);LogPlane contribution;
   try{auto plane=verify_old_log_plane(d,Rat(1),w,slopeP,slopeQ);contribution={weight*plane.p,weight*plane.q,weight*plane.h};}
   catch(const std::runtime_error& e){if(std::string(e.what())!="fixed-point range guard")throw;contribution=verify_old_weighted_log_plane(d,Rat(1),w,slopeP,slopeQ,weight);}
   constant+=contribution.h;g[j]+=contribution.p;addr(start,j,contribution.q*H);
  }
 }
 i128 activeLo=0,activeHi=0;
 for(int j=0;j<N;++j){activeLo+=diffLo[j];activeHi+=diffHi[j];g[j]+=I(checked(activeLo),checked(activeHi));}
 i128 suffix=0,minimum=0;int index=N;int64_t gradientWidth=0;
 for(int j=N-1;j>=0;--j){suffix+=g[j].lo;gradientWidth=std::max(gradientWidth,g[j].hi-g[j].lo);if(suffix<minimum){minimum=suffix;index=j;}}
 return {checked(i128(constant.lo)+minimum),constant.lo,constant.hi,checked(minimum),gradientWidth,index};
}
int main(int argc,char** argv){
 try{
  if(argc!=5)throw std::runtime_error("usage: clean_verify TASKS TARGET_NUM TARGET_DEN");
  std::ifstream tasks(argv[1]);int64_t tn=std::stoll(argv[2]),td=std::stoll(argv[3]);
  // argv[4] is the fixed format version, preventing accidental old CLI use.
  if(!tasks||td<=0||std::string(argv[4])!="RECTANGLE_V1")throw std::runtime_error("tasks, target or version");
  std::map<std::string,std::unique_ptr<Proposal>> sources;std::string source;int N,lo,width;bool all=true;
  while(tasks>>source>>N>>lo>>width){
   auto it=sources.find(source);if(it==sources.end())it=sources.emplace(source,std::make_unique<Proposal>(source)).first;
   Bound r=evaluate(*it->second,N,lo,width);bool pass=i128(r.lower)*td>=i128(tn)*Q;all=all&&pass;
   std::cout<<"{\"status\":\""<<(pass?"PASS":"FAIL")<<"\",\"source\":\""<<source<<"\",\"N\":"<<N<<",\"lo\":"<<lo<<",\"width\":"<<width<<",\"lower_numerator\":"<<r.lower<<",\"denominator\":"<<Q<<",\"bound\":"<<std::setprecision(17)<<double(r.lower)/Q<<",\"constant_lo\":"<<r.constLo<<",\"constant_hi\":"<<r.constHi<<",\"min_suffix\":"<<r.suffix<<",\"min_suffix_index\":"<<r.minIndex<<",\"max_gradient_width\":"<<r.gradientWidth<<"}\n"<<std::flush;
  }
  if(!tasks.eof())throw std::runtime_error("malformed tasks");return all?0:1;
 }catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<"\n";return 2;}
}
