#pragma once
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
using i128=__int128_t;
using u128=__uint128_t;
constexpr int SHIFT=56;
constexpr int64_t Q=int64_t(1)<<SHIFT;
constexpr int64_t GUARD=int64_t(1)<<62;
inline i128 floordiv(i128 a,i128 b){
 if(b<=0)throw std::runtime_error("nonpositive divisor");
 i128 t=a/b,r=a%b;return t-(r<0);
}
inline i128 ceildiv(i128 a,i128 b){
 if(b<=0)throw std::runtime_error("nonpositive divisor");
 i128 t=a/b,r=a%b;return t+(r>0);
}
inline int64_t checked(i128 x){
 if(x<=-GUARD||x>=GUARD)throw std::runtime_error("fixed-point range guard");
 return int64_t(x);
}
struct I{
 int64_t lo,hi;
 I():lo(0),hi(0){}
 I(int64_t l,int64_t h):lo(l),hi(h){
  if(l>h||l<=-GUARD||h>=GUARD)throw std::runtime_error("invalid interval");
 }
 static I integer(int v){return I(checked(i128(v)*Q),checked(i128(v)*Q));}
 static I rat(i128 n,i128 d){
  i128 scaled;
  if(__builtin_mul_overflow(n,i128(Q),&scaled))throw std::runtime_error("rational scaling overflow");
  return I(checked(floordiv(scaled,d)),checked(ceildiv(scaled,d)));
 }
 double mid()const{return (double(lo)+double(hi))/(2*Q);}
};
inline I operator+(I a,I b){return I(checked(i128(a.lo)+b.lo),checked(i128(a.hi)+b.hi));}
inline I operator-(I a,I b){return I(checked(i128(a.lo)-b.hi),checked(i128(a.hi)-b.lo));}
inline I operator-(I a){return I(-a.hi,-a.lo);}
inline I operator*(I a,I b){
 i128 p[4]={i128(a.lo)*b.lo,i128(a.lo)*b.hi,i128(a.hi)*b.lo,i128(a.hi)*b.hi};
 i128 l=*std::min_element(p,p+4),h=*std::max_element(p,p+4);
 return I(checked(floordiv(l,Q)),checked(ceildiv(h,Q)));
}
inline I operator/(I a,I b){
 if(b.lo<=0)throw std::runtime_error("division interval touches zero");
 i128 l=std::min({floordiv(i128(a.lo)*Q,b.lo),floordiv(i128(a.lo)*Q,b.hi),floordiv(i128(a.hi)*Q,b.lo),floordiv(i128(a.hi)*Q,b.hi)});
 i128 h=std::max({ceildiv(i128(a.lo)*Q,b.lo),ceildiv(i128(a.lo)*Q,b.hi),ceildiv(i128(a.hi)*Q,b.lo),ceildiv(i128(a.hi)*Q,b.hi)});
 return I(checked(l),checked(h));
}
inline I& operator+=(I& a,I b){a=a+b;return a;}
inline I& operator-=(I& a,I b){a=a-b;return a;}
inline I times(I a,int b){return a*I::integer(b);}
inline I meet_min(I a,I b){return I(std::min(a.lo,b.lo),std::min(a.hi,b.hi));}
inline std::string str128(i128 v){
 if(v==0)return "0";
 bool neg=v<0;u128 x=neg?u128(-v):u128(v);std::string s;
 while(x){s.push_back(char('0'+x%10));x/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;
}
struct Rat{
 int64_t n,d;
 static int64_t gcd(int64_t a,int64_t b){a=a<0?-a:a;while(b){int64_t r=a%b;a=b;b=r;}return a;}
 Rat(int64_t a=0,int64_t b=1):n(a),d(b){
  if(d<0){n=-n;d=-d;}if(d==0)throw std::runtime_error("zero rational denominator");
  int64_t g=gcd(n,d);n/=g;d/=g;
 }
 double dbl()const{return double(n)/d;}
 I box()const{return I::rat(n,d);}
};
inline Rat operator+(Rat a,Rat b){
 int64_t g=Rat::gcd(a.d,b.d);
 i128 n=i128(a.n)*(b.d/g)+i128(b.n)*(a.d/g),d=i128(a.d/g)*b.d;
 i128 an=n<0?-n:n,ad=d;
 while(ad){i128 r=an%ad;an=ad;ad=r;}
 if(an){n/=an;d/=an;}
 if(n<=INT64_MIN||n>=INT64_MAX||d>=INT64_MAX)throw std::runtime_error("rational overflow");
 return Rat(int64_t(n),int64_t(d));
}
inline Rat operator-(Rat a,Rat b){return a+Rat(-b.n,b.d);}
inline Rat operator*(Rat a,Rat b){
 int64_t g1=Rat::gcd(a.n,b.d),g2=Rat::gcd(b.n,a.d);
 i128 n=i128(a.n/g1)*(b.n/g2),d=i128(a.d/g2)*(b.d/g1);
 if(n<=INT64_MIN||n>=INT64_MAX||d>=INT64_MAX)throw std::runtime_error("rational overflow");
 return Rat(int64_t(n),int64_t(d));
}
inline Rat operator/(Rat a,Rat b){return a*Rat(b.d,b.n);}
inline bool operator<(Rat a,Rat b){return i128(a.n)*b.d<i128(b.n)*a.d;}
inline bool operator<=(Rat a,Rat b){return !(b<a);}
inline Rat rmin(Rat a,Rat b){return a<b?a:b;}
inline Rat rmax(Rat a,Rat b){return a<b?b:a;}
inline I positive_log_series(int64_t n,int64_t d){
 // This is log(n/d), with 1 <= n/d <= 2.
 if(n<d||i128(n)>2*i128(d))throw std::runtime_error("log range reduction");
 I z=I::rat(n-d,n+d),z2=z*z,term=z,sum;
 for(int q=0;q<24;++q){sum+=term/I::integer(2*q+1);term=term*z2;}
 // For z <= 1/3, the omitted atanh tail is < 2^(-79), hence < 1/Q.
 I out=times(sum,2);out.hi=checked(i128(out.hi)+1);return out;
}
inline I log_dyadic(int64_t x){
 if(x<=0)throw std::runtime_error("log nonpositive argument");
 static const I log2=positive_log_series(2,1);
 static std::unordered_map<int64_t,I> cache;
 auto it=cache.find(x);if(it!=cache.end())return it->second;
 int64_t n=x,d=Q;int power=0;
 while(n<d){n*=2;--power;}
 while(i128(n)>=2*i128(d)){d*=2;++power;}
 I out=positive_log_series(n,d)+times(log2,power);
 if(cache.size()>100000)cache.clear();
 cache.emplace(x,out);return out;
}
inline I log_box(I x){
 if(x.lo<=0)throw std::runtime_error("log interval touches zero");
 I a=log_dyadic(x.lo),b=log_dyadic(x.hi);return I(a.lo,b.hi);
}
inline I log_rat(Rat r){return log_box(r.box());}
// Enclose exp(x) for x <= 0. First evaluate exp(-x/2^k), with
// 0 <= -x/2^k <= 1/2. The degree-20 Taylor remainder is less
// than e^(1/2) (1/2)^21 / 21! < 2^-84 < 1/Q.
inline I exp_nonpositive_point(int64_t x){
 if(x>0)throw std::runtime_error("positive exp argument");
 if(x==0)return I::integer(1);
 int k=0;int64_t divisor=1;
 while(i128(-x)>i128(Q/2)*divisor){divisor*=2;++k;}
 I z=I::rat(-i128(x),i128(Q)*divisor),term=I::integer(1),sum=term;
 for(int n=1;n<=20;++n){term=term*z/I::integer(n);sum+=term;}
 sum.hi=checked(i128(sum.hi)+1);
 I out=I::integer(1)/sum;
 for(int n=0;n<k;++n)out=out*out;
 return out;
}
inline I exp_nonpositive(I x){
 if(x.hi>0)throw std::runtime_error("positive exp interval");
 I lo=exp_nonpositive_point(x.lo),hi=exp_nonpositive_point(x.hi);
 return I(lo.lo,hi.hi);
}
inline I pow_unit(Rat base,Rat exponent){
 if(base.n<0||base.d<base.n||exponent.n<0||exponent.d<exponent.n)
  throw std::runtime_error("unit power domain");
 if(exponent.n==0||base.n==base.d)return I::integer(1);
 if(base.n==0)return I();
 return exp_nonpositive(log_rat(base)*exponent.box());
}
