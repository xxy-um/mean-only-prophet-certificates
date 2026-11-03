#pragma once
#include "interval.hpp"
#include <cmath>
#include <array>
struct LogPlane { I p,q,h; };
inline LogPlane verify_old_log_plane(Rat d,Rat a,Rat w,Rat p,Rat q){
 if(d.n==0)return {I(),I(),I()};
 if(d.n<0||d.d<=d.n||a.n<=0||w.n<0||w.d<=w.n||
    !((d+a*w)<Rat(1)))throw std::runtime_error("log-plane triangle");
 if(!(a.n==a.d || i128(a.n)*d.d==i128(a.d)*(d.d-d.n)))
  throw std::runtime_error("log-plane unsupported a");
 I D=d.box(),A=a.box(),W=w.box(),P=p.box(),Qq=q.box(),one=I::integer(1);
 auto proposal_derivative=[&](double u,int edge){
  double dd=d.dbl(),aa=a.dbl(),ww=w.dbl(),pp=p.dbl(),qq=q.dbl();
  if(edge==0)return dd/(1-dd*u)-dd-pp;
  if(edge==1)return -aa*ww/(1-ww*u)+(dd+aa*ww)/(1-(dd+aa*ww)*u)-dd-pp-qq*ww;
  return -aa*ww/(1-ww*u)+aa*ww/(1-dd-aa*ww*u)-qq*ww;
 };
 int64_t lower=GUARD-1;
 for(int edge=0;edge<3;++edge){
  double lo=0,hi=1,u;
  if(proposal_derivative(0,edge)>=0)u=0;
  else if(proposal_derivative(1,edge)<=0)u=1;
  else {for(int k=0;k<60;++k){u=(lo+hi)/2;if(proposal_derivative(u,edge)<0)lo=u;else hi=u;}u=(lo+hi)/2;}
  constexpr int64_t UQ=int64_t(1)<<30;
  Rat ur(std::clamp<int64_t>(std::llround(u*UQ),0,UQ),UQ);I U=ur.box();
  I X=edge==2?one:U,Y=edge==0?I():W*U;
  I Z=one-D*X-A*Y;
  I phi=A*log_box(one-Y)-log_box(Z)-D*X-P*X-Qq*Y;
  I derivative;
  if(edge==0)derivative=D/Z-D-P;
  else if(edge==1)derivative=-A*W/(one-Y)+(D+A*W)/Z-D-P-Qq*W;
  else derivative=-A*W/(one-Y)+A*W/Z-Qq*W;
  I edgeLower=phi-U*derivative+meet_min(I(),derivative);
  lower=std::min(lower,edgeLower.lo);
 }
 return {P,Qq,I(lower,lower)};
}

// Scale an interval by an exact rational without boxing the potentially
// large rational itself. Endpoints stay at the original 2^56 scale.
inline I interval_times_rat(I value,Rat ratio){
 i128 x=i128(value.lo)*ratio.n,y=i128(value.hi)*ratio.n;
 return I(checked(std::min(floordiv(x,ratio.d),floordiv(y,ratio.d))),
          checked(std::max(ceildiv(x,ratio.d),ceildiv(y,ratio.d))));
}

// The returned p,q,h are ALREADY weighted. This certifies omega*Phi directly
// when a valid unweighted proposal slope exceeds the interval range.
inline LogPlane verify_old_weighted_log_plane(Rat d,Rat a,Rat w,Rat p,Rat q,I weight){
 if(weight.lo<0)throw std::runtime_error("negative log-plane weight");
 if(d.n==0||weight.hi==0)return {I(),I(),I()};
 if(d.n<0||d.d<=d.n||a.n<=0||w.n<0||w.d<=w.n||
    !((d+a*w)<Rat(1)))throw std::runtime_error("weighted log-plane triangle");
 if(!(a.n==a.d || i128(a.n)*d.d==i128(a.d)*(d.d-d.n)))
  throw std::runtime_error("weighted log-plane unsupported a");
 I D=d.box(),A=a.box(),W=w.box(),P=interval_times_rat(weight,p),
   Qq=interval_times_rat(weight,q),one=I::integer(1);
 auto proposal_derivative=[&](double u,int edge){
  double dd=d.dbl(),aa=a.dbl(),ww=w.dbl(),pp=p.dbl(),qq=q.dbl();
  if(edge==0)return dd/(1-dd*u)-dd-pp;
  if(edge==1)return -aa*ww/(1-ww*u)+(dd+aa*ww)/(1-(dd+aa*ww)*u)-dd-pp-qq*ww;
  return -aa*ww/(1-ww*u)+aa*ww/(1-dd-aa*ww*u)-qq*ww;
 };
 int64_t lower=GUARD-1;
 for(int edge=0;edge<3;++edge){
  double lo=0,hi=1,u;
  if(proposal_derivative(0,edge)>=0)u=0;
  else if(proposal_derivative(1,edge)<=0)u=1;
  else {for(int k=0;k<60;++k){u=(lo+hi)/2;if(proposal_derivative(u,edge)<0)lo=u;else hi=u;}u=(lo+hi)/2;}
  constexpr int64_t UQ=int64_t(1)<<30;
  Rat ur(std::clamp<int64_t>(std::llround(u*UQ),0,UQ),UQ);I U=ur.box();
  I X=edge==2?one:U,Y=edge==0?I():W*U,Z=one-D*X-A*Y;
  I phi=weight*(A*log_box(one-Y)-log_box(Z)-D*X)-P*X-Qq*Y;
  I derivative;
  if(edge==0)derivative=weight*D/Z-weight*D-P;
  else if(edge==1)derivative=weight*(D+A*W)/Z-weight*A*W/(one-Y)-weight*D-P-Qq*W;
  else derivative=weight*A*W/Z-weight*A*W/(one-Y)-Qq*W;
  I edgeLower=phi-U*derivative+meet_min(I(),derivative);
  lower=std::min(lower,edgeLower.lo);
 }
 return {P,Qq,I(lower,lower)};
}


