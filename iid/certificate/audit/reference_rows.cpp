// Independent literal Appendix A.1 row scan, for diagnostics only.
// Long-double calculations over the validated frozen inputs; not the final proof checker.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#ifndef COEFFICIENT_FILE
#define COEFFICIENT_FILE "coeff_N640_Q131093.h"
#endif
#include COEFFICIENT_FILE
using I=std::int64_t; using R=long double;
constexpr I U=I(1)<<40,L=I(1)<<30;
std::vector<I> knots,values;
int m;
std::size_t idx(int i,int j){return std::size_t(i)*(i+1)/2+j;}
R interpolate(R y,R x){
 if(x>=y)return 0;
 R Y=y*L,X=x*L;
 int i=int(std::lower_bound(knots.begin(),knots.end(),Y)-knots.begin())-1;
 int j=int(std::upper_bound(knots.begin(),knots.end(),X)-knots.begin())-1;
 i=std::clamp(i,0,m-1);j=std::clamp(j,0,i);
 R a=(Y-knots[i])/(knots[i+1]-knots[i]);
 R b=(X-knots[j])/(knots[j+1]-knots[j]);
 if(b<=a)return ((1-a)*values[idx(i,j)]+(a-b)*values[idx(i+1,j)]+b*values[idx(i+1,j+1)])/U;
 return ((1-b)*values[idx(i,j)]+(b-a)*values[idx(i,j+1)]+a*values[idx(i+1,j+1)])/U;
}
R row(int k,int i,int j,int leaf){
 R d=R(1)/N,t=R(k)/N,b=bs[k],c=cap;
 R y=R(knots[i])/L,x=R(knots[j])/L;
 R p=c*x+y*(R(DD[k])-b),u=b*(1-y+t*y);
 R E=ee[leaf],Ival=ii[leaf];
 R phi=1-y+y*(t+Ival+E*(1-t-d));
 R yp=y*E/phi;
 R xp=std::min(yp,E*(c*x+y*b*d*R(vv[leaf]))/(c*phi));
 R cost=std::min(p*d,R(0));
 int h=leaf+1;
 if(h<Q){
  R alpha=(y*b*R(ee[h])+std::max(-p,R(0))*R(ii[h])/d)*d*d;
  R beta=d*(u+y*b*R(ii[h]));
  R base=std::min(p*d,p*R(ii[h]));
  cost=std::max(cost,base+beta-beta*beta/(4*alpha));
  cost=std::max(cost,base+alpha/4*(3+2*std::log(beta/alpha)));
  R z=beta*R(er[h]);
  if(z>alpha/2)cost=std::max(cost,base+z+alpha*R(rr[h])/2-alpha*alpha/(24*(z-alpha/2)));
 }
 return cost+phi*interpolate(yp,xp);
}
int main(int argc,char **argv){
 assert(argc==4);
 std::ifstream g(argv[1]);g>>m;knots.resize(m+1);for(auto &v:knots)g>>v;assert(g);
 values.resize((m+1)*(m+2)/2);std::ifstream f(argv[2],std::ios::binary);f.read(reinterpret_cast<char*>(values.data()),values.size()*sizeof(I));assert(f);
 int k=std::atoi(argv[3]),i,j;I post;
 while(std::cin>>i>>j>>post){
  R best=0;int leaf=0;
  for(int l=0;l<Q;l++){R rhs=row(k,i,j,l);if(rhs<best){best=rhs;leaf=l;}}
  std::cout<<std::setprecision(24)<<"{\"k\":"<<k<<",\"i\":"<<i<<",\"j\":"<<j<<",\"leaf\":"<<leaf<<",\"post_integer\":"<<post<<",\"rhs\":"<<best<<",\"post_minus_rhs\":"<<R(post)/U-best<<"}\n";
 }
}
