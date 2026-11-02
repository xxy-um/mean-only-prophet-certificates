"""Rigorous selected literal FRLP row check using rational interval arithmetic.

Usage: python3 audit/check_exact_row.py candidate.json
The JSON object (or list) specifies k, i, j, leaf, post_integer, table_path.
A relative table_path is resolved relative to the candidate JSON file.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,bisect,array
BASE=Path(__file__).resolve().parents[1]/'data'
SCALE=2**200
class Iv:
 def __init__(self,lo,hi=None):
  lo=F(lo);hi=F(lo if hi is None else hi)
  self.lo=F((lo*SCALE).__floor__(),SCALE);self.hi=F((hi*SCALE).__ceil__(),SCALE);assert self.lo<=self.hi
 def __add__(a,b):
  if not isinstance(b,Iv):b=Iv(b)
  return Iv(a.lo+b.lo,a.hi+b.hi)
 __radd__=__add__
 def __neg__(a):return Iv(-a.hi,-a.lo)
 def __sub__(a,b):return a+-asiv(b)
 def __rsub__(a,b):return asiv(b)+-a
 def __mul__(a,b):
  b=asiv(b);v=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi];return Iv(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(a,b):
  b=asiv(b);assert not b.lo<=0<=b.hi
  return a*Iv(1/b.hi,1/b.lo)
 def __rtruediv__(a,b):return asiv(b)/a
 def __pow__(a,n):
  assert n>=0;v=Iv(1)
  for _ in range(n):v=v*a
  return v
 def desc(a):return {'lo':float(a.lo),'hi':float(a.hi),'width':float(a.hi-a.lo)}
def asiv(x):return x if isinstance(x,Iv) else Iv(x)
def expneg(z):
 z=F(z);assert z>=0
 if z>=512:return Iv(0,F(1,2**int(z)))
 k=0
 while z>1:z/=2;k+=1
 term=s=F(1)
 for j in range(1,81):
  term*=-z/j;s+=term
  if j==79:lo=s
  if j==80:hi=s
 ans=Iv(lo,hi)
 for _ in range(k):ans=Iv(ans.lo*ans.lo,ans.hi*ans.hi)
 return ans
LOGCACHE={}
def logf(z):
 z=F(z);assert z>0
 original=z
 if z in LOGCACHE:return LOGCACHE[z]
 if z==1:return Iv(0)
 if z<1:return -logf(1/z)
 k=0
 while z>2:z/=2;k+=1
 w=(z-1)/(z+1);v=w;s=F(0)
 n=80
 for j in range(n):s+=v/(2*j+1);v*=w*w
 # positive atanh remainder <= w^(2n+1)/((2n+1)(1-w²))
 tail=2*v/((2*n+1)*(1-w*w))
 ans=Iv(2*s,2*s+tail)
 if k:ans+=k*logf(F(2))
 LOGCACHE[original]=ans
 return ans
def logiv(z):return Iv(logf(z.lo).lo,logf(z.hi).hi)
def maxiv(vals):return Iv(max(v.lo for v in vals),max(v.hi for v in vals))
def miniv(a,b):return Iv(min(a.lo,b.lo),min(a.hi,b.hi))
def main(cand):
 d=json.loads((BASE/'threshold40.json').read_text());a=[F(v,d['a_denominator']) for v in d['a_numerators']];avg=sum(a)/40;b=[v/(1+avg) for v in a]
 N=640;dt=F(1,N);k=cand['k'];t=F(k,N);s=b[k//16];theta=sum(b[k//16:])/40-F(k%16,N)*s;cap=F(d['cap'])
 gj=json.loads((BASE/'grid1200.json').read_text());g=[F(v,gj['denominator']) for v in gj['numerators']];i,j=cand['i'],cand['j'];y,x=g[i],g[j]
 rr=json.loads((BASE/'rates_refined_Q131093.json').read_text());rates=[F(v,rr['denominator']) for v in rr['numerators']];leaf=cand['leaf'];rl=rates[leaf];rh=rates[leaf+1] if leaf+1<len(rates) else None
 pp=cap*x+y*(1-theta-s);u=s*(1-y+t*y)
 def EI(r):
  e=expneg(dt*r);ii=(1-e)/r if r else Iv(dt)
  return e,ii
 ee,ii=EI(rl);phi=1-y+y*(t+ii+ee*(1-t-dt));yp=y*ee/phi
 num=cap*x+y*s*dt*(rl+expneg(rl)-1);clipped=num.lo>=cap*y;xp=ee*num/(cap*phi);xp=miniv(xp,yp)
 costs=[Iv(min(pp*dt,0))]
 if rh is not None:
  ehi,ihi=EI(rh);alpha=(y*s*ehi+max(-pp,0)*ihi/dt)*dt*dt;beta=dt*(u+y*s*ihi);base=miniv(Iv(pp*dt),pp*ihi)
  costs += [base+beta-beta*beta/(4*alpha),base+alpha/4*(3+2*logiv(beta/alpha))]
  z=beta*expneg(rh)
  if z.lo>alpha.hi/2:costs.append(base+z+alpha*rh/2-alpha*alpha/(24*(z-alpha/2)))
  elif z.hi>alpha.lo/2:raise ValueError('shifted guard unresolved')
 cost=maxiv(costs)
 # table read
 tab=array.array('q');file=Path(cand['table_path'])
 with file.open('rb') as f:tab.fromfile(f,file.stat().st_size//8)
 def tv(i,j):return F(tab[i*(i+1)//2+j],2**40)
 # Find cell containing entire exact interval
 if clipped or xp.lo>=yp.hi:
  interp=Iv(0);vertices=[]
 else:
  Y=(yp.lo+yp.hi)/2;X=(xp.lo+xp.hi)/2;ci=bisect.bisect_left(g,Y)-1;cj=bisect.bisect_right(g,X)-1;ci=max(0,min(1199,ci));cj=max(0,min(ci,cj))
  assert g[ci]<=yp.lo<=yp.hi<=g[ci+1],('ycell',ci,yp.desc())
  assert g[cj]<=xp.lo<=xp.hi<=g[cj+1],('xcell',cj,xp.desc())
  dy=g[ci+1]-g[ci];dx=g[cj+1]-g[cj];ay=yp-g[ci];bx=xp-g[cj];cross=bx*dy-ay*dx
  if cross.hi<=0:
   weights=[(dy-ay)/dy,(ay*dx-bx*dy)/(dx*dy),bx/dx];vertices=[(ci,cj),(ci+1,cj),(ci+1,cj+1)]
  elif cross.lo>=0:
   weights=[(dx-bx)/dx,(bx*dy-ay*dx)/(dx*dy),ay/dy];vertices=[(ci,cj),(ci,cj+1),(ci+1,cj+1)]
  else:raise ValueError('triangle unresolved')
  interp=sum([w*tv(*v) for w,v in zip(weights,vertices)],Iv(0))
 rhs=cost+phi*interp;post=F(cand['post_integer'],2**40);gap=Iv(post)-rhs
 out={'candidate':cand,'rate_interval':[str(rl),str(rh)],'theta':str(theta),'slope':str(s),'p':str(pp),'phi':phi.desc(),'successor_y':yp.desc(),'successor_x':xp.desc(),'vertices':vertices,'vertex_integers':[tab[ii*(ii+1)//2+jj] for ii,jj in vertices],'cost':cost.desc(),'rhs':rhs.desc(),'post':float(post),'post_minus_rhs':gap.desc(),'VIOLATION':gap.lo>0,'rational_rhs_upper':str(rhs.hi),'rational_post':str(post)}
 print(json.dumps(out,indent=2))
 return out
if __name__=='__main__':
 f=Path(sys.argv[1]).resolve();d=json.loads(f.read_text())
 rows=d if isinstance(d,list) else [d]
 for c in rows:
  table=Path(c['table_path'])
  if not table.is_absolute():c['table_path']=str((f.parent/table).resolve())
  main(c)
