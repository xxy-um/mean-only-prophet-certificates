"""Exact selected-row audit of the literal Appendix A.1 FRLP.
Uses the prior interval checker, with the exact r=0 successor identity handled
analytically to avoid an interval straddling an unchanged grid vertex.
"""
from pathlib import Path
from fractions import Fraction as F
import json,array,sys,importlib.util,contextlib,io
sp=importlib.util.spec_from_file_location('exact_row',str(Path(__file__).resolve().parent/'check_exact_row.py'))
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
def exact_zero(c):
 d=json.loads((m.BASE/'threshold40.json').read_text());a=[F(v,d['a_denominator']) for v in d['a_numerators']];avg=sum(a)/40;b=[v/(1+avg) for v in a]
 k=c['k'];dt=F(1,640);t=F(k,640);s=b[k//16];theta=sum(b[k//16:])/40-F(k%16,640)*s;cap=F(d['cap'])
 grid=json.loads((m.BASE/'grid1200.json').read_text());y=F(grid['numerators'][c['i']],grid['denominator']);x=F(grid['numerators'][c['j']],grid['denominator'])
 rr=json.loads((m.BASE/'rates_refined_Q131093.json').read_text());rh=F(rr['numerators'][1],rr['denominator'])
 p=cap*x+y*(1-theta-s);u=s*(1-y+t*y);ee=m.expneg(dt*rh);ii=(1-ee)/rh
 alpha=(y*s*ee+max(-p,0)*ii/dt)*dt*dt;beta=dt*(u+y*s*ii);base=m.miniv(m.Iv(p*dt),p*ii)
 costs=[m.Iv(min(p*dt,0)),base+beta-beta*beta/(4*alpha),base+alpha/4*(3+2*m.logiv(beta/alpha))]
 z=beta*m.expneg(rh)
 if z.lo>alpha.hi/2:costs.append(base+z+alpha*rh/2-alpha*alpha/(24*(z-alpha/2)))
 elif z.hi>alpha.lo/2:raise ValueError('Unresolved shifted condition')
 cost=m.maxiv(costs);tab=array.array('q');tab.frombytes(Path(c['table_path']).read_bytes());vnext=F(tab[c['i']*(c['i']+1)//2+c['j']],2**40)
 rhs=cost+m.Iv(vnext);post=F(c['post_integer'],2**40);gap=m.Iv(post)-rhs
 return {'candidate':c,'method':'exact r=0 identity: Phi=1 and successor=(y,x)','rate_interval':['0',str(rh)],'cost':cost.desc(),'next_grid_integer':int(vnext*2**40),'rhs':rhs.desc(),'post':float(post),'post_minus_rhs':gap.desc(),'VIOLATION':gap.lo>0,'strictly_feasible':gap.hi<0,'rational_rhs_upper':str(rhs.hi),'rational_post':str(post)}
def run(results_path,table_path,outpath):
 rows=[json.loads(v) for v in Path(results_path).read_text().splitlines()];rows=sorted(rows,key=lambda r:r['post_minus_rhs'],reverse=True)
 selected=rows[:3];target=next(r for r in rows if (r['i'],r['j'])==(1200,0))
 if target not in selected:selected.append(target)
 output=[]
 for r in selected:
  c={key:r[key] for key in ['k','i','j','leaf','post_integer']};c['table_path']=str(Path(table_path).resolve())
  if c['leaf']==0:o=exact_zero(c)
  else:
   with contextlib.redirect_stdout(io.StringIO()):o=m.main(c)
   o['strictly_feasible']=o['post_minus_rhs']['hi']<0
  output.append(o);print(c['k'],c['i'],c['j'],c['leaf'],o['post_minus_rhs'],'strictly_feasible',o['strictly_feasible'])
 Path(outpath).write_text(json.dumps(output,indent=2))
if __name__=='__main__':run(*sys.argv[1:])
