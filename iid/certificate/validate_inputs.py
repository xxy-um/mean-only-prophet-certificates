"""Validate every frozen Bellman coefficient by exact rational enclosures.

No optimizer or platform exp/log/cbrt implementation is used.  The input data
are the frozen integer slope/rate/state knots and the hexadecimal header.
Usage: python3 validate_inputs.py [header] [rates] [knots]
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re,sys,time,os
if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
    raise SystemExit('Assertions must remain enabled: do not use python -O or PYTHONOPTIMIZE.')
HERE=Path(__file__).resolve().parent/'data'
header=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'coeff_N640_Q131093.h'
rates_file=Path(sys.argv[2]) if len(sys.argv)>2 else HERE/'rates_refined_Q131093.json'
grid_file=Path(sys.argv[3]) if len(sys.argv)>3 else HERE/'grid1200.json'
text=header.read_text(); data=json.loads((HERE/'threshold40.json').read_text())
N,Q=map(int,re.search(r'N=(\d+),Q=(\d+)',text).groups())
assert N>=320 and N%data['steps']==0
eta=F(1,2**48);U=2**40;sub=N//data['steps'];a=[F(v,data['a_denominator']) for v in data['a_numerators']]
assert len(a)==data['steps']==40
P=sum(a)/40;b=[v/(1+P) for v in a];height=P/(1+P)
threshold=[sum(b[k//sub:])/40-F(k%sub,N)*b[k//sub] for k in range(N)]+[F(0)]
cap=max(b[j]+threshold[j*sub]-1 for j in range(40))
assert F(85,100)<cap<F(9,10) and all(F(1,4)<v<2 for v in b)
rates=json.loads(rates_file.read_text()); rr=[F(v,rates['denominator']) for v in rates['numerators']]
assert Q==len(rr)==rates['Q'] and rr[0]==0 and all(v<w for v,w in zip(rr,rr[1:])) and rr[-1]<8193
knots=json.loads(grid_file.read_text());kk=knots['numerators'];assert len(kk)==knots['m']+1 and kk[0]==0 and kk[-1]==knots['denominator']==2**30 and all(v<w for v,w in zip(kk,kk[1:]))
grid_header=grid_file.with_suffix('.h').read_text()
assert int(re.search(r'GRID_M=(\d+)',grid_header).group(1))==knots['m']
assert [int(v) for v in re.search(r'GRID_KNOTS\[\]=\{([^}]*)\}',grid_header).group(1).split(',')]==kk
# Header constants are converted exactly from their finite hexadecimal strings.
def exacthex(s):
    return F.from_float(float.fromhex(s))
arrays={name:[exacthex(v) for v in seq.split(',')] for name,seq in re.findall(r'const double (\w+)\[\]=\{([^}]*)\};',text)}
constants={name:exacthex(v) for name,v in re.findall(r'(dt|cap|H)=(0x[^,;]+)',text)}
assert abs(constants['dt']-F(1,N))<eta/N
assert int(re.search(r'HU=(\d+)LL',text).group(1))==(height*U).__floor__()
maxerr=F(0)
def verify(v,lo,hi):
    global maxerr
    error=max(abs(v-lo),abs(v-hi));assert error<eta
    maxerr=max(maxerr,error)
def exp_bounds(z):
    assert z>=0
    if z>=64:return F(0),F(1,2**int(z))
    k=0
    while z>1:z/=2;k+=1
    term=s=F(1)
    for j in range(1,26):
        term*=-z/j;s+=term
        if j==24:hi=s
        if j==25:lo=s
    for _ in range(k):lo*=lo;hi*=hi
    return lo,hi
start=time.monotonic()
for j,r in enumerate(rr):
    elo,ehi=exp_bounds(r/N);rlo,rhi=exp_bounds(r);zlo,zhi=exp_bounds((2+F(1,N))*r)
    bounds={'rr':(r,r),'ee':(elo,ehi),'ii':((1-ehi)/r,(1-elo)/r) if r else (F(1,N),F(1,N)),'vv':(r+rlo-1,r+rhi-1),'er':(rlo,rhi),'e2':(zlo,zhi)}
    for key,(lo,hi) in bounds.items():verify(arrays[key][j],lo,hi)
for key in ['rr','ee','ii','vv','er','e2']:
    assert len(arrays[key])==Q+1 and arrays[key][-1]==0
for key,values in [('bs',[b[k//sub] for k in range(N)]),('bl',[b[k//sub] for k in range(N)]),('bp',[F(0)]*N),('DD',[1-v for v in threshold[:-1]])]:
    assert len(arrays[key])==N
    for supplied,true in zip(arrays[key],values):verify(supplied,true,true)
verify(constants['H'],height,height);verify(constants['cap'],cap,cap)
print(json.dumps({'validation':'PASS: exact rational coefficient and knot checks','N':N,'Q':Q,'m':knots['m'],'coefficient_max_error':float(maxerr),'coefficient_allowance':'2^-48','H':str(height),'cap':str(cap),'header_sha256':hashlib.sha256(header.read_bytes()).hexdigest(),'seconds':time.monotonic()-start},indent=2))
