"""Actual paired resolvent evaluator; every unclosed error is kept explicit."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import json,sys,traceback,time
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import sympy as sp
import current_full_domain_contact_assembly as prior
root=Path(__file__).resolve().parents[1]
STEM='current-signed-resolvent-projection';PREVIOUS=prior.STEM
OWN=['notes/'+STEM+'.md','checks/current_signed_resolvent_projection.py','evidence/'+STEM+'-review.json']
OWN.append('checks/current_signed_resolvent_certificate.py')
ROWS=[];NUMERICS=[];encoded=prior.encoded;decode=prior.decode
t=sp.Symbol('t');a=sp.symbols('a0:15');v=sp.symbols('v0:16');variables=(t,*a,*v)

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS},name
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name+': '+str(detail))

prior.check=check;prior.gy.check=check;prior.bh.prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def hull(lo,hi):return arb(lo.lower()).union(arb(hi.upper()))

def Dgeo(expr):
 return sp.expand(sum((sp.diff(expr,a[i])*a[i+1] for i in range(len(a)-1)),sp.S.Zero)
                 +sum((sp.diff(expr,v[i])*v[i+1] for i in range(len(v)-1)),sp.S.Zero))
def BB(expr,j):return sp.expand(Dgeo(expr)+a[0]*((j+2)*t-1)*expr-2*a[0]*t*(1-t)*sp.diff(expr,t))

def coefficients(N):
 P=[sp.S.Zero,a[0]/4];V=[sp.S.Zero,(v[1]+a[0]*(3*t-2)*v[0])/4]
 for j in range(1,N):
  ss=sum((P[k]*P[j-k] for k in range(1,j)),sp.S.Zero)
  sv=sum((V[k]*P[j-k]+P[k]*V[j-k] for k in range(1,j)),sp.S.Zero)
  P.append(sp.expand(-BB(P[j],j)/2+a[0]*t*(1-t)*ss/4))
  V.append(sp.expand(-BB(V[j],j)/2+(1-t)*v[0]*BB(P[j],j)/2
    +t*(1-t)*((v[1]+a[0]*(3*t-2)*v[0])*ss+a[0]*sv)/4))
  print(json.dumps(dict(event='radial_coefficient',order=j+1,P_terms=len(sp.Add.make_args(P[-1])),V_terms=len(sp.Add.make_args(V[-1])))),flush=True)
 return P,V

def local_identities(P,V):
 # Formal finite series; the rational remainder is retained separately.
 NN=6
 def plus(x,y):return [sp.expand(xx+yy) for xx,yy in zip(x,y)]
 def sc(x,c):return [sp.expand(c*xx) for xx in x]
 def mul(x,y):return [sp.expand(sum((x[j]*y[k-j] for j in range(k+1)),sp.S.Zero)) for k in range(NN+1)]
 A=P[:NN+1];B=[(-1)**j*c for j,c in enumerate(A)]
 AV=V[:NN+1];BV=[(-1)**j*c for j,c in enumerate(AV)]
 Q=mul(A,B);QV=plus(mul(AV,B),mul(A,BV));S=plus(A,B);SV=plus(AV,BV)
 N=plus(sc(Q,4*(1-t)),sc(S,2))
 NV=plus(sc(plus(QV,sc(Q,v[0]*(1+t))),4*(1-t)),sc(plus(SV,sc(S,v[0]*t)),2))
 q=sc(Q,t*(1-t));qv=sc(QV,t*(1-t));qq=mul(q,q)
 one=[sp.S.One]+[sp.S.Zero]*NN
 inv=plus(plus(one,sc(q,-1)),qq);inv2=plus(plus(one,sc(q,-2)),sc(qq,3))
 H=plus(mul(NV,inv),sc(mul(mul(N,qv),inv2),-1))
 check('local_thermal_polynomial_degree',all(sp.degree(H[j],t)<=12 for j in [2,4,6]))
 coeffs=[]
 for k in [1,2,3]:
  hh=sp.Poly(H[2*k],t);val=sp.S.Zero
  for (j,),c in hh.terms():
   moment=sum((sp.Rational((-1)**i*comb(k-1,i),2*j+2*i+3) for i in range(k)),sp.S.Zero)
   val+=c*moment
  coeffs.append(sp.expand(val))
 Delta=lambda f:sp.expand(Dgeo(Dgeo(f))+a[0]*Dgeo(f))
 R=-2*(a[1]+a[0]**2);VV=v[0];DV=Delta(VV);LL=lambda f:Delta(f)+R*f
 O4=sp.expand(Delta(DV)+R*DV+Delta(R*VV)+R**2*VV/2)
 O6=sp.expand(-(LL(Delta(LL(VV)))+Delta(R)*LL(VV)+Delta(VV*Delta(R)))/105
  -(16*R**3*VV+24*R**2*DV+24*Delta(R**2*VV)+24*Delta(R*DV))/1008)
 check('continuous_frequency_curvature_conversion',sp.expand(coeffs[0]+DV/6)==0)
 check('continuous_frequency_geometric_a2',sp.expand(coeffs[1]+O4/60)==0)
 check('continuous_frequency_geometric_a6',sp.expand(coeffs[2]-O6/4)==0)
 NUMERICS.append(dict(kind='exact_continuous_frequency_contact_identities',
  coefficients=[str(c) for c in coeffs],geometric_a2=str(O4),geometric_a6=str(O6),
  physical_frequency_sum_not_replaced_by_integral=True,finite_temperature_error_still_required=True))

def bernstein(cs):
 n=len(cs)-1
 return [sum((cs[k]*arb(comb(j,k))/comb(n,k) for k in range(j+1)),arb(0)) for j in range(n+1)]
def poly_range(cs,lo,hi):
 n=len(cs)-1;d=hi-lo
 shifted=[sum((cs[k]*comb(k,j)*lo**(k-j)*d**j for k in range(j,n+1)),arb(0)) for j in range(n+1)]
 bs=bernstein(shifted);return hull(min(c.lower() for c in bs),max(c.upper() for c in bs))

def jet_mul(x,y,n):return [sum((x[j]*y[k-j] for j in range(k+1)),arb(0)) for k in range(n)]
def jet_D(x,F):
 der=[(j+1)*x[j+1] for j in range(len(x)-1)]
 return [-c for c in jet_mul(der,F,len(der))]
def geometry_jets(p,lo,hi,order):
 beta,u,ep,tau=p['beta'],p['u'],p['eps'],p['tau']
 # Retain the exact normalization 6 beta+u+epsilon=1 before interval evaluation.
 cs=[1-3*tau,-8*beta-2*u-4*ep+tau,3*beta+u+6*ep,-4*ep,ep]
 yy=hull(lo,hi);tj=[];n=order+3
 for k in range(n):
  normal=poly_range([c*(-j)**k for j,c in enumerate(cs)],lo,hi)
  thermal=-tau*2**k/yy**2+3*tau/yy
  tj.append((normal+thermal)/factorial(k))
 assert tj[0]>0
 F=[tj[0].sqrt()]
 for k in range(1,n):F.append((tj[k]-sum((F[j]*F[k-j] for j in range(1,k)),arb(0)))/(2*F[0]))
 aa=[F[j]-(j+1)*F[j+1] for j in range(n-1)]
 vv=[poly_range([arb(-1),arb(4)*(-1)**k,arb(-6)*(-2)**k,arb(4)*(-3)**k,arb(-1)*(-4)**k],lo,hi)/factorial(k) for k in range(n)]
 # The constant coefficient of v vanishes under x derivatives of positive order.
 for k in range(1,n):vv[k]+=arb(1)/factorial(k)
 AJ=[];VJ=[]
 for j in range(order+1):
  AJ.append(aa[0]);VJ.append(vv[0]);aa=jet_D(aa,F);vv=jet_D(vv,F)
 return AJ,VJ,F[0]

def la_add(x,y):
 out=dict(x)
 for k,c in y.items():out[k]=out.get(k,arb(0))+c
 return {k:c for k,c in out.items() if not c.is_zero()}
def la_mul(x,y):
 out={}
 for i,c in x.items():
  for j,d in y.items():out[i+j]=out.get(i+j,arb(0))+c*d
 return {k:c for k,c in out.items() if not c.is_zero()}
def la_scale(x,c):return {k:c*v for k,v in x.items() if not (c*v).is_zero()}
def la_D(x):return {k-1:k*c for k,c in x.items() if k}
def la_y(x):return {k+1:c for k,c in x.items()}
def la_compile(x):
 pos={k:c for k,c in x.items() if k>=0};neg={k:c for k,c in x.items() if k<0}
 if not pos:return (0,[arb(0)],neg)
 lo=min(pos);hi=max(pos)
 return lo,bernstein([pos.get(k,arb(0)) for k in range(lo,hi+1)]),neg
def split_bernstein(bs,x):
 rows=[list(bs)]
 while len(rows[-1])>1:rows.append([(1-x)*a+x*b for a,b in zip(rows[-1],rows[-1][1:])])
 return [r[0] for r in rows],[r[-1] for r in reversed(rows)]
def la_range(comp,lo,hi):
 power,bs,neg=comp
 left,_=split_bernstein(bs,hi)
 _,restricted=split_bernstein(left,lo/hi)
 yy=hull(lo,hi)
 out=hull(min(c.lower() for c in restricted),max(c.upper() for c in restricted))*yy**power
 for k,c in neg.items():out+=c*hull(hi**k,lo**k)
 return out
def geometry_numerators(p,order):
 beta,u,ep,tau=p['beta'],p['u'],p['eps'],p['tau']
 T={-2:-tau,-1:3*tau,0:1-3*tau,1:-8*beta-2*u-4*ep+tau,2:3*beta+u+6*ep,3:-4*ep,4:ep}
 DT=la_D(T);AA=[la_add(T,la_scale(la_y(DT),arb(1)/2))]
 BB=[{0:arb(-1),1:arb(4),2:arb(-6),3:arb(4),4:arb(-1)}]
 for j in range(order):
  AA.append(la_y(la_add(la_mul(T,la_D(AA[-1])),la_scale(la_mul(AA[-1],DT),-arb(j+1)/2))))
  BB.append(la_y(la_add(la_mul(T,la_D(BB[-1])),la_scale(la_mul(BB[-1],DT),-arb(j)/2))))
 return la_compile(T),[la_compile(c) for c in AA],[la_compile(c) for c in BB]
def rational_geometry_jets(comp,lo,hi):
 TT,AA,BB=comp;T=la_range(TT,lo,hi);assert T>0
 F=T.sqrt();AJ=[la_range(c,lo,hi)/F**(j+1) for j,c in enumerate(AA)]
 VJ=[la_range(c,lo,hi)/F**j for j,c in enumerate(BB)]
 return AJ,VJ,F

def compile_poly(expr):
 pp=sp.Poly(expr,*variables);out=[];degree=pp.degree(t)
 for powers,c in pp.terms():
  out.append((powers[0],arb(int(sp.numer(c)))/int(sp.denom(c)),tuple((i,powr) for i,powr in enumerate(powers[1:]) if powr)))
 return degree,out
def eval_poly(comp,values):
 degree,terms=comp;cs=[arb(0)]*(degree+1);cache={}
 for k,c,mm in terms:
  term=c
  for i,powr in mm:
   if (i,powr) not in cache:cache[i,powr]=values[i]**powr
   term*=cache[i,powr]
  cs[k]+=term
 bs=bernstein(cs)
 return arb(max(abs(b).upper() for b in bs))

def bi_add(x,y):
 out=dict(x)
 for k,c in y.items():out[k]=out.get(k,arb(0))+c
 return {k:c for k,c in out.items() if not c.is_zero()}
def bi_scale(x,c):return {k:c*v for k,v in x.items() if not (c*v).is_zero()}
def bi_mul(x,y):
 out={}
 for (i,j),c in x.items():
  for (k,l),d in y.items():out[i+k,j+l]=out.get((i+k,j+l),arb(0))+c*d
 return {k:c for k,c in out.items() if not c.is_zero()}
def bi_D(x,index):
 out={}
 for (i,j),c in x.items():
  deg=(i,j)[index]
  if deg:out[i-(index==0),j-(index==1)]=deg*c
 return out
def bi_shift(x,iy=0,it=0):return {(i+iy,j+it):c for (i,j),c in x.items()}
def bi_value(x,y,tt):return sum((c*y**i*tt**j for (i,j),c in x.items()),arb(0))
def joint_coefficients(p,N):
 beta,u,ep,tau=p['beta'],p['u'],p['eps'],p['tau']
 T={(k,0):c for k,c in {-2:-tau,-1:3*tau,0:1-3*tau,1:-8*beta-2*u-4*ep+tau,2:3*beta+u+6*ep,3:-4*ep,4:ep}.items()}
 Ty=bi_D(T,0);A0=bi_add(T,bi_scale(bi_shift(Ty,iy=1),arb(1)/2))
 vv={(k,0):arb(c) for k,c in enumerate([-1,4,-6,4,-1])}
 B1=bi_shift(bi_mul(T,bi_D(vv,0)),iy=1)
 tl={(0,1):arb(1),(0,2):arb(-1)};ell={(0,0):arb(1),(0,1):arb(-1)}
 Bv=bi_add(B1,bi_mul(bi_mul(A0,vv),{(0,1):arb(3),(0,0):arb(-2)}))
 def op(x,j):
  radial=bi_shift(bi_add(bi_mul(T,bi_D(x,0)),bi_scale(bi_mul(x,Ty),-arb(j)/2)),iy=1)
  transport=bi_mul(bi_mul(A0,x),{(0,1):arb(j+2),(0,0):arb(-1)})
  return bi_add(bi_add(radial,transport),bi_scale(bi_mul(bi_mul(A0,tl),bi_D(x,1)),-2))
 P=[{},bi_scale(A0,arb(1)/4)];V=[{},bi_scale(Bv,arb(1)/4)]
 for j in range(1,N):
  ss={};sv={}
  for k in range(1,j):
   ss=bi_add(ss,bi_mul(P[k],P[j-k]));sv=bi_add(sv,bi_scale(bi_mul(V[k],P[j-k]),2))
  dp=op(P[j],j);dv=op(V[j],j)
  P.append(bi_add(bi_scale(dp,-arb(1)/2),bi_scale(bi_mul(bi_mul(A0,tl),ss),arb(1)/4)))
  V.append(bi_add(bi_add(bi_scale(dv,-arb(1)/2),bi_scale(bi_mul(bi_mul(ell,vv),dp),arb(1)/2)),
   bi_scale(bi_mul(tl,bi_add(bi_mul(Bv,ss),bi_mul(A0,sv))),arb(1)/4)))
  print(json.dumps(dict(event='joint_radial_coefficient',order=j+1,P_terms=len(P[-1]),V_terms=len(V[-1]))),flush=True)
 return T,P,V,op

def la_split(bs,x):
 rows=[list(bs)]
 while len(rows[-1])>1:rows.append([la_add(la_scale(a,1-x),la_scale(b,x)) for a,b in zip(rows[-1],rows[-1][1:])])
 return [r[0] for r in rows],[r[-1] for r in reversed(rows)]
def t_control_polynomials(poly):
 degree=max(j for i,j in poly);rows=[]
 for j in range(degree+1):
  row={}
  for (iy,it),c in poly.items():
   if it<=j:row[iy]=row.get(iy,arb(0))+c*arb(comb(j,it))/comb(degree,it)
  rows.append(row)
 segments=[rows]
 for _ in range(2):segments=[half for segment in segments for half in la_split(segment,arb(1)/2)]
 return [la_compile(row) for segment in segments for row in segment]
def positive_y_tree(bs):
 values=[];remaining=bs
 def leaves(row,level):
  if level==0:return [arb(max(abs(b).upper() for b in row))]
  left,right=split_bernstein(row,arb(1)/2)
  return leaves(left,level-1)+leaves(right,level-1)
 for k in range(12):
  remaining,band=split_bernstein(remaining,arb(1)/2)
  values.extend(leaves(band,3))
 values.append(arb(max(abs(b).upper() for b in remaining)))
 return values
def joint_geometry_cover(p,P,V,N,joint):
 T,PN,VN,op=joint;cover=[]
 for k in range(12):
  lo=arb(2)**(-k-1);hi=arb(2)**-k
  for j in range(8):cover.append((lo+(hi-lo)*j/8,lo+(hi-lo)*(j+1)/8))
 cover.append((arb(2)**-56,arb(2)**-12))
 geo=geometry_numerators(p,N+1);georows=[];a_bounds=[arb(0)]*(N+2);v_bounds=[arb(0)]*(N+2)
 for lo,hi in cover:
  AJ,VJ,F=rational_geometry_jets(geo,lo,hi)
  for k in range(N+2):a_bounds[k]=arb(max(a_bounds[k],abs(AJ[k]).upper()));v_bounds[k]=arb(max(v_bounds[k],abs(VJ[k]).upper()))
  georows.append((AJ,VJ,F))
 power_cache={}
 def power(i,k):
  if (i,k) not in power_cache:
   lo,hi=cover[i];power_cache[i,k]=hull(hi**k,lo**k) if k<0 else hull(lo**k,hi**k)
  return power_cache[i,k]
 def bound(poly,weight):
  out=[arb(0)]*len(cover)
  for factor,bs,neg in t_control_polynomials(poly):
   ranges=positive_y_tree(bs)
   for i,row in enumerate(ranges):
    val=row*arb(power(i,factor).upper())
    for k,c in neg.items():val+=abs(c)*arb(power(i,k).upper())
    val=arb(val.upper())/arb(georows[i][2].lower())**weight
    out[i]=arb(max(out[i],val.upper()))
  return out
 PR=[None];VR=[None]
 for j in range(1,N+1):
  PR.append(bound(PN[j],j));VR.append(bound(VN[j],j))
  print(json.dumps(dict(event='joint_complete_coefficient_bound',N=N,j=j,P=str(max(PR[-1])),V=str(max(VR[-1])))),flush=True)
 RR=bound(op(PN[N],N),N+1);RV=bound(op(VN[N],N),N+1)
 PB=[arb(0)]+[arb(max(row)) for row in PR[1:]];VB=[arb(0)]+[arb(max(row)) for row in VR[1:]]
 L=arb(0);cells=[];yc=arb(2)**-49
 for i,(lo,hi) in enumerate(cover):
  ll=arb(max(lo,yc));F=georows[i][2]
  if hi>ll:L+=(hi/ll).log()/arb(F.lower())
  cells.append(dict(y_lo=encoded(lo),y_hi=encoded(hi),F=encoded(F),
   P_bounds=[encoded(arb(0))]+[encoded(row[i]) for row in PR[1:]],
   V_bounds=[encoded(arb(0))]+[encoded(row[i]) for row in VR[1:]],D_P_N=encoded(RR[i]),D_V_N=encoded(RV[i])))
 # Independent abstract-jet and joint-numerator values at exact rational points.
 for y,tt in [(arb(1)/2,arb(1)/4),(arb(3)/4,arb(3)/4),(arb(1)/32,arb(1)/2)]:
  AJ,VJ,F=rational_geometry_jets(geo,y,y);vals=AJ+[arb(0)]*(len(a)-len(AJ))+VJ+[arb(0)]*(len(v)-len(VJ))
  def val(expr):
   degree,terms=compile_poly(expr);out=arb(0)
   for k,c,mm in terms:
    x=c*tt**k
    for ii,powr in mm:x*=vals[ii]**powr
    out+=x
   return out
  for j in [1,3,6,N]:
   check('joint_numerator_identity_'+str(N)+'_'+str(j)+'_'+str(y.mid()),
    val(P[j]).overlaps(bi_value(PN[j],y,tt)/F**j) and val(V[j]).overlaps(bi_value(VN[j],y,tt)/F**j))
 NUMERICS.append(dict(kind='complete_joint_radial_coefficient_bounds',N=N,
  P=[encoded(x) for x in PB],V=[encoded(x) for x in VB],D_P_N=encoded(max(RR)),D_V_N=encoded(max(RV)),
  a_jets=[encoded(x) for x in a_bounds],v_jets=[encoded(x) for x in v_bounds],proper_bulk_length_upper=encoded(L),
  physical_inner_start=encoded(arb(2)**-56),positive_polynomial_tail_interval_includes_zero=True,
  negative_thermal_powers_keep_positive_endpoint=True,complete_t_cover=[0,1],thermal_p_range='0..infinity',cells=cells))
 return PB,VB,arb(max(RR)),arb(max(RV)),a_bounds,v_bounds,arb(max(row[2].upper() for row in georows))

def geometry_bounds(p,P,V,N,subdivisions=8):
 cp=[None]+[compile_poly(c) for c in P[1:]];cv=[None]+[compile_poly(c) for c in V[1:]]
 dr=compile_poly(BB(P[N],N));dv=compile_poly(BB(V[N],N))
 PB=[arb(0)]*(N+1);VB=[arb(0)]*(N+1);RB=arb(0);RVB=arb(0)
 ageo=[arb(0)]*(N+2);vgeo=[arb(0)]*(N+2);cells=[];Fmax=arb(0)
 geometry=geometry_numerators(p,N+1)
 for k in range(56):
  bot=arb(2)**(-k-1);top=arb(2)**-k
  for j in range(subdivisions):
   lo=bot+(top-bot)*j/subdivisions;hi=bot+(top-bot)*(j+1)/subdivisions
   AJ,VJ,F=rational_geometry_jets(geometry,lo,hi)
   if j==0 and k in [0,1,5,20,55]:
    independent_a,independent_v,independent_f=geometry_jets(p,(lo+hi)/2,(lo+hi)/2,N+1)
    check('independent_geometry_jet_'+str(N)+'_'+str(k),all(x.overlaps(y) for x,y in zip(AJ,independent_a)) and
      all(x.overlaps(y) for x,y in zip(VJ,independent_v)) and F.overlaps(independent_f))
   values=AJ+[arb(0)]*(len(a)-len(AJ))+VJ+[arb(0)]*(len(v)-len(VJ))
   for i in range(N+2):ageo[i]=arb(max(ageo[i],abs(AJ[i]).upper()));vgeo[i]=arb(max(vgeo[i],abs(VJ[i]).upper()))
   pb=[arb(0)]+[eval_poly(cp[i],values) for i in range(1,N+1)]
   vb=[arb(0)]+[eval_poly(cv[i],values) for i in range(1,N+1)]
   rb=eval_poly(dr,values);rvb=eval_poly(dv,values)
   for i in range(1,N+1):PB[i]=arb(max(PB[i],pb[i]));VB[i]=arb(max(VB[i],vb[i]))
   RB=arb(max(RB,rb));RVB=arb(max(RVB,rvb));Fmax=arb(max(Fmax,F.upper()))
   cells.append(dict(y_lo=encoded(lo),y_hi=encoded(hi),a_jets=[encoded(x) for x in AJ],v_jets=[encoded(x) for x in VJ],
    P_bounds=[encoded(x) for x in pb],V_bounds=[encoded(x) for x in vb],D_P_N=encoded(rb),D_V_N=encoded(rvb)))
  if k in [0,1,3,7,15,31,55]:print(json.dumps(dict(event='complete_geometry_cover',dyadic_band=k,order=N,D_P_N=str(RB),D_V_N=str(RVB))),flush=True)
 NUMERICS.append(dict(kind='complete_radial_coefficient_bounds',N=N,subdivisions_per_dyadic_band=subdivisions,
  P=[encoded(x) for x in PB],V=[encoded(x) for x in VB],D_P_N=encoded(RB),D_V_N=encoded(RVB),
  a_jets=[encoded(x) for x in ageo],v_jets=[encoded(x) for x in vgeo],F_upper=encoded(Fmax),cells=cells,
  entire_original_thermal_spectrum_in_residual=True,not_a_physical_mode_truncation=True))
 return PB,VB,RB,RVB,ageo,vgeo,Fmax

def execute():
 p,h,zh,env=prior.gy.actual_parameters();p=prior.correlated_parameters(p,h);NUMERICS.extend(prior.NUMERICS)
 P,V=coefficients(11);local_identities(P,V)
 joint=joint_coefficients(p,11)
 for N in [9,11]:
  PB,VB,RB,RVB,AJ,VJ,Fmax=joint_geometry_cover(p,P[:N+1],V[:N+1],N,joint)
  M=((arb(p['q'])+1)/p['C']).sqrt()
  J=sum((PB[j]/M**j for j in range(1,N+1)),arb(0))/2
  JV=sum((VB[j]/M**j for j in range(1,N+1)),arb(0))/2
  R=RB;RV=RVB
  for i in range(1,N+1):
   for j in range(1,N+1):
    if i+j>=N:
     R+=AJ[0]*PB[i]*PB[j]/(2*M**(i+j-N))
     RV+=((VJ[1]+AJ[0]*VJ[0])*PB[i]*PB[j]+AJ[0]*(VB[i]*PB[j]+PB[i]*VB[j]))/(2*M**(i+j-N))
  NUMERICS.append(dict(kind='complete_residual_capacity',N=N,M1=encoded(M),approximation_upper=encoded(J),variation_upper=encoded(JV),
   residual_profile_coefficient=encoded(R),variation_residual_profile_coefficient=encoded(RV),
   residual_profile='(1+(p/(M*fmax))^2)^(-(N+1)/2)/M^N',
   propagated_projection_error_not_yet_computed=True))
  print(json.dumps(dict(event='complete_residual_capacity',N=N,X_upper=str(J),X_v_upper=str(JV),R=str(R),R_v=str(RV))),flush=True)
  import current_signed_resolvent_certificate as certificate
  certificate.certify(sys.modules[__name__],p,P,V,N,PB,VB,R,RV,J,JV,AJ,VJ,Fmax)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes());paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
  {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prior.bh.prior.scalar.prev.endpoint.dependencies()
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,dependencies=deps,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],invocation=sys.argv,
  scope='Complete actual-domain signed interval, subject to every displayed continuum hypothesis and all endpoint/spectral errors.')
 receipt['runs'].append(run)
 def save():out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 save();error=None
 try:
  check('previous_local_tail_run_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',prior.bh.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prior.bh.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 certificates=[r for r in NUMERICS if r['kind']=='whole_signed_projection_certificate']
 certified=status=='passed' and any(r['passes'] for r in certificates)
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  complete_signed_projection_evaluated=bool(certificates) and status=='passed',signed_projection_value_measured=False,
  full_projection_milestone_complete=certified,actual_controlled_solution=False,goal_complete=False)
 save()
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
