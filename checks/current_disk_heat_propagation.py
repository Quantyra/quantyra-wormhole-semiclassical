"""Full analytic transport majorants; actual-domain propagation is separately gated."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json,sys,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,arb_poly,ctx
import sympy as sp
import current_local_boundary_heat_remainder as previous
bh=previous.bh;gy=previous.gy
root=Path(__file__).resolve().parents[1]
STEM='current-disk-heat-propagation';PREVIOUS=previous.STEM
OWN=['notes/'+STEM+'.md','checks/current_disk_heat_propagation.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=bh.encoded
N=144;KEEP=120;CELLS=64
T0=arb('0.001');RA=arb('.75');RT=arb('.70');RC=arb('.625')
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS},name
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name+': '+str(detail))
gy.check=check;bh.prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def aa(x):return arb(abs(x).upper())
def upper(x):return arb(x.upper())
def tr(x):return x.truncate(N+1)
def mul(a,b):return tr(a*b)
def D(a):return tr(a.derivative().left_shift(1))
def df(a):return a.derivative()
def sh(a,k):return tr(a.left_shift(k)) if k>=0 else a.right_shift(-k)
def co(a,z):
 o=arb_poly([])
 for x in reversed(a):o=mul(o,z)+x
 return o
def power(b,c):
 assert b[0]==0
 u=[arb(1)]+[arb(0)]*N
 for n in range(1,N+1):
  u[n]=sum((((n-k)+c*k)*b[k]*u[n-k] for k in range(1,n+1)),arb(0))/n
 return arb_poly(u)
def pa(x,y):return (x[0]+y[0],x[1]+y[1])
def pm(x,y):return (mul(x[0],y[0]),mul(x[0],y[1])+mul(x[1],y[0]))
def ps(x,c):return (x[0]*c,x[1]*c)
def po(x,c):return (power(x[0],c),mul(power(x[0],c+1),x[1])*c)
def op(fn,x):return tuple(fn(a) for a in x)
def sa(x,y):return (x[0]+y[0],x[1]+y[1])
def sm(x,y):return (x[0]*y[0],x[0]*y[1]+x[1]*y[0])
def ss(x,c):return (x[0]*c,x[1]*c)
def pv(a,x):
 out=arb(0)
 for c in reversed(a):out=out*x+c
 return out

def analytic_bounds(F,K,p,zmax):
 d=None
 for ii in range(1,1001):
  trial=arb(ii)/1000
  if p*RA+RA**2*pv(F,trial)/2<trial and RA**2*pv(K,trial)/2<1:
   d=trial;break
 if d is None:raise AssertionError('No analytic Picard ball')
 k=pv(K,d);fz=pv(F,d);qa=[p*p]+[2*F[j-1]/j for j in range(1,7)]
 qz=pv(qa,d);v=(zmax+d)**4;vp=4*(zmax+d)**3
 hh=(RA*zmax**4*p+RA**2*(2*v*fz+vp*(qz+2*p*p))/2)/(1-RA**2*k/2)
 kd=pv([j*K[j] for j in range(1,len(K))],d)*hh+2*v*k+2*vp*fz+12*(zmax+d)**2*qz
 x=k.sqrt()*RA;b=x.sinh()/x-1
 db=RA**2*kd*(1+b)/(2*(1-RA**2*k/2))
 if not b<1:raise AssertionError('Jacobi majorant pole')
 bp=(b,db);bp1=ss(bp,100);Db=ss(bp1,arb('.74'))
 DDb=ss(bp,arb('.74')*100+2*arb('.74')**2*10000)
 A=ss(sa(bp,Db),arb('.5'));DA=ss(sa(Db,DDb),arb('.5'))
 U=(1/(1-b),db/(1-b)**2);U2=sm(U,U);U3=sm(U2,U)
 C=(arb(20000)+1/(arb('.7')*arb('.01')),arb(0))
 C=sa(C,ss(sm(bp1,U),100))
 C=sa(C,ss((U2[0]-1,U2[1]),(arb('.74')*100+2*arb('.74')**2*10000)/arb('.7')**2))
 C=sa(C,ss(sm(Db,U3),arb('.74')*100/arb('.7')**2))
 C=sa(C,ss(sm(A,U2),2*arb('.74')*100/arb('.7')**2))
 C=sa(C,ss(sm(DA,U2),1/arb('.7')**2))
 C=sa(C,ss(sm(sm(A,A),U2),1/arb('.7')**2))
 C=sa(C,ss(sm(sm(Db,A),U3),1/arb('.7')**2))
 C=sa(C,(k/2,kd/2))
 W0=((1-b)**(-arb('.5')),db/2*(1-b)**(-arb('1.5')))
 values=[W0];W=W0
 for j in range(1,4):W=ss(sm(sm(W0,W0),sm(C,W)),arb(1)/j);values.append(W)
 values.append(sm(C,W))
 return dict(radius=d,Kmax=k,kdelta=kd,bmax=b,dbmax=db,norms=values,
  ball_slack=d-p*RA-RA**2*pv(F,d)/2,contraction=RA**2*k/2)

def transport(Q,params,lo,hi):
 zball=(lo+hi)/2+arb(0,(hi-lo)/2)
 F=[aa(pv([arb(i)*Q[i]*comb(i-1,k)/2 for i in range(k+1,7)],zball)) for k in range(6)]
 K=[(k+1)*F[k+1] for k in range(5)]
 P=[params['beta'],2*params['beta'],3*params['beta']+params['u'],arb(0),params['eps']]
 p=(1-lo)*upper(pv(P,zball)).sqrt()
 ab=analytic_bounds(F,K,p,hi)
 Z=[arb(0)]*(N+1);Z[1]=p
 for n in range(N-1):Z[n+2]=co(F,arb_poly(Z).truncate(n+1))[n]/((n+1)*(n+2))
 Z=arb_poly(Z);kz=co(K,Z);fz=co(F,Z)
 qz=co([p*p]+[2*F[j-1]/j for j in range(1,7)],Z)
 base=Z+hi;base2=mul(base,base);vp=mul(base2,base)*4;v=mul(base2,base2)
 src=mul(v,fz)*2+mul(vp,qz+2*p*p)
 H=[arb(0)]*(N+1);H[1]=hi**4*p
 for n in range(N-1):H[n+2]=(sum((kz[j]*H[n-j] for j in range(n+1)),arb(0))+src[n])/((n+1)*(n+2))
 kd=mul(co([j*K[j] for j in range(1,len(K))],Z),arb_poly(H))+mul(v,kz)*2+mul(vp,fz)*2+mul(base2,qz)*12
 B=[arb(0)]*(N+1);B[1]=arb(1);C=[arb(0)]*(N+1)
 for n in range(N-1):
  B[n+2]=sum((kz[j]*B[n-j] for j in range(n+1)),arb(0))/((n+1)*(n+2))
  C[n+2]=sum((kz[j]*C[n-j]+kd[j]*B[n-j] for j in range(n+1)),arb(0))/((n+1)*(n+2))
 b=(arb_poly(B).right_shift(1)-1,arb_poly(C).right_shift(1));kp=(kz,kd)
 U=po(b,arb(1));U2=pm(U,U);U3=pm(U2,U);W0=po(b,arb('.5'));A=ps(pa(b,op(D,b)),arb('.5'))
 def L(w):
  flat=op(lambda a:arb_poly([(k+2)**2*a[k+2] for k in range(N-1)]),w)
  terms=[flat,pm(pm(op(df,b),U),op(df,w)),
   pm(op(lambda a:sh(a,-2),(U2[0]-1,U2[1])),op(lambda a:D(D(a)),w)),
   pm(pm(op(lambda a:sh(D(a),-2),b),U3),op(D,w)),
   pm(pm(op(lambda a:sh(a,-2),ps(A,2)),U2),op(D,w)),
   pm(pm(op(lambda a:sh(D(a),-2),A),U2),w),
   pm(pm(op(lambda a:sh(a,-2),pm(A,A)),U2),w),
   pm(pm(pm(op(lambda a:sh(D(a),-2),b),A),U3),w),ps(pm(kp,w),arb('.5'))]
  out=(arb_poly([]),arb_poly([]))
  for term in terms:out=pa(out,term)
  return out
 W=W0;values=[W]
 for j in range(1,4):
  W=pm(W0,op(lambda a:arb_poly([a[k]/(k+j) for k in range(N+1)]),pm(W0,L(W))))
  values.append(W)
 values.append(L(W))
 return values,ab

TAILSUM={}
def tailsum(p):
 if p not in TAILSUM:
  x=sp.symbols('x');g=x**(KEEP+1)/(1-x)
  for _ in range(p):g=x*sp.diff(g,x)
  a=sp.Rational((g/x**(KEEP+1)).subs(x,sp.Rational(25,28)))
  TAILSUM[p]=arb(int(a.p))/int(a.q)
 return TAILSUM[p]

MOMENTS={}
def moment(n,k,which):
 key=(n,k,which)
 if key not in MOMENTS:
  a=arb(n)/2+k
  val=(4*T0)**(arb(n)/2)
  if which=='row':val*= (a+1).gamma()
  else:val*=1 if a==0 else a**a*(-a).exp()
  MOMENTS[key]=val
 return MOMENTS[key]

def norm(W,M,p,k,which,shift=0):
 val=sum((W[n]*(n**p)*moment(n+shift,k,which) for n in range(KEEP+1)),arb(0))
 val+=M/RT**(KEEP+1)*tailsum(p)*moment(KEEP+1+shift,k,which)
 return val

def kernel_table(values,ab,defect=False):
 tab={key:{which:arb(0) for which in ['inf','row']} for key in ['base','metric','time','mixed']}
 entries=[(3,values[4],ab['norms'][4])] if defect else [(j,values[j],ab['norms'][j]) for j in range(4)]
 for j,W,M in entries:
  for which in ['inf','row']:
   def nm(s=0,p=0,k=0,shift=0):return norm(W[s],M[s],p,k,which,shift)
   def qw(k):return nm(1,0,k)+nm(0,1,k)+nm(0,1,k,1)+2*nm(0,0,k)+nm(0,0,k,1)
   h=abs(j-1);fac=arb(1) if defect else T0**j
   vals={'base':nm(),'time':h*nm()+nm(k=1),'metric':qw(0)+2*nm(k=1),
    'mixed':h*qw(0)+qw(1)+2*(h+1)*nm(k=1)+2*nm(k=2)}
   for key,val in vals.items():
    tab[key][which]+=fac*val*(arb('1.14') if which=='row' else 1/(4*arb.pi()))
 return tab

def encode_table(t):return {k:{w:encoded(x) for w,x in v.items()} for k,v in t.items()}

def decode(x):
 d=x['decimal_enclosure'];return arb(d['mid'],d['rad'])*arb(10)**int(d['exponent'])

def physical_geometry(Q,p,zh):
 qp=arb_poly(Q).derivative();qpp=qp.derivative();q3=qpp.derivative();q4=q3.derivative()
 fact=arb_poly([p['beta'],2*p['beta'],3*p['beta']+p['u'],0,p['eps']])
 maxima={k:arb(0) for k in ['K','gradK','Kyy','v','gradv']}
 for i in range(512):
  z=arb('-.4')+arb('1.4')*(arb(i)+arb('.5'))/512+arb(0,arb('1.4')/1024)
  qb=aa((1-z)**2*fact(z)-p['tau']*z**3);K=-qpp(z)/2
  vals={'K':aa(K),'gradK':aa(q3(z))*qb.sqrt()/2,
   'Kyy':aa(q4(z))*qb/2+aa(q3(z)*qp(z))/4,
   'v':aa(z**4),'gradv':4*aa(z**3)*qb.sqrt()}
  for k,x in vals.items():maxima[k]=arb(max(maxima[k].upper(),x.upper()))
 check('real_geodesic_curvature_guard',maxima['K']<2)
 check('real_geodesic_curvature_gradient_guard',maxima['gradK']<4 and maxima['Kyy']<100)
 # Exact |v|<=1 on [-.4,zh]; interval endpoint overhang is not a larger domain.
 check('real_geodesic_metric_gradient_guard',maxima['gradv']<arb('.5'))
 check('common_two_Schur_volume_guard',(arb(2).sqrt()*RC).sinh()/(arb(2).sqrt()*RC)<arb('1.14'))
 check('fixed_endpoint_Jacobi_angle_guard',arb('.5')/2/(1-2*RC**2/8)<1)
 check('fixed_endpoint_spin_transport_guard',arb('.25')+2*RC**2*arb('.5')/(8*(1-2*RC**2/8))<1)
 check('actual_finite_bolt_parameters',p['gamma']<arb('6e-34') and 1-zh<arb('1e-32'))
 period=2*arb.pi()*p['beta'].sqrt()/p['gamma']
 check('original_thermal_period_guard',period>arb('1e33'))
 check('complete_thermal_image_scale',period**2/(16*T0*arb('1e12'))>arb('1e55'))
 s5=p['gamma']/arb('1.01')*((arb('1.01').sqrt()*5).cosh()-1)
 check('negative_curvature_bolt_bootstrap',s5+1-zh<arb('1e-30'))
 zb=1-arb('5e-31')+arb(0,arb('5e-31'))
 check('actual_bolt_curvature_comparison',-qpp(zb)/2<arb('-.99') and qpp(zb)/2<arb('1.01'))
 check('large_circle_after_bolt_neighborhood',(arb('.99').sqrt()*arb('3.5')).sinh()/arb('.99').sqrt()>16)
 check('outer_monotone_circle_guard',12*p['beta']+4*p['u']-10*aa(p['eps'])>arb('1.5'))
 check('outer_turning_strip_circle_guard',p['beta']-100*p['u']**2>arb('.1'))
 check('conjugate_radius_exceeds_support',arb.pi()/arb(10).sqrt()>RA)
 check('short_loop_cannot_wind',2*RA/16<2*arb.pi())
 NUMERICS.append(dict(kind='actual_geodesic_cover_guards',maxima={k:encoded(x) for k,x in maxima.items()},
  near_bolt_s5=encoded(s5),common_volume_guard='1.14',physical_bolt_regular_domain_retained=True))

def cap_geometry(Q,p):
 # Y=5 lies at .9<z<.995: positive integral comparisons use Q=(1-z)^2 P-tau z^3.
 P=arb_poly([p['beta'],2*p['beta'],3*p['beta']+p['u'],0,p['eps']]);P1=P.derivative();P2=P1.derivative();P3=P2.derivative()
 check('finite_tau_perturbation_guard',0<p['tau']<arb('1e-60'))
 lower0=p['beta']-arb('1e-50');lower1=P(arb('.5'))-arb('1e-50')
 check('positive_cap_start_lower_z',arb(2).log()/lower0.sqrt()+arb(5).log()/lower1.sqrt()<5)
 check('positive_cap_start_upper_z',-arb('.005').log()>5)
 qp=arb_poly(Q).derivative();qpp=qp.derivative();q3=qpp.derivative();q4=q3.derivative()
 posmax=[arb(0)]*4;vmax=[arb(0)]*4;qmax=arb(0)
 for i in range(512):
  lo=arb('.995')*i/512;hi=arb('.995')*(i+1)/512;z=(lo+hi)/2+arb(0,(hi-lo)/2);delta=1-z
  q=delta**2*P(z)-p['tau']*z**3;assert q>0
  qmax=arb(max(qmax.upper(),q.upper()))
  q1=-delta*z*((12*p['beta']+4*p['u'])*z-2*p['u']+p['eps']*(6*z**3-4*z**2))-3*p['tau']*z*z
  k=-qpp(z)/2;ky=-q3(z)*q.sqrt()/2;kyy=-q4(z)*q/2-q3(z)*q1/4
  w=-q1/(2*q.sqrt());w1=k+w*w;w2=ky+2*w*w1;w3=kyy+2*w1*w1+2*w*w2
  fv=[w,w1+w*w,w2+3*w*w1+w**3,w3+4*w*w2+3*w1*w1+6*w*w*w1+w**4]
  vv=[-z**4,-4*z**3*q.sqrt(),-12*z*z*q-2*z**3*q1,q.sqrt()*(-24*z*q-18*z*z*q1-2*z**3*qpp(z))]
  for j in range(4):posmax[j]=arb(max(posmax[j].upper(),aa(fv[j]).upper()));vmax[j]=arb(max(vmax[j].upper(),aa(vv[j]).upper()))
 # Use a rational formula retaining the cancellation in w' at z close to one.
 for i in range(256):
  z=arb('.9')+arb('.095')*(arb(i)+arb('.5'))/256+arb(0,arb('.095')/512);d=1-z
  a=P(z);b=P1(z);c=P2(z);e=P3(z);fuzz=arb(0,arb('1e-40'))
  w=a.sqrt()-d*b/(2*a.sqrt())+fuzz
  w1=d*b-d*d*c/2+d*d*b*b/(4*a)+fuzz
  w2=d*a.sqrt()*(-b+2*d*c-d*d*e/2-d*b*b/(2*a)+d*d*b*c/(2*a)-d*d*b**3/(4*a*a))+fuzz
  q=d*d*a-p['tau']*z**3;q1=qp(z);v2=-12*z*z*q-2*z**3*q1
  assert arb('.8')<w<arb('1.01') and 0<w1<arb('.2') and arb('-.2')<w2<0
  assert 0<v2<arb('.5')
 check('positive_cap_endpoint_signed_guards',True)
 # Actual negative germ coefficients, with explicit Cauchy tails through order four.
 fhalf=p['beta'].sqrt();z=[arb(0)]*(N+3);z[1]=fhalf;F=[(i+1)*Q[i+1]/2 for i in range(6)]
 for n in range(N+1):z[n+2]=co(F,arb_poly(z).truncate(n+1))[n]/((n+1)*(n+2))
 zp=[(n+1)*z[n+1] for n in range(N+1)];ff=[arb(0)]*(N+1);ff[0]=arb(1)
 for n in range(1,N+1):ff[n]=-sum((zp[k]*ff[n-k]/fhalf for k in range(1,n+1)),arb(0))
 vp=-mul(mul(arb_poly(z),arb_poly(z)),mul(arb_poly(z),arb_poly(z)))
 def tail(M,j):
  x=sp.symbols('x');a=sp.Rational(sp.diff(x**(N+1)/(1-x),x,j).subs(x,sp.Rational(5,6)))
  return arb(M)/RA**j*(arb(int(a.p))/int(a.q))
 def jets(poly,M):
  out=[];norms=[]
  for j in range(5):
   out.append(poly(-RC)+arb(0,upper(tail(M,j))))
   norms.append(sum((aa(c)*RC**n for n,c in enumerate(poly.coeffs())),arb(0))+tail(M,j));poly=poly.derivative()
  return out,norms
 fj,fn=jets(arb_poly(ff),'1.72');vj,vn=jets(vp,'.0256')
 fmin=1-(sum((aa(ff[n])*RC**n for n in range(1,N+1)),arb(0))+tail('1.72',0))
 nw=-fj[1]/fj[0];nw1=fj[2]/fj[0]-nw*nw;nw2=-(fj[3]/fj[0]+3*nw*nw1+nw**3)
 # w_out=-w, w_out'=w', w_out''=-w''.
 negative=[aa(nw),aa(nw1),aa(nw2)];negv=[aa(vj[0]),aa(vj[1]),aa(vj[2])]
 allcaps=[]
 for label,L,w0,w1,w2 in [('negative',arb('.5'),*negative),('positive',arb(2),arb('1.01'),arb('.2'),arb('.2'))]:
  a=w0+L*w1+(L*L*w2/2 if label=='negative' else 0)
  ap=w1+L*w2;app=w2;e1=arb('1.875')/L;e2=arb(6)/L**2;e3=arb(60)/L**3
  ww=[a,ap+e1*a,app+2*e1*ap+e2*a,3*e1*app+3*e2*ap+e3*a]
  rel=[ww[0],ww[1]+ww[0]**2,ww[2]+3*ww[0]*ww[1]+ww[0]**3,
       ww[3]+4*ww[0]*ww[2]+3*ww[1]**2+6*ww[0]**2*ww[1]+ww[0]**4]
  for j,guard in enumerate([2,10,100,10000]):check(label+'_cap_relative_derivative_'+str(j+1),rel[j]<guard,rel[j])
  allcaps.append(dict(side=label,relative_derivative_bounds=[encoded(x) for x in rel]))
  if label=='negative':check('negative_cap_positive_value',fj[0]*(-L*a).exp()>arb('.6'))
 for j,guard in enumerate([2,10,100,10000]):check('positive_germ_relative_derivative_'+str(j+1),posmax[j]<guard)
 for j,guard in enumerate([2,10,100,10000]):check('negative_germ_relative_derivative_'+str(j+1),fn[j+1]/fmin<guard)
 for j,guard in enumerate([1,10,100,10000]):check('positive_germ_test_derivative_'+str(j),vmax[j]<guard)
 for j,guard in enumerate([1,10,100,10000]):check('negative_germ_test_derivative_'+str(j),vn[j]<guard)
 L=arb('.5');qa=negv[1]+L*negv[2]
 check('negative_cap_test_guards',negv[0]+L*negv[1]/2+L*L*negv[2]/7<1 and qa<10 and negv[2]+arb('1.875')/L*qa<100 and 2*arb('1.875')/L*negv[2]+6/L**2*qa<10000)
 check('positive_cap_test_guards',-arb('.9')**4+arb(4)/7*arb('.5')<0 and arb(4)*arb('.1')+2*arb('.5')<10)
 check('positive_cap_test_second_third_guards',arb('.5')+arb('1.875')/2*arb('1.4')<100 and 2*arb('1.875')/2*arb('.5')+arb(6)/4*arb('1.4')<10000)
 check('positive_cap_test_range_identity',2*p['beta']-4*aa(p['eps'])>0)
 check('positive_germ_warping_lower',p['beta']/qmax>arb('.6')**2)
 check('negative_germ_warping_values',fmin>arb('.6') and fn[0]<arb('1e6'))
 check('global_cap_warping_upper',p['beta'].sqrt()/(arb('.005')*(P(arb('.9'))-arb('1e-50')).sqrt())*(2*arb('1.41')).exp()<arb('1e6'))
 check('extended_cylinder_circles_large',p['beta'].sqrt()/(p['gamma']*arb('1e6'))>arb('1e20'))
 NUMERICS.append(dict(kind='enlarged_actual_germ_cap_guards',caps=allcaps,positive_germ_rel=[encoded(x) for x in posmax],
  negative_endpoint_f=[encoded(x) for x in fj],negative_endpoint_v=[encoded(x) for x in vj],
  negative_endpoint_w=[encoded(x) for x in negative],global_f_range=['.6','1e6'],
  global_relative_f_derivatives=['2','10','100','10000'],global_v_derivatives=['1','10','100','10000']))

UPPERGAMMA={}
def upper_gamma_half(n,x):
 # Gamma((n+1)/2,x), with integer n>=0, by the positive recurrence.
 if n%2==0:g=arb.pi().sqrt()*x.sqrt().erfc();a=arb('.5');start=0
 else:g=(-x).exp();a=arb(1);start=1
 for j in range(start,n,2):g=a*g+x**a*(-x).exp();a+=1
 return g

def cutmoment(n,which):
 key=(n,which)
 if key not in UPPERGAMMA:
  a=arb('.5')**2/(4*T0);power=arb(n)/2
  if which=='row':val=(4*T0)**power*upper_gamma_half(n+1,a)
  else:
   where=arb(max(a.upper(),power.upper()))
   val=(4*T0)**power*where**power*(-where).exp()
  UPPERGAMMA[key]=val
 return UPPERGAMMA[key]

def cutoff_norm(values,ab):
 check('cut_geometry_variation_guard_'+str(len(NUMERICS)),ab['bmax']<arb('.2') and ab['dbmax']<100 and ab['kdelta']<100)
 out={w:arb(0) for w in ['inf','row']}
 for j in range(4):
  W=values[j];M=ab['norms'][j]
  for which in out:
   val=arb(0)
   for s,powr in [(0,0),(0,1),(0,2),(1,0),(1,1)]:
    val+=sum((W[s][n]*n**powr*cutmoment(n,which) for n in range(KEEP+1)),arb(0))
    val+=M[s]/RT**(KEEP+1)*tailsum(powr)*cutmoment(KEEP+1,which)
   out[which]+=T0**j*val
 for which in out:out[which]*=arb('1e8')*T0**-6*(arb('1.14') if which=='row' else 1/(4*arb.pi()))
 return out

def commutator_polynomial_check():
 # Positive square-free jet algebra in metric,time,endpoint-coordinate.
 a,l,u,w,xi=sp.symbols('a l u w xi',nonnegative=True);d=a+l
 def jp(v):return [sp.expand(x) for x in v]
 def jm(A,B):
  out=[sp.Integer(0)]*8
  for i in range(8):
   for j in range(8):
    if not i&j:out[i|j]+=A[i]*B[j]
  return jp(out)
 def ja(A,B):return jp([x+y for x,y in zip(A,B)])
 def js(A,c):return jp([x*c for x in A])
 def fun(X,cs):
  out=[sp.Integer(0)]*8;power=[sp.Integer(1)]+[sp.Integer(0)]*7
  for c in cs:out=ja(out,js(power,c));power=jm(power,X)
  return out
 f=jp([1,1+2*d,l,l*(7+5*d),2,16+10*d,5*l,l*(100+50*d)])
 g=jp([2,14+10*d,5*l,l*(95+50*d),10,200+100*d,50*l,l*(5950+5000*d)])
 F=jm(f,f);F[0]=0
 inv=fun(F,[1,1,1,1]);log=fun(F,[0,1,sp.Rational(1,2),sp.Rational(1,3)])
 tjet=[1,0,1,0,0,0,0,0] # t<=1, D t=t.
 invt=[1,0,1,0,0,0,0,0] # Absolute derivatives of 1/t.
 logjet=js(log,sp.Rational(1,2))
 # Endpoint prefactor, -log(t), and radial phase absolute derivatives.
 endpoint=jp([0,1,1,0,1,5,0,0])
 radial=jp([0,2*u,u,2*u,a*xi/2,a*xi,a*xi/2,a*xi])
 tangent=js(jm(invt,inv),w);tangent[0]=0
 potential=js(jm(tjet,jm(jm(g,g),inv)),sp.Rational(1,4));potential[0]=0
 phase=js(jm(g,inv),1+w);phase[0]=0 # sqrt(t*uT)<=1+uT.
 for piece in [endpoint,radial,tangent,potential,phase]:logjet=ja(logjet,piece)
 kernel=fun(logjet,[1,1,sp.Rational(1,2),sp.Rational(1,6)])
 # [L,chi] has two derivatives of chi. Its fixed-coordinate conformal
 # metric variation has bounded coefficients; 100 dominates all of them.
 alljets=sp.Poly(sp.expand(100*sum(kernel)),a,l,u,w,xi)
 worst=sp.Rational(0);degree=0
 for powers,c in alljets.terms():
  n=sum(powers[:4]);assert powers[4]<=3 and n<=12 and c>=0
  multinomial=sp.factorial(12)/sp.factorial(12-n)
  for power in powers[:4]:multinomial/=sp.factorial(power)
  worst=max(worst,c/multinomial);degree=max(degree,n)
 check('full_differentiated_commutator_polynomial',worst<10**12,str(worst))
 NUMERICS.append(dict(kind='commutator_complete_jet_majorant',jet_order=['metric','time','endpoint'],
  normalized_f_jets=[str(x) for x in f],normalized_G_jets=[str(x) for x in g],
  polynomial=str(alljets.as_expr()),maximum_coefficient_ratio=str(worst),degree=degree,
  full_majorant='1e12 t^-3 (1+a+ell+uY+uT)^12'))

def commutator_norm():
 x,y,s,rt=sp.symbols('x y s rt',nonnegative=True)
 poly=sp.Poly(sp.expand((1+2*rt*x+rt*s+x*x+y*y)**12),x,y,s,rt)
 eta=arb(1)/20;low=1/(2*T0.sqrt());out={w:arb(0) for w in ['inf','row']};thermal_prefactor=arb(0)
 def gm(n,c=arb(1)):
  if n==0:return arb(1)
  a=arb(n)/2;return (a/c)**a*(-a).exp()
 for (nx,ny,ns,nt),coef in poly.terms():
  assert (1-eta)*low**2>arb(nx)/2
  cf=arb(int(coef))*T0.sqrt()**nt
  smom=arb(1) if ns==0 else 2*(arb(ns)/2+1).gamma()
  tilted=1+2*eta/(1-eta) if ns==0 else smom/(1-eta)**(arb(ns)/2+1)
  out['inf']+=cf*low**nx*(-(1-eta)*low**2).exp()*gm(ny)*tilted
  out['row']+=cf*upper_gamma_half(nx,(1-eta)*low**2)*(arb(ny+1)/2).gamma()/arb.pi()/(1-eta)**(arb(nx+1)/2)*smom
  extra=2*(arb(ns+1)/2+1).gamma()/(1-eta)**(arb(ns+1)/2+1)
  thermal_prefactor+=cf*gm(nx,1-eta)*gm(ny,arb('.5'))*(tilted+extra)
 out['inf']*=2*arb('1e12')*T0**-6*(T0+2*T0/eta).exp()/(4*arb.pi())
 out['row']*=arb('1e12')*T0**-6*(T0+T0/eta).exp()
 check('commutator_full_time_monotonicity',(1-eta)/(4*T0)>18)
 check('all_commutator_tails',all(x<arb('1e-25') for x in out.values()),str(out))
 thermal_prefactor*=arb('1e14')*(T0+2*T0/eta).exp()/(4*arb.pi())
 check('explicit_complete_thermal_prefactor',thermal_prefactor<arb('1e100'))
 # Entire thermal image sum, using the original period and f<=1e6.
 thermal=arb('1e100')*T0**-6*2**101*arb(101).gamma()*arb('1e55')**-100
 check('all_thermal_image_jets',thermal<arb('1e-100'))
 NUMERICS.append(dict(kind='complete_commutator_and_thermal_tails',commutator={k:encoded(v) for k,v in out.items()},
  full_thermal_image_guard=encoded(thermal),spatial_gap='1',global_cap_f_upper='1e6'))
 NUMERICS[-1]['complete_thermal_prefactor']=encoded(thermal_prefactor)
 return out,thermal

def cap_table():
 x,y,s=sp.symbols('x y s',nonnegative=True);rt=sp.symbols('rt',positive=True)
 ell=rt*s;d=rt*(2*x+s);t=rt**2
 A=2+2*d+t*(16+14*d)+2*x*x+(2+4*d)*y*y+18*rt*y*(1+d)
 B=1+ell+t*(1+7*ell)+x*x+(1+2*ell)*y*y+9*ell*rt*y
 C=ell*(10+11*d)+t*(16+14*d+ell*(192+145*d))+2*x*x+y*y*(2+4*d+ell*(24+30*d))+ell*(181+150*d)*rt*y
 polys={'base':sp.Integer(1),'metric':A,'time':B,'mixed':A*B+C}
 tab={k:{w:arb(0) for w in ['inf','row']} for k in polys}
 eta=arb(1)/20
 def gaussmax(n,c):
  if n==0:return arb(1)
  a=arb(n)/2;return (a/c)**a*(-a).exp()
 for key,expr in polys.items():
  for (nx,ny,ns,nt),coef in sp.Poly(sp.expand(expr),x,y,s,rt).terms():
   assert coef>0
   cf=arb(int(coef.p))/int(coef.q)*T0.sqrt()**nt
   smom=arb(1) if ns==0 else 2*(arb(ns)/2+1).gamma()
   tilted=1+2*eta/(1-eta) if ns==0 else smom/(1-eta)**(arb(ns)/2+1)
   tab[key]['inf']+=cf*gaussmax(nx,1-eta)*gaussmax(ny,arb(1))*tilted
   tab[key]['row']+=cf*(arb(nx+1)/2).gamma()*(arb(ny+1)/2).gamma()/arb.pi()/(1-eta)**(arb(nx+1)/2)*smom
  tab[key]['inf']*=(T0+2*T0/eta).exp()/(4*arb.pi())
  tab[key]['row']*=(T0+T0/eta).exp()
 NUMERICS.append(dict(kind='conditional_complete_cap_kernel_table',table=encode_table(tab),
  polynomials={k:str(sp.expand(v)) for k,v in polys.items()},
  enlarged_cap_guards_verified=True,full_periodization_verified=True))
 return tab

def convolve_jets(P,F):
 def cv(p,f):return p['inf']*f['row']/32+7*p['row']*f['inf']/24
 return {'base':cv(P['base'],F['base']),
  'metric':cv(P['metric'],F['base'])+cv(P['base'],F['metric']),
  'time':cv(P['base'],F['base'])+cv(P['time'],F['base'])+cv(P['base'],F['time']),
  'mixed':cv(P['metric'],F['base'])+cv(P['base'],F['metric'])+cv(P['mixed'],F['base'])+
          cv(P['time'],F['metric'])+cv(P['metric'],F['time'])+cv(P['base'],F['mixed'])}

def complete_transfer(worst,cut):
 cap=cap_table();P={}
 comm,thermal=commutator_norm()
 for k in cap:
  cap[k]['inf']+=thermal
  for w in ['inf','row']:
   worst['F'][k][w]+=cut[w]+comm[w]+thermal
   worst['P'][k][w]+=cut[w]*T0**3
 for k in cap:
  P[k]={'inf':arb(max((2*cap[k]['inf']).upper(),worst['P'][k]['inf'].upper())),
        'row':cap[k]['row']+worst['P'][k]['row']}
 error=convolve_jets(P,worst['F']);smooth=convolve_jets(cap,worst['F'])
 rows=worst['F'];bbar=2*rows['base']['row']+2*rows['metric']['row']+rows['time']['row']+rows['mixed']['row']
 x=T0**4*bbar/140
 check('conditional_Volterra_ratio_small',x<arb('1e-9'))
 higher=sum(error.values(),arb(0))*x/(1-x)
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 oldbound=next(x for x in old['numerics'] if x['kind']=='whole_local_path_bound')
 boundary_compare=2*decode(oldbound['all_auxiliary_excursions_D_bound'])+decode(oldbound['all_nonzero_winding_D_bound'])+thermal
 boundary_new=decode(oldbound['zero_winding_D_bound'])+decode(oldbound['all_auxiliary_excursions_D_bound'])+thermal
 check('enlarged_boundary_component_still_within_149',boundary_new<149)
 check('original_boundary_split_comparison',boundary_compare<arb('4.041e-11'))
 bound=4*(error['mixed']+higher+smooth['mixed'])+1/(600*arb.pi())+boundary_compare
 check('complete_actual_complement_target',bound<952 and bound<arb('1628.7385'),bound)
 Dguard=1101;Bguard=arb('1.2e-8')*Dguard
 heatcap=arb('2.5e-43')/arb('1.0666666e-11')*(arb('2.3454197e-13')*T0+Bguard*T0*T0/2)
 check('complete_heat_cap_target',heatcap<arb('1.5483364e-43') and heatcap<arb('2.5e-43'))
 NUMERICS.append(dict(kind='complete_actual_domain_propagation_budget',first_actual={k:encoded(x) for k,x in error.items()},
  smooth_cap_comparison={k:encoded(x) for k,x in smooth.items()},higher_Volterra=encoded(higher),
  old_boundary_split_comparison=encoded(boundary_compare),new_boundary_D=encoded(boundary_new),
  actual_D_with_cutoffs=encoded(bound),required_D='1628.7385',cutoff_jet_norm={k:encoded(v) for k,v in cut.items()},
  complement_D_guard=952,total_D_guard=Dguard,B_heat_guard=encoded(Bguard),full_heat_cap_upper=encoded(heatcap),
  actual_domain_remainder_certified=True,full_heat_cap_certified=True,complete_signed_projection_evaluated=False))
 print(json.dumps(dict(actual_complement_D=str(bound),heat_cap=str(heatcap))),flush=True)

def main_calc():
 p,h,zh,env=gy.actual_parameters()
 Q=[p['beta'],arb(0),p['u'],-4*p['beta']-2*p['u']-p['tau'],3*p['beta']+p['u']+p['eps'],-2*p['eps'],p['eps']]
 check('positive_actual_tau',p['tau']>0)
 check('full_series_retention_margin',KEEP+12<N)
 physical_geometry(Q,p,zh);cap_geometry(Q,p);commutator_polynomial_check()
 worst={kind:{k:{w:arb(0) for w in ['inf','row']} for k in ['base','metric','time','mixed']} for kind in ['P','F']}
 cuts={w:arb(0) for w in ['inf','row']}
 for i in range(CELLS):
  lo=arb(i)/CELLS;hi=arb(i+1)/CELLS
  values,ab=transport(Q,p,lo,hi)
  check('analytic_ball_'+str(i),ab['ball_slack']>0 and ab['contraction']<1 and ab['bmax']<1)
  tabs={'P':kernel_table(values,ab),'F':kernel_table(values,ab,True)}
  cut=cutoff_norm(values,ab)
  for w,v in cut.items():cuts[w]=arb(max(cuts[w].upper(),v.upper()))
  for kind,tab in tabs.items():
   for k,v in tab.items():
    for w,x in v.items():worst[kind][k][w]=arb(max(worst[kind][k][w].upper(),x.upper()))
  NUMERICS.append(dict(kind='complete_transport_cell',index=i,z_interval=[str(lo),str(hi)],
   analytic_ball=encoded(ab['radius']),analytic_contraction=encoded(ab['contraction']),
   Kmax=encoded(ab['Kmax']),kdelta=encoded(ab['kdelta']),bmax=encoded(ab['bmax']),dbmax=encoded(ab['dbmax']),
   analytic_norms=[[encoded(x) for x in pair] for pair in ab['norms']],
   low_coefficients=[[[encoded(a[n]) for n in range(KEEP+1)] for a in pair] for pair in values],
   kernel_tables={k:encode_table(t) for k,t in tabs.items()}))
  if i%8==7:print(json.dumps(dict(progress_cells=i+1,total=CELLS)),flush=True)
 product=1/(600*arb.pi())
 check('complete_product_contact_remainder',product<arb('.000531'))
 NUMERICS.append(dict(kind='full_uncut_transport_norms',tables={k:encode_table(t) for k,t in worst.items()},
  product_contact_D=encoded(product),actual_disk_propagation_certified=False,
  explanation='Full analytic tails included. Cutoffs, boundary-compatible kernel and propagation not yet included.'))
 complete_transfer(worst,cuts)
 print(json.dumps({kind:{k:{w:str(x) for w,x in v.items()} for k,v in tab.items()} for kind,tab in worst.items()}),flush=True)

def main():
 global T0,RA,RT,RC
 ctx.prec=512
 T0=arb(1)/1000;RA=arb(3)/4;RT=arb(7)/10;RC=arb(5)/8
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=bh.prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  invocation=sys.argv,dependencies=deps,frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 def save():out.write_text(json.dumps(receipt,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
 save()
 try:
  check('previous_scoped_execution_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',bh.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  main_calc()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depsafter=bh.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depsafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depsafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  continuum_argument_scope='Ordinary paired-domain construction, regular-bolt chart cover, full analytic transport defect, differentiated cap/cutoff tails and Volterra propagation. Checks verify explicit inequalities.',
  complete_actual_domain_remainder_bounded=status=='passed',full_heat_cap_certified=status=='passed',complete_signed_projection_evaluated=False,
  controlled_solution_certified=False,goal_complete=False)
 save()
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
