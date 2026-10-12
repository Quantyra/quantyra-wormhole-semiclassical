"""Stationary electric/axial RN fields and a prescribed reflection-symmetric junction.
This does not construct a smooth semiclassical binary or its quantum state.
Run with numpy, scipy and sympy installed. No Lean or cloud command is used.
"""
from pathlib import Path
import hashlib,json,platform,os,gc
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy
from scipy.integrate import solve_ivp,solve_bvp,quad
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def tensors(kind):
 r,c=s.symbols('r c',real=True)
 V,W=[s.Function(n)(r) for n in ['V','W']]
 f=(r-1)**2/r**2 if kind=='RN' else 1+r*r
 b=r if kind=='RN' else s.S.One
 g=s.diag(-f,1/f,b*b/(1-c*c),b*b*(1-c*c)); inv=g.inv()
 h=s.zeros(4); h[0,3]=h[3,0]=-b*b*W*(1-c*c)
 di=-inv*h*inv
 def d(x,k): return s.diff(x,r) if k==1 else s.diff(x,c) if k==2 else s.S.Zero
 G,D={},{}
 for i in range(4):
  for j in range(4):
   for k in range(4):
    G[i,j,k]=s.factor(sum(inv[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))/2 for l in range(4)))
    D[i,j,k]=s.factor(sum((di[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))+inv[i,l]*(d(h[l,k],j)+d(h[l,j],k)-d(h[j,k],l)))/2 for l in range(4)))
 checks={}
 Ric,dRic=s.zeros(4),s.zeros(4)
 for i in range(4):
  for j in range(i,4):
   base=delta=s.S.Zero
   for k in range(4):
    base+=d(G[k,i,j],k)-d(G[k,i,k],j)
    delta+=d(D[k,i,j],k)-d(D[k,i,k],j)
    for l in range(4):
     base+=G[k,i,j]*G[l,k,l]-G[l,i,k]*G[k,j,l]
     delta+=D[k,i,j]*G[l,k,l]+G[k,i,j]*D[l,k,l]-D[l,i,k]*G[k,j,l]-G[l,i,k]*D[k,j,l]
   Ric[i,j]=Ric[j,i]=s.factor(base); dRic[i,j]=dRic[j,i]=s.factor(delta)
 scalar=s.trace(inv*Ric); ds=s.factor(s.trace(di*Ric+inv*dRic))
 F,dF=s.zeros(4),s.zeros(4); F[2,3]=-1;F[3,2]=1
 dF[1,0]=s.diff(V,r)*c;dF[0,1]=-dF[1,0];dF[2,0]=V;dF[0,2]=-V
 FU=inv*F*inv; dFU=di*F*inv+inv*dF*inv+inv*F*di
 square=sum(F[i,j]*FU[i,j] for i in range(4) for j in range(4))
 dsquare=sum(dF[i,j]*FU[i,j]+F[i,j]*dFU[i,j] for i in range(4) for j in range(4))
 stress=2*(F*inv*F.T-g*square/4)
 dstress=2*(dF*inv*F.T+F*di*F.T+F*inv*dF.T-h*square/4-g*dsquare/4)
 E=dRic-h*scalar/2-g*ds/2-dstress
 for i in range(4):
  for j in range(i,4):
   assert s.factor(Ric[i,j]-g[i,j]*scalar/2-stress[i,j])==0
   checks[f'background_E{i}{j}']='exact_zero'
 vol=s.trace(inv*h)/2
 raw={f'E{i}{j}':s.factor(E[i,j]) for i in range(4) for j in range(i,4)}
 for j in range(4):
  assert s.factor(sum(d(b*b*FU[i,j],i) for i in range(4)))==0
  checks[f'background_M{j}']='exact_zero'
  raw[f'M{j}']=s.factor(sum(d(b*b*(vol*FU[i,j]+dFU[i,j]),i) for i in range(4)))
 phi=V-W
 maxwell=s.diff(b*b*s.diff(V,r),r)-2*phi/f
 einstein=s.diff(b**4*s.diff(W,r),r)+4*phi/f
 assert s.factor(raw['M0']+c*maxwell)==0
 assert s.factor(raw['E03']-(1-c*c)*f*einstein/(2*b*b))==0
 checks['off_shell_equation_reduction']='exact_zero'
 e,J,k,L,beta=s.symbols('e J k L beta',real=True)
 H=r*(1+r*r)*s.atan(r)+r*r+s.Rational(2,3)
 if kind=='RN':
  vs=e*(r-2+3/r)+(J-e)/r**2+beta
  ws=2*e/r+(J-e)*(2*r-1)/r**4+beta
 else:
  vs=k*H+L*r+beta;ws=-2*k*H+L*r+beta
 sub={s.diff(field,r,n):s.diff(val,r,n) for field,val in [(V,vs),(W,ws)] for n in range(3)}
 for key,expr in raw.items():
  assert s.simplify(expr.xreplace(sub))==0,(kind,key)
  checks['regular_'+key]='exact_zero'
 if kind=='RN':
  u1=(2*r-1)*(5*r*r-5*r+2)/(30*r*(r-1)**3)
  z=s.symbols('frame_value',real=True)
  sub={V:z+f*s.diff(u1,r)/2,W:z,s.diff(V,r):u1/r**2,
       s.diff(V,r,2):s.diff(u1/r**2,r),s.diff(W,r):-2*u1/r**4,
       s.diff(W,r,2):s.diff(-2*u1/r**4,r)}
  for key,expr in raw.items():
   assert s.factor(expr.xreplace(sub))==0,key
   checks['decaying_'+key]='exact_zero'
 print(kind,'full_tensor_checks',len(checks),flush=True)
 return checks

def basis_checks():
 r,z,e,J=s.symbols('r z e J',positive=True)
 f=(r-1)**2/r**2
 u0=(r-1)**2*(r+2)/r
 u1=(2*r-1)*(5*r*r-5*r+2)/(30*r*(r-1)**3)
 op=lambda u:s.diff(f*s.diff(u,r),r)-(2/r**2+4/r**4)*u
 H=r*(1+r*r)*s.atan(r)+r*r+s.Rational(2,3)
 g=s.diff(H,r)
 identities={
  'growing_solution':op(u0),'decaying_solution':op(u1),
  'spin_particular':op(-2*J/r)-12*J/r**4,
  'Wronskian':f*(u0*s.diff(u1,r)-s.diff(u0,r)*u1)+1,
  'numerical_basis_derivative':s.diff(u1,r)-u_basis(r)[3],
  'throat_potential':(1+r*r)*s.diff(H,r,2)-6*H,
  'throat_flux':s.diff((1+r*r)*s.diff(g,r),r)-6*g,
  'RN_near_regular':s.limit(u0/(r-1)**2,r,1)-3,
  'RN_near_decaying':s.limit(u1*(r-1)**3,r,1)-s.Rational(1,15),
  'throat_growing':s.limit(g/r**2,r,s.oo)-3*s.pi/2}
 for key,expr in identities.items():assert s.simplify(expr)==0,key
 # u1>0 and u1'<0 for r=1+z, z>0, directly from positive polynomials.
 minus_derivative=s.factor(-s.diff(u1,r)).subs(r,1+z)
 numerator,denominator=s.fraction(s.factor(minus_derivative))
 assert all(v>0 for v in s.Poly(s.expand(numerator),z).all_coeffs())
 assert all(v>=0 for v in s.Poly(s.expand(denominator),z).all_coeffs())
 # General-mass/charge reference branches and the neutral Maxwell limit.
 M,P=s.symbols('M P',positive=True)
 fgen=1-2*M/r+P*P/r**2
 V=e*(r-2*M+3*P*P/r)+P*(J-e*P**3)/(M*r*r)
 W=2*P*e/r+(J-e*P**3)*(2*M*r-P*P)/(M*r**4)
 assert s.factor(s.diff(r*r*s.diff(V,r),r)-2*(V-P*W)/fgen)==0
 assert s.factor(r**4*s.diff(W,r)+2*P*r*r*s.diff(V,r)+6*J)==0
 G,q,S=s.symbols('G q S',positive=True)
 lam=G*q/(16*P*P)
 spin=G*S/(P*P)
 energy=(s.sqrt((P*P+s.sqrt(P**4+4*G*G*S*S))/2)-P)/G
 identities['spin_granularity']=q*spin/lam-16*S
 identities['LLL_spin_norm']=(q-1)*spin/(2*lam)-8*(1-1/q)*S
 identities['KN_rotation_energy']=s.series(energy,S,0,4).removeO()-G*S*S/(2*P**3)
 identities['scaled_rotation_energy']=P/lam*G*S*S/(2*P**3)-8*S*S/q
 for key,expr in identities.items():assert s.simplify(expr)==0,key
 return dict(count=len(identities)+3,identities={k:'exact_zero' for k in identities},
             decreasing_u1_numerator=str(s.expand(numerator)),
             neutral_photon_potential='2*(1-2M/r)/r^2; not the l=1 Regge-Wheeler potential')

def u_basis(r):
 z=r-1
 u0=z*z*(r+2)/r
 u0p=2*r-2/r**2
 u1=(2*r-1)*(5*r*r-5*r+2)/(30*r*z**3)
 u1p=-(5*z**4+15*z**3+21*z*z+13*z+3)/(15*r*r*z**4)
 return u0,u0p,u1,u1p

def throat_basis(r):
 H=r*(1+r*r)*np.arctan(r)+r*r+2/3
 g=(1+3*r*r)*np.arctan(r)+3*r
 gp=6*H/(1+r*r)
 return H,g,gp

def matched_checks():
 rows=[]
 for lam in [.01,.003,.001,.0003,.0001]:
  for scale in [.7,1.,1.4]:
   R=scale*lam**(-1/3);lo=lam*R;r0=1+lo;rout=21.
   H,g,gp=throat_basis(R);u0,u0p,u1,u1p=u_basis(r0);f=lo*lo/r0**2
   denominator=lam*(1+R*R)*gp*u1-g*f*u1p
   k=1/denominator;B=(k*g-u0)/u1
   assert denominator>0
   assert abs(lam*(1+R*R)*k*gp-f*(u0p+B*u1p))<2e-12
   def values(r):
    a,ap,b,bp=u_basis(r)
    return a+B*b,ap+B*bp
   points=np.geomspace(lo,rout-1,100)
   def rhs(t,y):
    z=np.exp(t);r=1+z
    return [y[1],(2+4/r**2)*y[0]-(2/(1+z)-1)*y[1]]
   outputs=[]
   for method,tol in [('DOP853',2e-12),('Radau',2e-11)]:
    initial=[k*g,lo*(u0p+B*u1p)]
    sol=solve_ivp(rhs,(np.log(lo),np.log(rout-1)),initial,t_eval=np.log(points),method=method,rtol=tol,atol=tol/100)
    assert sol.success,sol.message
    outputs.append(sol.y[0])
   exact=np.array([values(1+z)[0] for z in points])
   relative=max(np.max(np.abs(y-exact)/(1+np.abs(exact))) for y in outputs)
   assert relative<1e-8,relative
   integral=quad(lambda r:values(r)[0]/r**4,r0,rout,epsabs=2e-12,epsrel=2e-12)[0]
   beta=2*integral/lam+2*k*H
   regular_outer_W=2/rout-(2*rout-1)/rout**4
   beta_reference=(1-regular_outer_W)/lam
   geometric_error=abs(k/(2*lam*lam/np.pi)-1)
   # Odd spin sector, J_plus=1, J_minus=-1, no additional growing electric field.
   hp=1+3*R*R
   matrix=np.array([[u1,-hp],[f*u1p,-6*lam*R*(1+R*R)]])
   bs,ks=np.linalg.solve(matrix,np.array([2*(1/r0-1),-2*f/r0**2]))
   def spin_rhs(t,y):
    result=rhs(t,y);result[1]+=12/(1+np.exp(t))**2;return result
   exact_spin=np.array([bs*u_basis(1+z)[2]-2/(1+z) for z in points])
   outer_u1,outer_u1p=u_basis(rout)[2:]
   def spin_boundary(left,right):
    return np.array([left[1]-lo/f*6*lam*R*(1+R*R)/hp*(left[0]+2),
      right[1]-(rout-1)*((right[0]+2/rout)*outer_u1p/outer_u1+2/rout**2)])
   boundary_residual=spin_boundary(np.array([bs*u1-2/r0,lo*(bs*u1p+2/r0**2)]),
      np.array([bs*outer_u1-2/rout,(rout-1)*(bs*outer_u1p+2/rout**2)]))
   assert np.max(np.abs(boundary_residual))<2e-11
   # A boundary solver avoids amplifying a roundoff-sized growing mode
   # while shooting the small decaying spin solution outward.
   mesh=np.log(points);guess=np.array([-2/(1+points),2*points/(1+points)**2])
   spin_errors=[]
   for tol in [1e-8,2e-10]:
    sol=solve_bvp(spin_rhs,spin_boundary,mesh,guess,tol=tol,max_nodes=5000)
    assert sol.success,sol.message
    spin_errors.append(float(np.max(np.abs(sol.sol(mesh)[0]-exact_spin)/(1+np.abs(exact_spin)))))
   spin_error=spin_errors[-1]
   assert spin_error<1e-8 and spin_error<spin_errors[0]/5,spin_errors
   core_W=-2*ks*R*(1+R*R)-2*R
   outer_W=lam*core_W+quad(lambda r:(-6-2*(bs*u_basis(r)[2]-2/r))/r**4,r0,rout,epsabs=2e-12,epsrel=2e-12)[0]
   resistance=-outer_W
   assert resistance>0 and ks>0
   rows.append(dict(lam=lam,R=R,z_inner=lo,r_outer=rout,k=k,decaying_coefficient=B,
                    source_amplitude=1.,minimum_denominator=denominator,
                    growing_coefficient_relative_error=geometric_error,
                    growing_error_over_geometric_scale=geometric_error/(lo+R**-2),
                    beta=beta,beta_reference=beta_reference,
                    beta_relative_error=abs(beta/beta_reference-1),
                    spin_core_coefficient=ks,spin_decaying_coefficient=bs,
                    spin_rotation_resistance=resistance,
                    spin_resistance_reference=1-(2*rout-1)/rout**4,
                    spin_ODE_relative_error=spin_error,
                    spin_BVP_refinement_errors=spin_errors,
                    maximum_ODE_relative_error=float(relative),
                    core_Phi_center=2*k,core_W_center=beta-4*k/3))
 return rows

def scales():
 rows=[]
 for q in [1e3,1e6,1e9,1e12]:
  lam=1/(16*np.pi*q**3);delta=np.log(q)**-2;chi=lam/delta
  v=np.sqrt(chi);sigma=2*v/delta;e=2*v*chi*chi
  beta=e/lam
  assert np.isclose(beta,chi*sigma,rtol=1e-14,atol=0)
  rows.append(dict(q=q,electric_drive=e,chi=chi,Omega_ell=sigma,
                   zero_rotation_charge_frame_beta=beta,
                   LLL_constant_rotation_bound=(q-1)*abs(beta)/2,
                   LLL_unit_physical_spin_norm=8*(1-1/q),
                   classical_unit_spin_energy_times_ell=8/q,
                   spin_parameter_unit_rotation_threshold=lam/q,
                   asymptotic_core_flux_coefficient=2*e*lam*lam/np.pi))
 return rows

def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
 tensor_checks={kind:tensors(kind) for kind in ['RN','throat']}
 basis=basis_checks()
 s.core.cache.clear_cache()
 gc.collect()
 result=dict(status='passed_stationary_electric_axial_components_prescribed_junction',
  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=s.__version__,
  tensors=tensor_checks,basis=basis,
  matched_cases=matched_checks(),scale_cases=scales(),script_sha256=digest(Path(__file__)),
  input_sha256={n:digest(ROOT/n) for n in ['notes/stationary-electric-dipole.md','notes/binary-mouth-drive.md',
    'notes/weak-coupling-window.md','notes/lll-return-holonomy.md']},
  complete_four_dimensional_junction=False,quantum_state_or_response_derived=False,
  physical_spin_sector_selected=False,novelty_established=False,lean_invocations=0,goal_complete=False)
 (ROOT/'evidence/stationary-electric-dipole.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(dict(status=result['status'],tensor_checks=sum(len(x) for x in result['tensors'].values()),
  basis_checks=result['basis']['count'],matched_cases=len(result['matched_cases']),
  max_ODE_error=max(x['maximum_ODE_relative_error'] for x in result['matched_cases']))))

if __name__=='__main__':main()
