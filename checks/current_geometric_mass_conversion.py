"""Finite-mode geometric conversion and the quantitative cancellation test."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'checks'))
import current_product_matched_source_projection as prior
from flint import arb,ctx
encoded=prior.encoded
STEM='current-geometric-mass-conversion';PREVIOUS='current-product-matched-source-projection'
OWN=['notes/'+STEM+'.md','checks/current_geometric_mass_conversion.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[]
def check(name,ok,detail=''):
 assert name not in {x['name'] for x in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
prior.scalar.prev.york.check=check;prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def algebra():
 t,b,R,db,M,mu=sp.symbols('t b R db M mu',real=True,positive=True)
 check('invariant_mass_area',sp.simplify(sp.exp(2*t*b)*(M*sp.exp((1-t)*b))**2-M**2*sp.exp(2*b))==0)
 integral=sp.integrate(b*(R-2*t*db),(t,0,1))
 check('integrated_bulk_curvature',sp.simplify(integral-(b*R-b*db))==0)
 K,bn=sp.symbols('K bn')
 # Integration by parts contributes +b*bn/(24pi), exactly cancelling
 # the normal-derivative part of the integrated boundary coefficient.
 bdy=-sp.integrate(b*(K+t*bn)/(12*sp.pi),(t,0,1))+b*bn/(24*sp.pi)
 check('bulk_boundary_cross_term_cancels',sp.simplify(bdy+b*K/(12*sp.pi))==0)
 check('rank_two_a2_bulk',sp.simplify((2*R/6-R/2-2*M**2)/(4*sp.pi)+R/(24*sp.pi)+M**2/(2*sp.pi))==0)
 check('balanced_a2_boundary',sp.simplify((2*K/3-K)/(4*sp.pi)+K/(12*sp.pi))==0)
 A,Ad,V=sp.symbols('A Ad V',commutative=False)
 dA=-(V*A+A*V)/2;dAd=-(V*Ad+Ad*V)/2
 dP=sp.expand(dAd*A+Ad*dA)
 check('two_partner_square_variation',sp.expand(dP+(V*Ad*A+Ad*A*V)/2+Ad*V*A)==0)
 cL,cM,h,v,hp,vp,Q=sp.symbols('cL cM h v hp vp Q')
 u,w=sp.symbols('u w');shape=u*h+w*v;deriv=u*hp+w*vp
 density=cL*(R*shape+Q*deriv**2)+cM*shape*sp.exp(2*shape)
 H=sp.diff(density,u,w).subs({u:0,w:0})
 check('normalized_conversion_Hessian',sp.simplify(H-(2*cL*Q*hp*vp+4*cM*h*v))==0)
 check('product_subtracted_conversion_contact',sp.simplify(H-2*(-2*cL+cM)*h*v-(2*cL*(Q*hp*vp+2*h*v)+2*cM*h*v))==0)
 mass=M**2*sp.exp(2*b)/(4*sp.pi)*(sp.log(M**2/mu**2)+2*b-1)
 conv=-M**2*b*sp.exp(2*b)/(2*sp.pi)
 geo=M**2*sp.exp(2*b)/(4*sp.pi)*(sp.log(M**2/mu**2)-1)
 check('flat_density_identity_all_b',sp.simplify(mass+conv-geo)==0)
 probe=lambda f:sp.simplify((sp.diff(f,b,2)-2*sp.diff(f,b)).subs(b,0))
 check('flat_mass_contact_exact',sp.simplify(probe(mass)-M**2/sp.pi)==0)
 check('flat_conversion_contact_exact',sp.simplify(probe(conv)+M**2/sp.pi)==0)
 check('flat_geometric_product_contact_zero',probe(geo)==0)
 mom=sp.symbols('p',nonnegative=True);kappa=sp.symbols('kappa',positive=True)
 gradient=-M**2/sp.sqrt(M**2+mom**2)
 mass_H=sp.diff(gradient,M)*M
 integrand=sp.simplify(mass_H-2*gradient)
 check('flat_frequency_tadpole_contact',sp.simplify(integrand-M**4/(M**2+mom**2)**sp.Rational(3,2))==0)
 primitive=M**2*mom/sp.sqrt(M**2+mom**2)
 check('continuous_thermal_contact_primitive',sp.simplify(sp.diff(primitive,mom)-integrand)==0)
 check('continuous_thermal_contact_full_tail',sp.limit(primitive,mom,sp.oo)==M**2)
 rational_tail=1/(sp.sqrt(1+kappa**2)*(sp.sqrt(1+kappa**2)+kappa))
 check('continuous_thermal_contact_tail_rationalization',sp.simplify(1-kappa/sp.sqrt(1+kappa**2)-rational_tail)==0)
 n,N,q=sp.symbols('n N q',integer=True,positive=True)
 count=sp.summation(2*(q+2*n),(n,1,N))
 moment=sp.summation(2*(q+2*n)*n*(n+q),(n,1,N))
 check('finite_angular_count',sp.factor(count-2*N*(N+q+1))==0)
 check('finite_angular_mass_moment',sp.factor(moment-N*(N+1)*(N+q)*(N+q+1))==0)
 x=sp.symbols('x',positive=True)
 f=2*(q+2*x)/(x*x*(x+q)**2)
 check('angular_tail_primitive',sp.simplify(f+2*sp.diff(1/(x*(x+q)),x))==0)
 # The derivative numerator is negative with strictly positive coefficients.
 check('angular_tail_monotonicity',sp.simplify(sp.diff(f,x)+4*(q*q+3*q*x+3*x*x)/(x**3*(x+q)**3))==0)
 z=sp.symbols('z',positive=True);aa=sp.Function('a')(z);ww=sp.Function('w')(z);vv=-z**4
 bh=aa/2-z*ww-z*z*sp.diff(ww,z)/2
 weak=aa*vv/2+z*z*sp.diff(vv,z)*ww/2-sp.diff(z*z*vv*ww/2,z)
 check('actual_shape_pairing_integration_by_parts',sp.simplify(bh*vv-weak)==0)
 # Integral z^k*(-log(1-z)-z) = H_(k+1)/(k+1)-1/(k+2).
 mom=lambda k:sp.harmonic(k+1)/(k+1)-sp.Rational(1,k+2)
 check('log_moment_four',mom(4)==sp.Rational(29,100))
 check('weighted_log_moment_four',mom(4)-mom(5)==sp.Rational(103,4200))
 check('positive_weight_moment',sp.integrate(z**4*(1-z),(z,0,1))==sp.Rational(1,30))
 s=sp.symbols('s');a0,a1,a2=sp.symbols('a0 a1 a2');L=sp.symbols('L')
 # Gamma(s) cancellation leaves rational Mellin factors; derive their
 # derivatives at zero rather than independently summing mode finite parts.
 term0=a0*M**2*sp.exp(s*(2*L-sp.log(M**2)))/(8*sp.pi*(s-1))
 term1=a1*sp.exp(s*(2*L-sp.log(M**2)))/(8*sp.pi)
 term2=a2*s*sp.exp(s*(2*L-sp.log(M**2)))/(8*sp.pi*M**2)
 check('one_mass_geometric_a0_polynomial',sp.simplify(sp.diff(term0,s).subs(s,0)-a0*M**2*(sp.log(M**2)-2*L-1)/(8*sp.pi))==0)
 check('one_mass_geometric_a1_polynomial',sp.simplify(sp.diff(term1,s).subs(s,0)+a1*(sp.log(M**2)-2*L)/(8*sp.pi))==0)
 check('one_mass_geometric_a2_polynomial',sp.simplify(sp.diff(term2,s).subs(s,0)-a2/(8*sp.pi*M**2))==0)
 check('four_error_allocations',sp.Rational('2.5e-43')+sp.Rational('2.5e-43')+sp.Rational('5e-43')==sp.Rational('1e-42'))
 check('signed_plus_error_target',sp.Rational('1.277e-40')+sp.Rational('1e-42')==sp.Rational('1.287e-40'))
def execute():
 p,h,zh,env=prior.scalar.prev.york.actual_parameters()
 NUMERICS.extend(prior.scalar.prev.york.NUMERICS)
 B2=p['B2'];C=p['C'];q=arb(p['q']);lam=p['g'];Fl=1-384*arb.pi()*arb('1e22')*lam
 zc=1-arb(2)**-49;A=arb('3.768e-12');a=arb('4e-28');f=arb('4.34e-11');sr=arb('3.33e-10')*h
 ell=-(1-zc).log()-zc
 w=sr+2*(1-zc)/zc*(A*ell+f+a)
 endpoint=zc**6*w/2
 pair=A*(arb(29)/200+4*arb(103)/4200)+4*f/30+a*(arb(1)/10+arb(4)/30)+sr/3+endpoint
 check('actual_pairing_bound',pair<arb('6.704e-12'),str(pair))
 check('finite_cut_in_certified_bulk',(1-zc)/zc>arb('1e-15'))
 check('finite_cut_endpoint_retained',endpoint>0 and endpoint<arb('1e-23'))
 NUMERICS.append(dict(kind='actual_contact_pairing_envelope',B_squared=encoded(B2),C=encoded(C),q=str(p['q']),h=encoded(h),F_lambda=encoded(Fl),cold_cut=encoded(zc),
  endpoint_upper=encoded(endpoint),pairing_upper=encoded(pair),remaining_projection_target='1.287e-40'))
 rows=[]
 for N in [1,2,16,128]:
  moment=arb(N)*(N+1)*(N+q)*(N+q+1)
  coefficient=2*moment/(arb.pi()*C*B2)
  allowance=coefficient*pair/Fl;factor=allowance/arb('1.287e-40')
  threshold=arb('2.5e-43')*B2*Fl*N*(N+q)/(4*C*C)
  kappa=(allowance/(2*arb('2.5e-43'))).sqrt()
  check('separate_contact_allowance_fails_'+str(N),factor>1,str(factor))
  check('positive_angular_tail_threshold_'+str(N),threshold>0)
  rows.append(dict(N=N,finite_mass_moment=encoded(moment),contact_coefficient=encoded(coefficient),
   separate_contact_projection_allowance=encoded(allowance),allowance_target_factor=encoded(factor),
   required_full_time_B_hv_upper=encoded(threshold),
   sufficient_kappa_if_flat_contact_were_truncated=encoded(kappa),
   correlated_flat_contact_after_exact_combination='0',actual_signed_projection_evaluated=False))
 NUMERICS.append(dict(kind='finite_cutoff_contact_strategy_test',cutoffs=rows,physical_source_failure_proved=False))
 heat_factor=B2*Fl/(4*C);threshold=arb('2.5e-43')*heat_factor
 Aexample=threshold/arb('.001')
 check('positive_heat_cap_threshold',threshold>0)
 NUMERICS.append(dict(kind='required_differentiated_remainder_thresholds',physical_heat_factor=encoded(heat_factor),
  cap_error_allocation='2.5e-43',required_A_times_t0_power_over_alpha_upper=encoded(threshold),
  example_alpha='1',example_t0='0.001',required_A_hv_example_upper=encoded(Aexample),
  A_hv_known=False,alpha_known=False,B_hv_known=False,thermal_subtraction_matched=False,complete_projection_evaluated=False))
 print(json.dumps(dict(event='conversion_contact_strategy',pairing_upper=str(pair),first_mode_separate_allowance=str(prior.decode(rows[0]['separate_contact_projection_allowance'])),
  alpha1_t0_001_A_threshold=str(Aexample),full_projection_evaluated=False)),flush=True)
def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],dependencies=deps,invocation=sys.argv,previous_goal_turn='progress')
 try:
  check('previous_projection_components_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',prior.scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  finite_mode_conversion_algebra_verified=status=='passed',continuum_argument_scope='Ordinary note, not certified by checker passes.',
  complete_signed_projection_evaluated=False,thermal_subtraction_matched=False,differentiated_infinite_tail_proved=False,
  actual_source_failure_proved=False,total_source_improvement_certified=False,nonlinear_solution_constructed=False,
  physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
