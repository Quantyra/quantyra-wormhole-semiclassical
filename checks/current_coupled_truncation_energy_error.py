"""Supporting arithmetic for complete coupled transport error.

No retained signed center, physical correction or mass-column certificate.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1];STEM='current-coupled-truncation-energy-error'
OWN=['notes/'+STEM+'.md','checks/current_coupled_truncation_energy_error.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-geodesic-two-sided-remainder.md',
 'notes/current-global-far-return-energy-bound.md',
 'notes/current-coincident-dirac-action-assembly.md',
 'notes/current-global-mass-conversion.md','evidence/current-global-mass-conversion.json',
 'notes/current-adaptive-dirac-trace-remainder.md','evidence/current-adaptive-dirac-trace-remainder.json',
 'notes/current-original-graph-dirac-tail.md','evidence/current-original-graph-dirac-tail.json',
 'notes/current-canonical-interior-schur-bound.md','evidence/current-canonical-interior-schur-bound.json',
 'notes/current-coupled-input-energy-bound.md','notes/current-mass-response-regulator-audit.md']
CHECKS=[]
def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(55)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed';return r['result']
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    p=prior('current-global-mass-conversion')['parameters']
    par=prior('current-adaptive-dirac-trace-remainder')['actual_parameters']
    beta,Fl,u,eps,tau=[dec(p[k]) for k in ['beta','Fl','u','eps','tau']]
    C,B2,gamma=[dec(par[k]) for k in ['C','B2','gamma']]
    check('actual_parameter_signs',beta>arb('.1666') and u>0 and tau>0 and abs(eps)<arb('2e-29') and gamma>0 and gamma<arb('1e-30'))
    z=arb('.99')
    Q99=(1-z)**2*(beta*(1+2*z+3*z*z)+u*z*z+eps*z**4)-tau*z**3
    check('actual_root_lower_endpoint',Q99>0)
    check('actual_root_upper_endpoint',-tau<0)
    check('Q_monotone_on_final_half',12*beta>4*abs(eps))
    check('outer_half_Q_lower',beta/4-abs(eps)-tau>arb('.04'))
    check('inner_factored_H1_lower',beta-32*abs(eps)>arb('.15') and 3*beta+u+eps>0)
    check('outer_half_dominates_bolt_linear_lower',arb('.04')>gamma/2)
    CL=(500/(arb('.99')*Fl)).sqrt()
    CD=(100*arb.pi()*B2*beta/C+2*arb(500)**4*CL**2).sqrt()
    CI=CL/arb('.99').sqrt()+4*CD/gamma
    check('global_shape_L2_upper',CL<arb('22.474'))
    check('global_shape_Laplacian_upper',CD<arb('9.454e16'))
    check('actual_global_shape_supremum_upper',CI<arb('7.477e50'))
    q,s=arb('.3961'),arb('.401');rw=q**128;rs=s**128
    ep0=arb('9.204e-77')*CI**2/beta
    epsilon=ep0*rw;source=arb('9.204e-79')*rw*CI/beta
    check('old_energy_upper_certificate_excessive',ep0>arb('3.0865e26') and ep0<arb('3.0866e26'))
    check('new_rule_dominates_weak_rule',8/(1/s).log()>4/(1/q).log() and 6/(1/s).log()>2/(1/q).log())
    check('base_198_not_certified_by_this_upper',ep0*q**70>arb('.01'))
    check('base_199_certified_by_this_upper',ep0*q**71<arb('.01'))
    check('complete_energy_defect',epsilon<arb('1.020e-25'))
    check('complete_source_dual_error',source<arb('1.365e-78'))
    K=arb(7);Kt=K/(1-K*epsilon);H=arb('1.338').sqrt();F=H*arb('1.1e-10')
    check('full_Neumann_inverse_admissible',epsilon<1)
    check('interior_Neumann_inverse_admissible',K*epsilon<1 and Kt<arb('7.000001'))
    W4=arb('1.354e-53')*rs
    check('direct_graph_fourth_jet_error',W4<arb('2.159e-104'))
    rows=[('H3','6.861e-45','4.583e-13','1.094e-95','2.977e-36','2.914e-46'),
     ('quartic','1.571e-51','7.949e-26','2.505e-102','5.314e-49','5.221e-59'),
     ('mass','1.571e-52','7.949e-27','2.505e-103','5.314e-50','5.221e-60'),
     ('slow','1.896e-51','1.141e-22','3.023e-102','7.410e-46','7.255e-56')]
    outputs={}
    for tag,a,c,at,bt,ft in rows:
        aa=arb(a)*rs;cc=arb(c);dd=aa*CI
        error=dd*K*H+(cc+dd)*Kt*epsilon*(1+K*H)
        affine=aa+dd*K*F+(cc+dd)*Kt*(source+K*epsilon*F)
        check(tag+'_direct_original_output_error',aa<arb(at))
        check(tag+'_complete_interior_product_error',error<arb(bt))
        check(tag+'_complete_affine_error',affine<arb(ft))
        outputs[tag]={'direct_U3':enc(aa),'interior_product_per_energy':enc(error),'affine_total':enc(affine),'interior_tail_per_energy':enc(dd)}
    # Exact algebra supports, rather than replaces, the ordinary form proof.
    zS,v,ga,hh=S.symbols('z v gamma H1',positive=True)
    tt=1/zS;LL=1/v
    Qfact=(tt-v)*(2*ga*v+(tt-v)*hh)/tt**3
    check('exact_bolt_linear_term',S.simplify(Qfact-2*ga*v*v*zS*zS*(LL-zS)-(tt-v)**2*hh/tt**3)==0)
    r=S.symbols('r',positive=True)
    check('integrable_regular_bolt_derivative',S.integrate(2/(ga*S.sqrt(r)),(r,0,LL))==4*S.sqrt(LL)/ga)
    B,A,P,Pt,dB,dA,dD=S.symbols('B A P Pt dB dA dD',commutative=False)
    raw=B*P*A-(B-dB)*Pt*(A-dA)
    expected=dB*P*A+(B-dB)*(P-Pt)*A+(B-dB)*Pt*dA
    check('complete_noncommuting_Schur_difference',S.expand(raw-expected)==0)
    check('retained_fast_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'parameters':{'beta':enc(beta),'C':enc(C),'B2':enc(B2),'gamma':enc(gamma),'Fl':enc(Fl)},
     'global_shape_L2_constant':enc(CL),'global_shape_Laplacian_constant':enc(CD),'global_shape_supremum_constant':enc(CI),
     'old_base_128_energy_upper_expression':enc(ep0),'minimum_base_for_point01_upper_certificate':199,
     'new_spectral_order_rule':'256+ceil((8log(M/18)+6log(1+abs(p)))/log(1/.401))',
     'complete_energy_defect':enc(epsilon),'complete_source_dual_error':enc(source),
     'approximate_full_energy_inverse_upper':enc(1/(1-epsilon)),'approximate_interior_energy_inverse_upper':enc(Kt),
     'strong_tail_factor':enc(rs),'weak_tail_factor':enc(rw),'direct_graph_W4_coefficient':enc(W4),'outputs':outputs,
     'scope':{'actual_paired_mass_log_domain':True,'common_regulator_and_complete_spectrum_preserved':True,
      'original_energy_to_shape_supremum_proved':True,'complete_coupled_energy_approximation_error_certified':True,
      'all_interior_Schur_error_terms_certified':True,'original_outer_and_phase_error_bounds_certified':True,
      'signed_retained_center_computed':False,'strong_complementary_quotient_inverse_computed':False,
      'actual_input_or_border_map_norms_computed':False,'actual_kappa_m_computed':False,
      'complete_coupled_feasibility_bound_certified':False,'nonlinear_admission':False,'nonlinear_contraction':False,
      'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
      'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'energy_error':run.get('result',{}).get('complete_energy_defect'),
       'graph_error':run.get('result',{}).get('outputs',{}).get('H3'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
