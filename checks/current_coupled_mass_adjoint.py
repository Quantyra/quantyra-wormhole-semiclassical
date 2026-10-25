"""Current coupled mass-flux algebra and conditional budgets, NOT a solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import flint
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-coupled-mass-adjoint'
OWN = ['notes/'+STEM+'.md', 'checks/current_coupled_mass_adjoint.py',
       'evidence/'+STEM+'-review.json']
INPUTS = OWN + [
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-coupled-outer-return-test.md',
    'notes/current-coupled-input-energy-bound.md',
    'notes/current-global-return-energy-attempt.md',
    'notes/current-physical-york-dual-assembly.md',
    'notes/current-seed-boundary-phase-and-trace-defect.md',
    'notes/recentered-reference-inverse.md',
    'notes/recentered-mass-boundary.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-global-weighted-map.md',
    'notes/current-global-mass-conversion.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-global-mass-conversion.json',
]
CHECKS = []

def check(name, value):
    CHECKS.append({'name': name, 'passed': bool(value)})
    if not value:
        raise AssertionError(name)

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(35), 'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def zero_matrix(matrix):
    return all(sp.simplify(v)==0 for v in matrix)

def execute():
    ctx.prec=512
    s,be,Cv,r,D0=sp.symbols('s beta Cv r D0')
    N=sp.Function('N')
    aa=lambda t:(t+3)/(3*(t+1))
    J=-6*be*(Cv-r*s*(s-3))*(s-4)*(s+1)
    A=sp.Matrix([[J,J*aa(s)],
                 [aa(3-s)*J,aa(3-s)*J*aa(s)+4*D0*N(s)]])
    p=sp.Matrix([-1,1]); vm=sp.Matrix([-sp.Rational(1,2),1])
    check('current_normal_root_zero_column',zero_matrix((A.subs(s,0)*p).subs(N(0),0)))
    check('current_normal_mass_column',zero_matrix((A.subs(s,3)*vm).subs(N(3),0)))
    check('correct_reflected_left_null',zero_matrix((p.T*A.subs(s,3)).subs(N(3),0)))
    check('current_normal_adjoint_reflection',zero_matrix((A.T-A.subs(s,3-s)).subs(N(3-s),N(s))))
    flux=sp.simplify((p.T*A.diff(s)*vm)[0].subs(s,3))
    expected=4*D0*sp.Subs(sp.Derivative(N(s),s),s,3)
    check('complete_scalar_and_mixed_flux_cancellation',sp.simplify(flux-expected)==0)
    wrong=(p.T*A.subs(s,0)).subs(N(0),0)
    check('prior_wrong_weight_test_is_not_used',not zero_matrix(wrong))

    g,Fl,eg,km,L=sp.symbols('g Fl EulerGamma kappa log_mu')
    Fstar=Fl+23*g/(180*sp.pi)
    h3=1-eg+km
    d0=1/(480*sp.pi**2); AT=Fstar/(48*sp.pi*g)
    nd=3*(2*d0*h3-AT)
    K=sp.simplify(4*(-8*sp.pi*g*be)*nd/be)
    Zc=Fl+g/(5*sp.pi)*(L+eg+sp.Rational(151,180))
    check('complete_flux_equals_twice_physical_charge',
          sp.simplify(K.subs(km,-sp.Rational(6,5)-L)-2*Zc)==0)

    G,C,B,gamma,nh=sp.symbols('G C B gamma nh',positive=True)
    Wabs=nh*sp.pi*B**2/(G*gamma)
    betaT=2*sp.pi*B/(C*gamma)
    check('physical_action_radial_measure',sp.simplify(Wabs/(nh*4*sp.pi*betaT)-C*B/(8*sp.pi*G))==0)
    lam=C*be/B**2
    factor=8*sp.pi*G/(lam*B**3)/(nh*4*sp.pi*betaT)
    check('energy_action_conversion',sp.simplify(factor-1/(Wabs*be))==0)

    # The abstract residual identity is checked algebraically; these symbols
    # are NOT numerical physical matrix entries or a finite-mode surrogate.
    Fp,Fva,qUa,hvaUa=sp.symbols('Fp Fva qUa hvaUa')
    ia1=Fp-Fva+hvaUa-qUa
    ia2=Fp-qUa+hvaUa-Fva
    check('two_primal_dual_centers_identical',sp.expand(ia1-ia2)==0)

    receipt=json.loads((ROOT/'evidence/current-global-mass-conversion.json').read_bytes())
    prev=receipt['runs'][-1]
    check('parameter_source_was_successful_acceptance',prev['status']=='passed' and prev['mode']=='acceptance')
    for f in ['notes/current-global-mass-conversion.md']:
        check('parameter_source_argument_unchanged',sha256((ROOT/f).read_bytes()).hexdigest()==prev['input_sha256_before'][f])
    pp=prev['result']['parameters']
    ag=dec(pp['g']); afl=dec(pp['Fl'])
    # The entire inherited kappa interval is retained, not a fitted contact.
    ak=arb('9.3','1.5')
    euler=arb.const_euler()
    ah=1-euler+ak
    Knum=2*(afl+23*ag/(180*arb.pi()))-2*ag*ah/(5*arb.pi())
    Klo=arb('1.9999999764')
    check('complete_flux_positive_lower',Knum>Klo)
    check('complete_flux_upper_below_two',Knum<2)
    alog=-arb(6)/5-ak
    Znum=afl+ag/(5*arb.pi())*(alog+euler+arb(151)/180)
    check('numerical_full_contact_charge_agrees',Knum.overlaps(2*Znum))
    check('physical_charge_guard_preserved',Znum>arb('.9999999882'))

    mass_target=arb('2.29e-15')
    Ihat_cap=arb('3.3e-15'); eval_cap=arb('1e-16')
    residual_cap=arb('1e-5'); source_cap=arb('1.1e-10')
    numerator=Ihat_cap+eval_cap+residual_cap*source_cap
    mass_cap=numerator/Klo
    available_product=Klo*mass_target-Ihat_cap-eval_cap
    residual_threshold=available_product/source_cap
    graph_cap=arb('5e-12')+arb(2174)*mass_cap
    check('conditional_residual_product_fits',residual_cap*source_cap<available_product)
    check('conditional_mass_certificate_fits_split',mass_cap<mass_target)
    check('conditional_mass_rounded_guard',mass_cap<arb('2.250000027e-15'))
    check('conditional_graph_conversion_fits_original_radius',graph_cap<arb('1e-11'))
    check('allowable_adjoint_residual_exceeds_rounded_target',residual_threshold>arb('1.0727e-5'))
    check('unchanged_fast_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'normal_mass_flux_symbolic':str(flux),
            'physical_flux_symbolic':str(K),'physical_flux_equals':'2 Z_c',
            'complete_mass_flux':enc(Knum),'physical_charge':enc(Znum),
            'used_mass_flux_lower':enc(Klo),
            'conditional_flux_center_allowance':enc(Ihat_cap),
            'conditional_flux_evaluation_error_allowance':enc(eval_cap),
            'conditional_adjoint_residual_allowance':enc(residual_cap),
            'known_forward_source_dual_upper':enc(source_cap),
            'conditional_residual_product_ceiling':enc(available_product),
            'conditional_adjoint_residual_ceiling':enc(residual_threshold),
            'conditional_mass_upper':enc(mass_cap),
            'conditional_original_graph_upper':enc(graph_cap),
            'scope':{'current_fixed_seed_mass_adjoint_argument':True,
                     'actual_mass_computed':False,'actual_signed_pairing_computed':False,
                     'actual_adjoint_residual_computed':False,'actual_kappa_m_computed':False,
                     'complete_physical_coupled_return':'unresolved',
                     'physical_feasibility_test_executed':False,
                     'nonlinear_admission':False,'construction_obstruction':False,
                     'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
         'input_sha256_before':hashes,
         'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
         'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                    'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary={'mode':mode,'status':run['status'],'checks':len(CHECKS),
             'result':run.get('result'), 'traceback':run.get('traceback')}
    message=json.dumps(summary,indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
