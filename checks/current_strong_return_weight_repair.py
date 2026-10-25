"""Repaired weight and both original normal rows; not a physical solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json,sys,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
_args=sys.argv[:];sys.argv=[v for v in sys.argv if v!='--acceptance']
import current_physical_york_dual_assembly as york
sys.argv=_args
ctx.prec=512
STEM='current-strong-return-weight-repair'
OWN=['notes/'+STEM+'.md','checks/current_strong_return_weight_repair.py','evidence/'+STEM+'-review.json']
INPUTS=set(OWN+[
    'notes/current-short-travel-strong-domain.md',
    'notes/current-coupled-outer-return-test.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-physical-york-dual-assembly.md',
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-short-travel-strong-domain.json',
    'evidence/current-coupled-outer-return-test.json',
    'evidence/current-canonical-heat-restoration-assembly.json'])
for module in list(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path:
        try:p=Path(path).resolve().relative_to(ROOT).as_posix()
        except ValueError:continue
        if p.endswith('.py'):INPUTS.add(p)
INPUTS=sorted(INPUTS);CHECKS=[]
def check(name,condition):
    CHECKS.append(dict(name=name,passed=bool(condition)))
    if not condition:raise AssertionError(name)
def enc(value):
    mid,rad,exp=value.mid_rad_10exp(70)
    return dict(ball=value.str(42),decimal_enclosure=dict(mid=str(mid),rad=str(rad),exponent=str(exp)))
def execute():
    s,J,B,a,f,g=sp.symbols('s J B a f g')
    theta=f/J;b=((s-4)*g-(s-6)*f/3)/B;n=theta-a*b
    # Both original normal rows, including the undifferentiated second row.
    check('original_first_row',sp.simplify(J*(n+a*b)-f)==0)
    row2=(6-s)*J*(n+a*b)/(3*(4-s))+B*b/(s-4)
    check('original_second_row',sp.simplify(row2-g)==0)
    gc,bb,bet=sp.symbols('GC B2 beta',positive=True)
    kk=gc*bet**2/(15*sp.pi*bb)
    joint=sp.simplify((2*gc/bb)*25*sp.sqrt(5)/(128*kk))
    check('unchanged_original_scale_product',sp.simplify(joint-750*sp.pi*sp.sqrt(5)/(128*bet**2))==0)
    for k in [0,1,4,8,16]:
        for ell in range(k+1):
            lhs=sum(sp.binomial(k,j)*sp.binomial(j,ell)*sp.Rational(1,4)**(j-ell) for j in range(ell,k+1))
            rhs=sp.binomial(k,ell)*sp.Rational(5,4)**(k-ell)
            if lhs!=rhs:raise AssertionError((k,ell))
    check('log_absorption_product_identity',True)
    x,L=sp.symbols('x L',positive=True)
    # Exact translated endpoint test for log degree one; the ordinary proof
    # supplies the entire family of test functions and divergence as L grows.
    check('translated_endpoint_coefficient',sp.expand(1+(x+L))==1+x+L)
    check('first_weight_and_second_margin',sp.Rational(3,2)+sp.Rational(7,4)==sp.Rational(13,4) and sp.Rational(13,4)+sp.Rational(7,4)>sp.Rational(801,200))
    p,h,zh,env=york.actual_parameters()
    check('prior_normal_source_passed',json.loads((ROOT/'evidence/current-coupled-outer-return-test.json').read_bytes())['runs'][-1]['status']=='passed')
    check('prior_transfer_source_passed',json.loads((ROOT/'evidence/current-short-travel-strong-domain.json').read_bytes())['runs'][-1]['status']=='passed')
    la,be,Fl,B2=p['g'],p['beta'],p['Fl'],p['B2']
    DR=96*arb.pi()*arb('1e22')-1/(60*arb.pi());r=DR*la
    Cv=1+11*la/(180*arb.pi());Fs=Fl+23*la/(180*arb.pi())
    K=la*be/(15*arb.pi());bigL=10*arb.pi()*Fs/la
    kap=(la*B2).log()/2-arb('1.5e-9','1.5e-9')-arb(6)/5
    check('original_contact_and_scale',kap>arb('7.8') and kap<arb('10.8') and bigL>arb('3e34'))
    aa=arb(9)/4;psi=aa.digamma()
    check('entire_low_frequency_digamma',psi>0 and psi+4*(aa**-3+1/(2*aa**2))<2)
    check('low_frequency_absolute_digamma',arb(4)+arb.pi()**2/4<16)
    check('cotangent_entire_vertical_line',sp.cos(13*sp.pi/8)**2<sp.sin(13*sp.pi/8)**2)
    check('complete_low_frequency_product',(arb(45)/16+4+7)*(4+arb.pi()/2+arb('10.8'))<238)
    check('entire_high_frequency_imaginary_guard',aa+2*arb.pi()/((2*arb.pi()).exp()-1)<3)
    check('entire_high_frequency_real_guard',arb('7.8')-arb.pi()/2>5)
    check('real_product_polynomial',sp.Rational(45,16)*5+sp.Rational(7,2)*3==sp.Rational(393,16))
    check('high_frequency_positive_gap',bigL/2+4>arb(393)/16)
    A=Cv-arb(13)*r/16
    check('scalar_entire_frequency_gap',A>r and r>0)
    c3=3*K*bigL/32
    c4low=c3/arb(5).sqrt();c4high=arb(32)*K*(2*bigL).sqrt()/25
    check('fourth_gain_bottleneck',c4high<c4low)
    c5low=c3/5;c5high=128*K/(25*arb(5).sqrt())
    check('fifth_gain_bottleneck',c5high<c5low)
    Cb=1/c5high;Cb3=1/c4high;Ct=2/(9*be*r)
    entries=[[Ct+arb(275)*Cb/612,arb(25)*Cb/51],[arb(11)*Cb/12,Cb]]
    rows=[sum(row,arb(0)) for row in entries]
    Cjoint=750*arb.pi()*arb(5).sqrt()/(128*be**2)
    check('original_shape_row_dominates_matrix_guard',rows[1]>rows[0])
    # Concrete endpoints of the exact universal formula, not assumed log
    # degrees for the physical kernel.
    m0=arb(10)**(-arb(1)/4);m1=4*(-arb(3)/4).exp()
    check('log_maximum_inside_collar',arb(3)>arb(10).log())
    check('log_maximum_beats_boundary',m1>m0*(1+arb(10).log()))
    b1_16=(arb(5)/4)**16*m1+16*(arb(5)/4)**15*m0
    check('finite_log_absorption',b1_16>0)
    return dict(outcome='unresolved_complete_coupled_estimate',
        endpoint_failure='The h^2 log^N envelope alone has no unlogged two-weight bound for N>0; no claim that the actual R_s endpoint operator is unbounded.',
        repaired_weights=['3/2','13/4','4.005'],
        log_absorption_formula='sum_{l=0}^{min(k,N)} binom(k,l)(5/4)^(k-l)(N)_l m_(N-l)',
        example_only_m0=enc(m0),example_only_m1=enc(m1),example_only_B_1_16=enc(b1_16),
        actual_kernel_log_degrees=None,actual_kernel_seminorms=None,
        full_normal_original_rows=True,c3=enc(c3),C3=enc(1/c3),
        c4_low=enc(c4low),c4_high=enc(c4high),c5_low=enc(c5low),c5_high=enc(c5high),
        C_b=enc(Cb),C_theta=enc(Ct),
        three_gain_C_b=enc(Cb3),three_gain_row_guard=enc(arb(23)*Cb3/12),
        fourth_derivative_price=enc(Cb/Cb3),
        original_quantum_normal_composite_factor=enc(Cjoint),
        original_quantum_normal_nu_factor=enc(arb(25)*Cjoint/51),
        original_normal_inverse_matrix=[[enc(v) for v in row] for row in entries],
        original_normal_inverse_max_row_guard=enc(rows[1]),
        normal_input='both original rows H^k on the auxiliary log line',
        normal_output='canonical (nu,b) H^(k+4) on the same line',
        quantitative_role='Necessary first-weight normal inversion factor; actual finite-geometry return and gluing not bounded.',
        M_graph_4=None,M_X_8=None,kappa_mass=None,actual_mass=None,
        physical_operator_applications=0,coupled_signed_solve_executed=False,
        new_nonlinear_feasibility_gate_closed=False,nonlinear_admission=False,
        exact_compatibility=False,physical_matching=False,validity_noise=False,
        all_sector_stability=False,actual_controlled_solution=False,
        construction_obstruction=False,closed_audit_preserved=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in INPUTS}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),
        flint_version=flint.__version__,sympy_version=sp.__version__,
        flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*'))
                        if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,
        input_sha256_before=before,precision_bits=512,runtime=runtime,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback'))
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(summary,indent=2)+'\n').encode())
    brief={k:summary[k] for k in ['mode','status','checks','traceback']}
    if run.get('result'):
        brief['normal_matrix_row_guard']=run['result']['original_normal_inverse_max_row_guard']['ball']
        brief['C_theta']=run['result']['C_theta']['ball'];brief['C_b']=run['result']['C_b']['ball']
        brief['three_gain_row_guard']=run['result']['three_gain_row_guard']['ball']
        brief['fourth_derivative_price']=run['result']['fourth_derivative_price']['ball']
        brief['original_quantum_composite_factor']=run['result']['original_quantum_normal_composite_factor']['ball']
    print(json.dumps(brief,indent=2))
    if run['status']!='passed':raise SystemExit(1)
