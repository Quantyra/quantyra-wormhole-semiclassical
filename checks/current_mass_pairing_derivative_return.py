"""Whole-cutoff jets and complete massive/LLL quantum mass pairing bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import argparse
import json
import sys
import traceback

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_mass_pairing_integrable_return as old
m=old.m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
fr=old.fr
J=fr.J
ROOT,BASELINE,runtime=m.ROOT,m.BASELINE,m.runtime
STEM='current-mass-pairing-derivative-return'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
up=lambda x:arb(abs(x).upper())
maxup=lambda a,b:arb(max(up(a).upper(),up(b).upper()))


def interval(lo,hi):
    lo,hi=arb(lo.lower()),arb(hi.upper())
    return arb((lo+hi)/2,(hi-lo)/2)


def polynomial(coeff,x):
    out=arb(0)
    for c in reversed(coeff): out=out*x+c
    return out


def cutoff_cells(n=1024):
    sp=m.sp
    t=sp.Symbol('t');a=t*(1-t)
    P=sp.Integer(1);polys=[]
    for j in range(4):
        poly=sp.Poly(sp.expand(P),t)
        polys.append([arb(int(poly.nth(k))) for k in range(poly.degree()+1)])
        P=sp.expand(a*a*sp.diff(P,t)+sp.diff(a,t)*(1-2*j*a)*P)
    cells=[]
    for i in range(n):
        t0,t1=arb(i)/n,arb(i+1)/n
        avals=[t0*(1-t0),t1*(1-t1)]
        alo=arb(min(x.lower() for x in avals))
        ahi=arb(1)/4 if 2*i<=n<=2*(i+1) else arb(max(x.upper() for x in avals))
        gl=arb(0) if i==0 or i==n-1 else (-1/alo).exp()
        gh=(-1/ahi).exp()
        cells.append(dict(t0=t0,t1=t1,alo=alo,ahi=ahi,glo=gl,ghi=gh))
    Zlo=sum((c['glo']/n for c in cells),arb(0))
    Zhi=sum((c['ghi']/n for c in cells),arb(0))
    Z=interval(Zlo,Zhi)
    assert Zlo>0 and Zhi<1
    cumulative_lo,cumulative_hi=arb(0),arb(0)
    out=[]
    for i,c in enumerate(cells):
        T=interval(c['t0'],c['t1'])
        lo=1-(cumulative_hi+c['ghi']/n)/Zlo
        hi=1-cumulative_lo/Zhi
        lo=arb(max(arb(0).lower(),lo.lower()))
        hi=arb(min(arb(1).upper(),hi.upper()))
        chi=[interval(lo,hi)]
        for j in range(4):
            def at(a):
                return arb(0) if not a>0 else a**(-2*j)*(-1/a).exp()
            candidates=[at(c['alo']),at(c['ahi'])]
            if j and c['alo']<=arb(1)/(2*j)<=c['ahi']:
                candidates.append(at(arb(1)/(2*j)))
            qlo=arb(min(v.lower() for v in candidates[:2]))
            qhi=arb(max(v.upper() for v in candidates))
            derivative=interval(qlo,qhi)*polynomial(polys[j],T)/Z
            chi.append(-(arb(10)**(j+1))*derivative/factorial(j+1))
        out.append((arb('.1')+c['t0']/10,arb('.1')+c['t1']/10,J(chi)))
        cumulative_lo+=c['glo']/n
        cumulative_hi+=c['ghi']/n
    assert cumulative_lo.overlaps(Zlo) and cumulative_hi.overlaps(Zhi)
    return out,Z,[str(p) for p in polys]


def mass_cell(z0,z1,chi,coeff,beta):
    z=J([interval(z0,z1),arb(1)]+[arb(0)]*4)
    Q=J([arb(0)]*6)
    for c in reversed(coeff): Q=Q*z+c
    qp=Q.d()
    z5=J(z.c[:5]);q5=J(Q.c[:5])
    g=z5**3/(q5*q5)
    primitive=interval(z0**4/(4*arb('.17')**2),z1**4/(4*arb('.15')**2))
    jm=J([primitive]+[g.c[k-1]/k for k in range(1,5)])
    b=beta*z5**3/q5+beta*qp*jm
    h=-arb('1.5')*chi*b
    hd=[factorial(j)*h.c[j] for j in range(5)]
    q,q1,q2,q3=[factorial(j)*Q.c[j] for j in range(4)]
    assert q>arb('.15') and q<arb('.17')
    K1=up(q.sqrt()*hd[1])
    A=q*hd[2];B=q1*hd[1]/2
    H1sq=arb.pi()*up(q)*up(hd[1])**2*(z1-z0)
    H2sq=(3*arb.pi()*up(A)**2/4+2*arb.pi()*up(A*B)+2*arb.pi()*up(B)**2)*(z1-z0)
    K3=up(q**arb('1.5'))*up(hd[3])+arb('1.5')*up(q.sqrt()*q1)*up(hd[2])+up(q.sqrt()*q2)*up(hd[1])/2
    K4=up(q*q)*up(hd[4])+3*up(q*q1)*up(hd[3])+(arb('.75')*up(q1)**2+2*up(q*q2))*up(hd[2])+(up(q1*q2)/4+up(q*q3)/2)*up(hd[1])
    return H1sq,H2sq,K1,K3,K4


def execute():
    ctx.prec=768
    p,geo,bs=fr.load()
    mp={k:dec(v) for k,v in m.result('current-global-mass-conversion')['parameters'].items()}
    beta=p['beta']
    coeff=[beta,arb(0),mp['u'],-4*beta-2*mp['u']-mp['tau'],3*beta+mp['u']+mp['eps'],-2*mp['eps'],mp['eps']]
    cuts,Z,polys=cutoff_cells()
    cells=[(arb(i)/5120,arb(i+1)/5120,J([arb(1)]+[arb(0)]*4)) for i in range(512)]+cuts
    H1sq,H2sq,K1,K3,K4=[arb(0)]*5
    for z0,z1,ch in cells:
        h1,h2,k1,k3,k4=mass_cell(z0,z1,ch,coeff,beta)
        H1sq+=h1;H2sq+=h2
        K1=maxup(K1,k1);K3=maxup(K3,k3);K4=maxup(K4,k4)
    H1,H2=H1sq.sqrt(),H2sq.sqrt()
    t0=arb(1)/32
    Bfrac=((8*arb.pi()/3)*K4*K4*t0**arb('1.5')+512*K3*K3*t0.sqrt()+16*arb.pi()*K3*K3/t0.sqrt()).sqrt()
    prev=m.constants('current-mass-pairing-integrable-return')
    adj=m.constants('current-filtered-mass-adjoint-bound-attempt2')
    bridge=m.constants('current-filtered-test-mass-trace-attempt2')
    green=m.constants('current-filtered-mixed-green-residual')
    v={k:old.ball(old.rational(x)) for k,x in m.result('current-correlated-signed-trial-test-attempt3')['exact_nominal_parameters'].items()}
    assert v['muY']<v['twice_a'] and v['muT']>1 and v['d1']/2>1
    ga,L,a=p['gamma'],p['L'],p['a'];muY=4*a/3
    A0=3*v['alpha']/(4*v['gamma']*v['A'])
    Kf=A0*(1+2*v['scalar_c']/v['muS'])
    L0=Kf*(ga/(2*v['muT'])).sqrt()
    L1=prev['complete_Fref_L1_bound']*L0
    P0=ga/arb('1e-110');Wmax=2*(arb(4)/3).exp()
    X10=Kf*(ga*((1+P0).log()/2+Wmax*P0/(1+P0))).sqrt()
    X20=Kf*(ga*Wmax*P0).sqrt()
    Astar=Kf*v['kappa']/(1-prev['complete_Fref_Neumann_ratio'])
    Cv,delta=v['Cv'],v['delta']
    X11=Astar*(ga*(Cv*Cv/8+Cv*delta/3+delta*delta/4)).sqrt()
    X21=Astar*(ga*(Cv*Cv/6+Cv*delta/2+delta*delta/2)).sqrt()
    frac=m.result('current-source-adapted-fractional-return')
    sc={k:[dec(x) for x in row] for k,row in frac['scalar_coefficient_derivative_bounds'].items()}
    detail=frac['derivative_accounting']
    m1=(2*ga*L+a*L*L).sqrt();m2=arb(2).sqrt()*(ga+a*L)
    Dmax=arb('6.3')+arb('8.1')*L
    d2=arb(max(arb(1).upper(),(4*a).upper())).sqrt()
    r1=1/(2*muY.sqrt())
    Klow=(sc['A'][1]+sc['B'][1]/muY+sc['B'][0]*r1
         +sc['C'][1]*m1*r1+sc['C'][0]*(m2*r1+m1*d2)
         +sc['E'][1]*L*Dmax/muY+sc['E'][0]*Dmax*m1/muY)
    N10=(2*arb.pi()).sqrt()*(sc['A'][0]*X10+Klow*L0)
    N11=(2*arb.pi()).sqrt()*(sc['A'][0]*X11+Klow*L1)
    C2=10*dec(detail['mixed_derivative_constants']['(2, 0)'])
    N20=(2*arb.pi()).sqrt()*C2*X20
    N21=(2*arb.pi()).sqrt()*C2*X21
    eps=arb('1e-120');K=arb('1e13');R0=arb(1)/400
    split=(eps/(4*K)).sqrt()
    assert split<R0
    grad=K*(arb('.25')+2*(R0/split).log())+arb('1e13')+arb('1e7')
    Gh=(eps*H2*H2/3+grad*H1*H1).sqrt()
    G0=(eps*N20*N20/3+grad*N10*N10).sqrt()
    G1=(eps*N21*N21/3+grad*N11*N11).sqrt()
    alpha=prev['quantum_action_alpha_q'];Kg=adj['actual_filtered_mass_functional_norm_upper'];q=arb('.338')
    Egeo=2*alpha*Kg*(q*N10+N11)*(H1+Bfrac)
    Egain=alpha*Kg*Gh*(q*G0+G1)
    ELLL=alpha*Kg*(arb(2000)**4/p['C'])*K1*(q*N10+N11)/(12*arb.pi()*arb.pi().sqrt())
    Ezero=prev['complete_zero_order_mass_block_upper']
    Eshape=Ezero+Egeo+Egain+ELLL
    conv=bridge['full_E_to_retained_Et_norm_comparison'];Eb=bridge['full_to_retained_j_energy_error']
    target=arb('1.49e34')
    remaining=(target-Eb)/conv-Eshape
    vals=dict(cutoff_normalization=Z,mass_shape_first_flow_L2_upper=H1,mass_shape_second_flow_L2_upper=H2,
        mass_shape_first_flow_sup_upper=K1,mass_shape_third_flow_sup_upper=K3,mass_shape_fourth_flow_sup_upper=K4,
        mass_shape_fractional_13_over_4_upper=Bfrac,
        original_source_amplitude=A0,complete_fj_spectral_factor=Kf,
        original_trial0_profile_L2_upper=L0,original_trial1_profile_L2_upper=L1,
        original_trial0_auxiliary_order_one_upper=X10,original_trial1_auxiliary_order_one_upper=X11,
        original_trial0_auxiliary_order_two_upper=X20,original_trial1_auxiliary_order_two_upper=X21,
        actual_return_coefficient_D_bound=Dmax,full_first_derivative_low_order_coefficient=Klow,
        actual_trial0_first_flow_norm_upper=N10,actual_trial1_first_flow_norm_upper=N11,
        actual_order_two_reconstruction_constant=C2,
        actual_trial0_second_flow_norm_upper=N20,actual_trial1_second_flow_norm_upper=N21,
        gain_two_epsilon=eps,gain_two_split=split,complete_gain_two_gradient_coefficient=grad,
        trial0_second_derivative_fraction_of_gain_norm=eps*N20*N20/(3*G0*G0),
        complete_principal_mass_block_upper=Egeo,complete_gain_two_mass_block_upper=Egain,
        complete_LLL_mass_block_upper=ELLL,unchanged_zero_order_mass_block_upper=Ezero,
        complete_massive_and_LLL_mass_pairing_upper=Eshape,
        shape_block_fraction_of_conditional_target=conv*Eshape/target,
        shape_block_original_bolt_charge=green['pref_t']*conv*Eshape,
        conditional_full_E_account_after_bridge_and_shape=remaining,
        complete_full_E_to_retained_Et_comparison=conv,unchanged_full_response_bridge=Eb)
    assert all(x.is_finite() for x in vals.values())
    assert all(x>0 for k,x in vals.items() if k!='conditional_full_E_account_after_bridge_and_shape')
    return dict(constants={k:enc(x) for k,x in vals.items()},whole_cell_count=len(cells),
        cutoff_cell_count=len(cuts),cutoff_derivative_polynomials=polys,
        numerical_contract_pass=bool(Eshape<arb('1e34')),
        complete_conditional_account_has_positive_remainder=bool(remaining>0),
        actual_local_and_scalar_mass_residual=None,actual_scalar_return_bolt_value=None,
        signed_r0_mass_pairing=None,signed_r1_mass_pairing=None,actual_signed_mass_vector_norm=None,
        actual_mass_free_residual=None,complete_application_error=None,actual_Dj=None,actual_signed_Gt=None,
        scope=dict(actual_mass_lift_fractional_constant_bounded=True,
            complete_massive_and_LLL_mass_pairings_bounded=True,
            all_local_Weyl_scalar_and_classical_rows_included=False,
            complete_signed_mass_term_certified=False,complete_signed_response_certified=False,
            complete_bolt_allowance=False,available_physical_margin=False,
            original_P32_defect=False,actual_kappa_m=False,complete_original_admission=False,
            nonlinear_neighborhood_transfer=False,nonlinear_existence=False,exact_compatibility=False,
            physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False,
            construction_obstruction=False,full_goal_complete=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
    out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable receipts'
    prior=json.loads((ROOT/'evidence/current-mass-pairing-integrable-return.json').read_bytes())
    inputs=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_mass_pairing_derivative_return.py',
        'evidence/current-mass-pairing-integrable-return.json',
        'evidence/current-mass-pairing-integrable-return-verification.json',*prior['input_sha256_after']]))
    before={p:digest(ROOT/p) for p in inputs};protected=json.loads(BASELINE.read_bytes())
    rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        assert all(digest(ROOT/p)==h for p,h in prior['input_sha256_after'].items())
        assert all(digest(Path(p))==h for p,h in protected.items())
        assert digest(Path(sys.executable))==runtime['python_sha256']
        for pkg in runtime['packages'].values():
            assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
        rec['result']=execute()
        after={p:digest(ROOT/p) for p in inputs};assert before==after
        assert all(digest(Path(p))==h for p,h in protected.items())
        rec.update(status='massive_and_LLL_pairings_enclosed_full_mass_residual_open',input_sha256_after=after,
            protected_paths_unchanged=len(protected),python_sha256=runtime['python_sha256'],
            pinned_packages={k:len(v['files']) for k,v in runtime['packages'].items()})
    except Exception:
        rec.update(status='failed',traceback=traceback.format_exc())
    rec['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
        numerical_contract_pass=rec.get('result',{}).get('numerical_contract_pass'),
        constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
    if rec['status']=='failed':raise SystemExit(1)
