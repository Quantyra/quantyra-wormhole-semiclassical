"""Analytic scalar trial and full original-domain weak residual enclosure."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc
from flint import acb

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-scalar-green-source-residual'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = ['notes/'+STEM+'.md', 'checks/current_scalar_green_source_residual.py',
          'checks/current_quartic_ritz_rank_test.py',
          'notes/current-first-inverse-source-screen.md',
          'notes/current-complete-bolt-energy-bound.md',
          'notes/current-seed-full-hessian-realization.md',
          'notes/product-centered-scalar-cold-response.md',
          'notes/current-global-green-trial-enclosure.md',
          'notes/current-filtered-mixed-green-residual.md',
          'evidence/current-global-green-trial-enclosure.json',
          'evidence/current-first-inverse-source-screen.json',
          'evidence/current-filtered-mixed-green-residual.json',
          'evidence/current-complete-bolt-energy-bound.json',
          'evidence/current-coupled-mass-adjoint-runtime.json']


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def prior(name):
    d = json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
    r = d.get('runs', [d])[-1]
    assert r['status'] not in ['failed', 'unresolved']
    return r['result']


def execute():
    ctx.prec = 512
    green = prior('current-global-green-trial-enclosure')
    old = prior('current-first-inverse-source-screen')['constants']
    old = {k:dec(v) for k,v in old.items()}
    filt = prior('current-filtered-mixed-green-residual')['constants']
    filt = {k:dec(v) for k,v in filt.items()}
    p = {k:dec(v) for k,v in green['trial_parameters'].items()}
    L,a,gamma,beta,d = [p[k] for k in ['L','a','gamma','beta','d']]
    At = dec(green['constants']['A_t'])
    S,SP,M3 = arb('1e-12'),arb('1e-16'),arb(1000)
    zm = L-2*S
    u,eps = p['u'],p['epsilon_seed']
    tau0 = (u+3*beta+eps)*L**4-eps*(5*L**5-arb(9)/2*L**6)
    taumax = arb('3.01')*beta
    taujet = 4*abs(u+3*beta+eps)+52*abs(eps)
    assert 0 < tau0 < taumax and taujet < 3 and L <= 1
    muY,muS = 4*a/3,4*beta/L**2
    c = 2*tau0/(3*L**2)
    degree = lambda mu: ((1+4*mu/a).sqrt()-1)/2
    ks = [p['ell'],p['nu'],degree(muY),degree(muS)]
    lam = [a*k*(k+1) for k in ks]
    assert ks[0] >= 1 and ks[1] > ks[0] and arb('.5') < ks[2] < 1
    assert 0 < ks[3] < arb('.5')
    t = [At,-At,arb(0),arb(0)]

    def resolvent(v,mu,index):
        ans = [arb(0)]*4
        for j in range(4):
            if j == index:
                assert v[j].is_zero()
                continue
            assert not (mu-lam[j]).contains(0)
            ans[j] = v[j]/(mu-lam[j])
        ans[index] = -sum(ans,arb(0))
        assert sum(ans,arb(0)).contains(0)
        for j in range(4):
            assert ((mu-lam[j])*ans[j]-v[j]).contains(0)
        return ans

    qY = resolvent(t,muY,2)
    xi = [-2*t[j]/3-muY*qY[j]/3 for j in range(4)]
    y = [c*v for v in resolvent(xi,muS,3)]
    assert sum(xi,arb(0)).contains(0)
    assert sum((v*l for v,l in zip(y,lam)),arb(0)).contains(0)
    xS = a*S/gamma
    X = 2*xS
    assert xS > 1
    B2 = (arb(3).log()**2/3+2*arb(3).log()/9+arb(2)/27)/4

    def tail_int(k):
        if k > arb('.5'):
            return 1/(2*k-1)
        assert k < arb('.5')
        return (X**(1-2*k)-1)/(1-2*k)

    def bounds(v):
        out = dict(N=arb(0),M=arb(0),H=arb(0),T=arb(0),U=arb(0))
        for aj,k in zip(v,ks):
            if aj.is_zero():
                continue
            aa = abs(aj)
            I = acb(k+1).polygamma(1).real/(2*k+1)
            tail = tail_int(k)
            out['N'] += aa*(gamma/a*I).sqrt()
            out['M'] += aa*((gamma/a)**3*(B2+arb.pi()**2*tail)).sqrt()
            out['H'] += aa*gamma*(gamma/a*(1+arb.pi()**2*(k+1)**2*tail)).sqrt()
            out['T'] += aa*arb.pi()*xS**(-k-1)
            out['U'] += aa*arb.pi()*gamma*(k+1)*xS**(-k)
        return out

    bt,bq,bx,by = [bounds(v) for v in [t,qY,xi,y]]
    r = M3*2*S/(6*a)
    Q2S = 2*gamma*S+a*S*S
    Q22S = 4*gamma*S+4*a*S*S
    Qmax = (1+r)*Q22S
    Qsmax = 2*gamma+4*a*S+M3*(2*S)**2/2
    assert r < arb('1e-8') and zm > arb('.99')
    assert 2*(2*a-2*M3*S)/3 > 1
    YN = 2*bq['H']/(3*L**2)
    YsM = 2*(muY*bq['M']+bt['M'])/(3*L**2)
    YT = 2*bq['U']/(3*L**2)
    YsT = 2*(muY*bq['T']+bt['T'])/(3*L**2)
    Dg = 16/(3*zm)
    Ds = (6*(1+r)*(2*gamma+2*a*S)/zm**2+2*M3/3
          +8*(2*a+M3*S)/(3*zm)+muY*M3/(6*a))
    york = dict(
        derivative_coefficient=4*(1+r).sqrt()*bt['H']/(3*L*zm**2*a.sqrt()),
        potential=(Dg*gamma*2*bt['N']/(3*L**2*muY.sqrt())+Ds*YN/a.sqrt())/(1-r).sqrt(),
        source_value=8*(1+r).sqrt()*((2*gamma).sqrt()*(bt['N']*bt['M']).sqrt()+a.sqrt()*bt['M'])/(3*zm**3),
        cutoff=(Qmax*S).sqrt()*(4*bt['T']/(3*zm**2*S)+4*YsT/S+15*YT/S**2))
    eY = sum(york.values(),arb(0))
    erec = YN/zm+L*YsM/zm**2+YT/S.sqrt()
    exi = (arb(1)/2+1/arb(5).sqrt())*eY+erec
    CE = 2*L+L**2*M3/(6*a)
    CB = 5*L**2*M3/(6*a)+2*L*(1+r)
    scalar = dict(
        actual_York_return=2*taumax*exi/3,
        curvature_variation=2*bx['M']/zm**2,
        geometry=(CE*(muS*by['M']+c*bx['M'])+CB*by['H'])/zm**2,
        cutoff=S.sqrt()/zm**2*(L**2*Qmax*(15*by['T']/S**2+4*by['U']/(S*Q2S))
                               +(L**2*Qsmax+2*L*Qmax)*2*by['T']/S))
    rP = sum(scalar.values(),arb(0))
    epsbar = p['C']/(60*arb.pi()*p['B2'])
    Cv = 1+11*epsbar*beta/3
    S0 = 1/arb('23.99').sqrt()
    S1 = (beta/arb('5.99')).sqrt()
    S2 = (beta/(6*d)).sqrt()
    SH = (S2*S2+13*beta*S1*S1).sqrt()
    w0,w2 = arb('1.1e-27')*beta,arb('4.3e-25')*beta**2
    KA = (2*d+11*epsbar/6)*w0+22*epsbar*taumax/3
    Gy = (1+r).sqrt()*((c*by['N']*bx['N']).sqrt()+2*(Q22S/S).sqrt()*by['T'])/zm
    EA = (KA*S1*Gy+2*w0*S0*by['N']/zm**2)/beta
    Rjets = 3*arb('1e6')*(4400+2*4*2*arb(2).sqrt()*20+4*2)
    Cjets = 100*Rjets
    Kcurv = arb(22).sqrt()*(100*Cjets*(SH+2*S1+S0)+arb('1e6')*S0)
    Kcurv_old = arb('1e20')*(SH+S1+S0)
    assert Rjets < arb('1e12') and Cjets < arb('1e14') and Kcurv < Kcurv_old
    brems = dict(
        curvature_density_contact=4*taumax*11*epsbar*w0*S0/36,
        scalar_curvature_Hessian=arb(22).sqrt()*d*w2*S0,
        full_curvature_anomaly=epsbar*Kcurv,
        magnetic_contact=20/(60*arb.pi())*arb('4e-27')*beta*S0,
        full_off_shell_return=arb('9e-28')*beta*5/beta.sqrt()*S1)
    Kbrem = sum(brems.values(),arb(0))
    xiH = filt['K_V']*bt['N']/zm**2
    Eb = Kbrem*xiH/beta
    Eb_old = (Kbrem-epsbar*Kcurv+epsbar*Kcurv_old)*xiH/beta
    inputdiff = At*gamma**2/(arb(27).sqrt()*a**2*SP**(arb(3)/2)*zm**2)
    charges = dict(
        accepted_input_to_trial_input=old['scalar_source_uniform_upper']*inputdiff,
        complete_core_and_cutoff=6*(Cv*S0+d*S2)*rP/beta,
        complete_scalar_remainder=EA,
        complete_mixed_remainder=Eb)
    rho = sum(charges.values(),arb(0))
    rho_old = rho-Eb+Eb_old
    Kscalar = old['first_scalar_residual_to_original_bolt_multiplier']
    physical = Kscalar*rho
    old_physical = Kscalar*rho_old
    subtotal = old['known_partial_subtotal_not_free_margin']+physical
    bolt = arb('1.667e-9')
    # The output records the outcome without redefining its target on failure.
    status = ('complete_scalar_source_residual_enclosed_signed_quantum_response_unresolved'
              if physical < bolt else 'complete_scalar_source_estimator_excessive_signed_response_unresolved')
    constants = dict(muY=muY,muS=muS,kY=ks[2],kS=ks[3],tau0=tau0,
                     cutoff_radius=S,relative_geometry_guard=r,
                     York_residual_energy_upper=eY,York_reconstruction_error=erec,
                     complete_xi_difference_Hmu_upper=exi,rP_Hmu_upper=rP,
                     curvature_tensor_jet_sum=Rjets,unrounded_curvature_adjoint=Kcurv,
                     old_curvature_adjoint=Kcurv_old,curvature_estimator_gain=Kcurv_old/Kcurv,
                     complete_mixed_remainder_constant=Kbrem,
                     full_scalar_residual_Es_dual_upper=rho,
                     previous_tensor_rounding_scalar_residual=rho_old,
                     scalar_residual_to_original_bolt_multiplier=Kscalar,
                     complete_scalar_source_physical_debit=physical,
                     old_rounding_complete_source_physical_debit=old_physical,
                     debit_fraction_of_whole_allowance=physical/bolt,
                     old_rounding_report_over_whole_allowance=old_physical/bolt,
                     known_partial_subtotal_not_free_margin=subtotal,
                     unspent_allowance_before_unknown_terms_not_certified_margin=bolt-subtotal)
    assert all(v.is_finite() for v in constants.values())
    return status, dict(
        precision_bits=ctx.prec,
        constants={k:enc(v) for k,v in constants.items()},
        profiles={name:[{'degree':enc(k),'coefficient':enc(v)} for k,v in zip(ks,vs)]
                  for name,vs in [('t0',t),('qY',qY),('xi0',xi),('y0',y)]},
        profile_bounds={name:{k:enc(v) for k,v in bs.items()}
                        for name,bs in [('t0',bt),('qY',bq),('xi0',bx),('y0',by)]},
        York_residual_parts={k:enc(v) for k,v in york.items()},
        scalar_second_order_residual_parts={k:enc(v) for k,v in scalar.items()},
        complete_mixed_remainder_parts={k:enc(v) for k,v in brems.items()},
        complete_scalar_dual_residual_parts={k:enc(v) for k,v in charges.items()},
        original_output_debits={k:enc(Kscalar*v) for k,v in charges.items()},
        parameter_use='correlated exact definitions enclosed; no independent midpoint replacement',
        scope=dict(actual_paired_form_domain=True,actual_global_York_return_charged=True,
                   complete_scalar_source_residual=True,independent_quantum_shape_application=False,
                   complete_signed_Zr=False,full_bolt_admission=False,other_coupled_budgets=False,
                   nonlinear_admission_contraction=False,exact_compatibility=False,
                   physical_matching=False,semiclassical_validity_noise=False,
                   all_sector_stability=False,construction_obstruction=False,full_goal_complete=False),
        unknowns=['signed b(y_hat,v)/N-e_t(x_hat,v) over the entire original E_t',
                  'subsequent independent quantum applications and their numerical errors',
                  'correlated j_r contribution and complete Z_r',
                  'other coupled source/graph/physical/mass/phase entries and nonlinear neighborhood'])


if __name__ == '__main__':
    output = ROOT/'evidence'/(STEM+'.json')
    if output.exists():
        raise RuntimeError('Preserve immutable evidence; select a new output before rerunning.')
    before = {p:digest(ROOT/p) for p in INPUTS}
    protected = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p)) == h for p,h in protected.items())
    rt = json.loads((ROOT/INPUTS[-1]).read_bytes())
    assert digest(Path(sys.executable)) == rt['python_sha256']
    assert all(digest(Path(pkg['root'])/p) == h for pkg in rt['packages'].values() for p,h in pkg['files'].items())
    record = dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
                  python_sha256=digest(Path(sys.executable)),
                  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
                  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        record['status'],record['result'] = execute()
        after = {p:digest(ROOT/p) for p in INPUTS}
        assert before == after
        assert all(digest(Path(p)) == h for p,h in protected.items())
        record.update(input_sha256_after=after,protected_paths_unchanged=len(protected))
    except Exception:
        import traceback
        record.update(status='failed',traceback=traceback.format_exc())
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    output.write_bytes((json.dumps(record,indent=2)+'\n').encode('utf-8'))
    cs = record.get('result',{}).get('constants',{})
    keys = ['full_scalar_residual_Es_dual_upper','complete_scalar_source_physical_debit',
            'old_rounding_complete_source_physical_debit','curvature_estimator_gain',
            'known_partial_subtotal_not_free_margin']
    print(json.dumps(dict(status=record['status'],error=record.get('traceback'),
                         constants={k:cs[k]['ball'] for k in keys if k in cs}),indent=2))
    if record['status']=='failed':
        raise SystemExit(1)
