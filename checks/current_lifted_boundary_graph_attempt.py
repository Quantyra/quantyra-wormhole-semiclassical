"""Arithmetic supporting the lifted mass/global-graph attempt, not a global quantum solver."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-lifted-boundary-graph-attempt'
OWN = ['notes/current-lifted-boundary-graph-attempt.md', 'checks/current_lifted_boundary_graph_attempt.py', 'evidence/current-lifted-boundary-graph-attempt-review.json']
INPUTS = ['notes/current-lifted-boundary-graph-attempt.md', 'checks/current_lifted_boundary_graph_attempt.py', 'evidence/current-lifted-boundary-graph-attempt-review.json', 'evidence/current-lifted-boundary-graph-attempt-exploration.json', 'notes/current-global-return-energy-attempt.md', 'notes/current-global-boundary-moment-return.md', 'notes/current-coupled-local-graph-bound.md', 'notes/current-complete-strong-response-attempt.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-finite-local-action-assembly.md', 'notes/current-symmetric-contact-feasibility.md', 'notes/current-controlled-solution-audit-closeout.md', 'notes/current-curved-paired-heat-bridge.md', 'evidence/current-global-boundary-moment-return.json', 'evidence/current-coupled-outer-return-test.json', 'evidence/current-coupled-local-graph-bound.json', 'evidence/current-complete-strong-response-attempt.json', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
CHECKS = []


def check(name, condition, detail=""):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def encoded(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(45), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def decode(value):
    p = value['decimal_enclosure']
    return arb(int(p['mid']), int(p['rad'])) * arb(10) ** int(p['exponent'])


from math import comb, factorial
sys.path.insert(0, str(ROOT/'checks'))
# A legacy parameter helper parses argv on import. Keep this wrapper's
# acceptance mode local while importing it, then restore the exact argv.
original_argv = sys.argv[:]
sys.argv = [v for v in sys.argv if v != '--acceptance']
import current_physical_york_dual_assembly as york
sys.argv = original_argv
york.check = check

def execute():
    ctx.prec=512
    a=arb
    rr=json.loads((ROOT/'evidence/current-global-boundary-moment-return.json').read_bytes())['runs'][-1]['result']
    dec=decode; enc=encoded
    e=dec(rr['exact_slow_exponent_offset']); pole=dec(rr['full_normal_shape_particular']) if 'full_normal_shape_particular' in rr else None
    prior=json.loads((ROOT/'evidence/current-coupled-outer-return-test.json').read_bytes())['runs'][-1]['result']
    B0=dec(prior['full_normal_shape_particular']); pole=B0.abs_upper()*e.abs_upper()
    g0=dec(rr['g0']); bp=dec(rr['Bprime_at_slow_root']); k3=dec(rr['known_source_H3_upper'])
    sig=a('4.005'); aa=a('.05'); eo=a('1.1e-10')
    co=1/(2*(a('.005')-e)).sqrt(); cf=1/(2*(a('.5')+e)).sqrt()
    mk=g0.abs_upper()*(1+a(10).log())*(e*(1+a(10).log())).exp()
    d0=2*pole+(mk+cf*a('1.02e-46'))/bp
    de=co/bp; dm=cf*a('6.50e-39')/bp
    md=(a(3)/5).sqrt()*aa
    me=1801*(a(3)/2).sqrt()/(2*sig-3).sqrt()*aa**(sig-3)
    ma=a(3).sqrt()*aa**a('-1.5')*a('2.518')*eo+pole*md*((-aa.log())**2+2*(-aa.log())/5+a(2)/25).sqrt()
    m0=(ma+md*d0+me*k3)/(1-md*dm); m1=(me+md*de)/(1-md*dm)
    etest=a('1e-11'); mm=m0+m1*etest; delta=d0+de*etest+dm*mm
    out={'mass_intercept':enc(m0),'mass_slope':enc(m1),'mass_at_target':enc(mm),'delta_at_target':enc(delta),'phase_A':enc(6*mm+delta+6*pole),'far_at_target':enc(a('1.02e-46')+a('6.50e-39')*mm)}
    pp,h,zh,env=york.actual_parameters(); B2=pp['B2']; C=pp['C']; pi=a.pi()
    mapping={'beta':pp['beta'],'u':pp['u'],'eps':pp['eps'],'tau':pp['tau'],'ell':a(2000)**-10*B2/C,'mag':pi*pp['q']**2*a(2000)**8/(C*B2),'cr':8*pi*a('1e22')*C/B2}
    def val(x):
        if x.is_Symbol:return mapping[str(x)]
        if x.is_Rational:return a(int(x.p))/int(x.q)
        if x.is_Add:return sum((val(v) for v in x.args),a(0))
        if x.is_Mul:
            r=a(1)
            for v in x.args:r*=val(v)
            return r
        if x.is_Pow:return val(x.base)**int(x.exp)
        raise ValueError(x)
    u,s=sp.symbols('u s'); cut=1-35*u**4+84*u**5-70*u**6+20*u**7
    def integralmoment(rate, pol):
        cs=sp.Poly(pol,u); top=cs.degree(); J=[(rate.exp()-1)/rate]
        for j in range(1,top+1):J.append((rate.exp()-j*J[-1])/rate)
        return sum((val(cs.nth(j))*J[j] for j in range(top+1)),a(0))
    def inner(k,l):
        pk=cut; pl=cut; kk=val(k); ll=val(l); v=a(0)
        for j in range(4):
            if j:pk=sp.expand(sp.diff(pk,u)+k*pk); pl=sp.expand(sp.diff(pl,u)+l*pl)
            v+=comb(3,j)*(kk**j*ll**j/(kk+ll)+integralmoment(kk+ll,sp.expand(pk*pl)))
        return a('.01')**(kk+ll)*(-(kk+ll)).exp()*v
    data=json.loads((ROOT/'evidence/current-coupled-local-graph-bound.json').read_bytes())['runs'][-1]['result']
    masscols=[]
    for i in range(2):
        terms=[]
        for term in data['matrix_bounds'][i][1]['terms']:
            coef=sp.sympify(term['symbol'], locals={k:sp.Symbol(k) for k in mapping}).subs(s,3)/4
            terms.append((sp.Rational(term['z_power']+3)-sp.Rational(801,200),val(coef)))
        normsq=sum((ci*cj*inner(ki,kj) for ki,ci in terms for kj,cj in terms),a(0))
        norm=normsq.sqrt();masscols.append(norm)
        out['local_bare_mass_row_'+str(i)]={'unit_H3':enc(norm),'energy_phase_allowance':enc(norm*mm),'threshold_mass_for_source':enc(a('7.49e-12')/norm),'coefficients':[{'power':str(k),'coef':enc(v)} for k,v in terms]}
    out['required_energy_reduction_phase_only']=enc(max(masscols)*mm/a('7.49e-12'))

    from math import factorial
    R=a('.1'); Rext=a('.2'); aa=a('.05'); ac=a('.04'); acore=a('.02')
    q=[a(1),a(0),pp['u']/pp['beta'],(-4*pp['beta']-2*pp['u']-pp['tau'])/pp['beta'],(3*pp['beta']+pp['u']+pp['eps'])/pp['beta'],-2*pp['eps']/pp['beta'],pp['eps']/pp['beta']]
    def qbounds(radius):
        dq=sum((abs(q[k])*radius**k for k in range(1,7)),a(0))
        d1=sum((k*abs(q[k])*radius**k for k in range(1,7)),a(0))
        return dq,d1,dq/(1-dq)+d1/(4*(1-dq)**2)
    out['mass_profile_complex_guards']={str(r):[enc(v) for v in qbounds(r)] for r in [aa,R,Rext]}
    check('zero_free_actual_profile_disks',qbounds(aa)[2]<a('.001') and qbounds(R)[2]<a('.008') and qbounds(Rext)[2]<a('.07'))
    out['profile_A_coefficient_upper']=enc(max(a('1.02')*R**3*factorial(j)/(R-ac)**j*(ac if j==4 else 1) for j in range(5)))
    check('mass_profile_nine_Euler_jets',all(a('1.02')*a('1.5').exp()*factorial(j)*2**j <= 6*40**j for j in range(10)))
    # Reclose the mass/slow/far system in the lifted phase decomposition.
    dm_new=dm*a('1.02'); den=1-a('.001')-md*dm_new
    mn0=(ma+md*d0+me*k3)/den;mn1=(me+md*de)/den
    mn=mn0+mn1*etest;dn=d0+de*etest+dm_new*mn
    out['lifted_boundary_block']={'mass_intercept':enc(mn0),'mass_slope':enc(mn1),'mass_at_target':enc(mn),'Delta_at_target':enc(dn),'phase_A_at_target':enc(100*mn+dn+6*pole),'Efar_at_target':enc(a('1.02e-46')+a('6.63e-39')*mn)}
    # Quantify the cost of treating the analytic lift correction as fast graph data.
    N=32
    inv=[a(1)]
    for n in range(1,N+1):inv.append(-sum((q[j]*inv[n-j] for j in range(1,min(6,n)+1)),a(0)))
    invsq=[sum((inv[j]*inv[n-j] for j in range(n+1)),a(0)) for n in range(N+1)]
    J={n+4:invsq[n]/(n+4) for n in range(N+1)}
    bcoef={n+3:inv[n] for n in range(N-2)}
    for j in range(1,7):
        for n,v in J.items():
            p=n+j-1
            if p<=N:bcoef[p]=bcoef.get(p,a(0))+j*q[j]*v
    bcoef[3]-=1
    KL=2*pp['beta']*(pp['Fl']+23*pp['g']/(180*a.pi()))/3
    def rawjet(n):
        k=a(n)-sig
        return acore**k*(1+k*k)**a('1.5')/(2*k).sqrt()
    lead=abs(KL*bcoef[6]*36)*rawjet(6)
    others=sum((abs(KL*coef*(n-4)*n*(n-3))*rawjet(n) for n,coef in bcoef.items() if n>=5 and n!=6),a(0))
    n=N+1;ratio=acore/R*(a(n+1)/n)**6
    tail=abs(KL)*a('1.02')*R**3*acore**(-sig)*a(n)**6*(acore/R)**n*3/(1-ratio)
    # psi has the beta(11,11) smoothstep on log(.02)..log(.04).
    wexp=a(5)-sig
    caps=[]
    for j in range(10):
        cap=a('1e9')*a('.5').exp()*ac**wexp/(2*wexp).sqrt()*sum((comb(j,l)*(21/a(2).log())**l*factorial(j-l)*2**(j-l) for l in range(j+1)),a(0))
        caps.append(cap)
    h9=sum((comb(9,j)*caps[j]**2 for j in range(10)),a(0)).sqrt()
    quantum=6000*pp['g']*pp['beta']/(15*a.pi())*h9
    out['lift_original_graph_cost']={'classical_leading_raw_jets':enc(lead),'other_series_raw_jets_upper':enc(others),'analytic_tail_upper':enc(tail),'profile_H9_upper':enc(h9),'complete_normal_quantum_upper':enc(quantum),'graph_cost_lower_per_mass':enc(lead-others-tail-quantum),'allowed_mass_interval_cost_at_upper':enc((lead-others-tail-quantum)*mn)}
    # Uniform admitted-family extension: only the finite phase H9 guard changes.
    bg=json.loads((ROOT/'evidence/current-complete-strong-response-attempt.json').read_bytes())['runs'][-1]['result']
    parts=bg['background_parts']; oldB=dec(bg['Bcap'])
    deltaB=512/a(2).log()*dec(parts['Wxi'])*(2+512*dec(parts['Wa']))*dec(parts['slow_H8log_over_h2'])*(a('.02')/a('.00105681')-1)
    out['extended_phase_Bcap']=enc(oldB+deltaB)
    out['extended_phase_C3_cap_allowance']=enc((oldB+deltaB)*a('8.123e-26'))
    # Direct analytic-coefficient majorant for -Ipot''+Isq'' applied to the lift.
    R=a('.1'); cr=mapping['cr']; eps=abs(pp['eps']); dqmag=a('3.4e-22')*pp['beta'];fl=pp['Fl'];beta=pp['beta']
    sumshift=[a(2),a(8),a(36),a(188)]
    nu=[2*R**3*v for v in sumshift];bb=nu[:]
    Qj=[pp['beta']*sum((abs(q[k])*k**j*R**k for k in range(7)),a(0)) for j in range(4)]
    # k^0 is 1 for k=0 in Python, retaining Q's constant.
    Hj=[eps*(4*3**j*R**3+6*4**j*R**4) for j in range(5)]
    def add(*vals):return [sum((v[j] for v in vals),a(0)) for j in range(min(map(len,vals)))]
    def scale(c,v):return [c*x for x in v]
    def mul(v,w):return [sum((comb(j,l)*v[l]*w[j-l] for l in range(j+1)),a(0)) for j in range(min(len(v),len(w)))]
    def D(v):return v[1:]
    def shift(p,v):return [R**p*sum((comb(j,l)*abs(p)**(j-l)*v[l] for l in range(j+1)),a(0)) for j in range(len(v))]
    A=scale(2*cr,mul(nu,D(Hj)))
    Hi=add(scale(2,mul(bb,Hj)),mul(nu,D(Hj)))
    lapA=shift(-2,add(mul(Qj,add(D(D(A)),D(A))),mul(D(Qj),D(A))))
    QA=mul(Qj,A)
    termgrad=mul(mul(Hj,Qj),D(nu))
    fn=add(scale(6,lapA),scale(12,shift(-2,add(D(QA),QA))),scale(24*beta,shift(-2,A)),scale(96*cr*beta,shift(-2,mul(Hj,add(nu,bb)))),scale(24*cr,shift(-2,add(D(termgrad),termgrad))))
    fb=add(scale(2,lapA),scale(4*pp['u'],A),scale(24*beta,shift(-2,A)),scale(cr,add(scale(4,mul(bb,mul(Hj,Hj))),scale(4,mul(Hj,Hi)))),scale(16*cr*pp['u'],mul(Hj,bb)),scale(96*cr*beta,shift(-2,mul(Hj,add(nu,bb)))))
    poly=lambda c1,c2:[abs(c1)*R+abs(c2)*2**j*R**2 for j in range(2)]
    fn=add(fn,scale(fl*eps,add(mul(poly(4,12),nu),mul(poly(8,12),bb))))
    fb=add(fb,scale(4*fl*dqmag,bb),scale(fl*eps,add(mul(poly(8,12),nu),mul(poly(16,12),bb))))
    Fn=shift(4,fn);Fb=shift(4,fb)
    T=Fb[1]+4*Fb[0]+(Fn[1]+6*Fn[0])/3
    W=max(Fn[0],T);kk=a(7)-sig
    Mpoly=lambda k:k**3+12*k*k+51*k+100
    localmass=W*R**-7*a('.01')**kk*Mpoly(kk)/(2*kk).sqrt()
    out['lifted_complete_gravitational_mass_column']={'analytic_coefficient_norm':enc(W),'outer_H3_upper_per_mass':enc(localmass),'source_at_mass_upper':enc(localmass*mn),'mass_feedback_slope':enc(localmass*mn1)}

    # Ordinary-argument guards, all full-series or exact arithmetic.
    check('conditional_mass_original',mm<a('4.341e-8'))
    check('bare_mass_positive_allowance_fails',max(masscols)*mm/a('7.49e-12')>a('34.5'))
    check('bare_mass_required_ceiling',a('7.49e-12')/max(masscols)<a('1.25720e-9'))
    check('lifted_mass_feedback_denominator',den>a('.998'))
    check('lifted_mass_ceiling',mn<a('4.345e-8'))
    check('lifted_slow_ceiling',dn<a('2.251e-10'))
    check('lifted_far_ceiling',a('1.02e-46')+a('6.63e-39')*mn<a('3.901e-46'))
    check('lifted_collar_A',100*mn+dn+6*pole+a('5.681e-5')<a('6.116e-5'))
    check('profile_A_guard',dec(out['profile_A_coefficient_upper'])<100)
    check('quartic_slow_phase_jet_guard',a('.001')+a('.00005681')+100*a('.001')/6+a('.001')<a('.02'))
    check('quartic_slow_ODE_recurrence',(1+a('4.01')/40)**4-1<1)
    check('classical_tail_ratio',ratio<1)
    check('classical_series_tail',tail<a('3.16e-11'))
    check('known_profile_old_graph_cost',lead-others-tail-quantum>a('.009068'))
    check('known_profile_old_graph_cost_at_allowed_mass',(lead-others-tail-quantum)*mn>a('3.93e-10'))
    check('complete_normal_quantum_cost',quantum<a('4.54e-11'))
    x=sp.symbols('x')
    smooth=1-sp.integrate(x**10*(1-x)**10,(x,0,x))*sp.factorial(21)/sp.factorial(10)**2
    coeffsum=sum(abs(v) for v in sp.Poly(smooth,x).all_coeffs())
    check('C10_cutoff_coefficient_guard',coeffsum<10**9)
    check('C10_cutoff_endpoints',smooth.subs(x,0)==1 and smooth.subs(x,1)==0 and all(sp.diff(smooth,x,j).subs(x,v)==0 for j in range(1,11) for v in [0,1]))
    are=a('3.005');tc=(are**2-4)/4
    check('whole_digamma_low_frequency_guard',(are/2).log()+(tc/are).atan()+1/are<2)
    check('whole_digamma_high_frequency_guard',a.pi()/2+1/are<2)
    check('whole_cotangent_guard',a.pi()/2/(a.pi()*a('.005')/2).tan()<200)
    check('whole_normal_polynomial_guard',a('4.005')*a('1.005')*a('3.005')*a('2.005')<25)
    check('whole_normal_bracket_guard',200+a('10.8')+3<222 and 25*222<6000)
    check('lifted_phase_cap',oldB+deltaB<a('2.062e30'))
    check('lifted_phase_C3_cap',(oldB+deltaB)*a('8.123e-26')<167500)
    check('conditional_quantum_fast_screen',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    check('full_gravitational_series_monotone',a('.1')*(1+1/(7-sig))**3<1)
    check('full_gravitational_coefficient_norm',W<a('1.089e-27'))
    check('full_gravitational_mass_column',localmass<a('1.762e-24'))
    check('conditional_gravitational_source',localmass*mn<a('7.66e-32'))
    check('conditional_gravitational_feedback',localmass*mn1<a('8.714e-23'))
    target=a('5e-5'); remainS=a('7.49e-12')-target*mn0;remainEta=a('.00368')-target*mn1
    check('specified_mass_map_target_affine',target*mn0<a('2.148e-12'))
    check('specified_mass_map_target_feedback',target*mn1<a('.002474'))
    check('specified_other_global_source_budget',remainS>a('5.342e-12'))
    check('specified_other_global_feedback_budget',remainEta>a('.001206'))
    # Independent exact curvature-variation identity on an arbitrary local Q.
    zz=sp.symbols('z',positive=True); be,uu=sp.symbols('beta u');QQ=sp.Function('Q')(zz);JJ=sp.Function('J')(zz)
    nn=-2*be*QQ*JJ/zz;bbm=be*zz**3/QQ+be*sp.diff(QQ,zz)*JJ
    lap=lambda f:sp.diff(QQ*sp.diff(f,zz),zz)
    HH=-sp.diff(QQ,zz,2)+6*sp.diff(QQ,zz)/zz-12*QQ/zz**2+2*uu+12*be/zz**2
    HHm=-2*lap(bbm)-6*lap(nn)+12*QQ*sp.diff(nn,zz)/zz+4*uu*bbm+24*be*(nn+bbm)/zz**2
    identity=HHm-2*bbm*HH+nn*zz*sp.diff(HH,zz)
    for j in [3,2,1]:identity=identity.subs(sp.diff(JJ,zz,j),sp.diff(zz**3/QQ**2,zz,j-1))
    check('exact_arbitrary_Q_mass_curvature_variation',sp.cancel(identity)==0)
    out.update(outcome='unresolved_complete_global_graph_estimate',bare_positive_estimator_failed=True,
     conditional_quantum_fast_screen_upper='.135973',conditional_quantum_fast_screen_certified=True,
     complete_quantum_matching_map_computed=False,actual_graph_regularity_proved=False,
     mass_map_target=enc(target),conditional_affine_debit=enc(target*mn0),conditional_feedback_debit=enc(target*mn1),
     other_global_source_remaining=enc(remainS),other_global_feedback_remaining=enc(remainEta),
     actual_controlled_solution=False,physical_X_comparisons_proved=False,full_goal_complete=False,
     physical_domain_changed=False,normalization_changed=False,next_calculation_executed=False)
    return out


if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    import flint
    flint_root = Path(flint.__file__).resolve().parent
    runtime = {'python': sys.version, 'python_executable_sha256': sha256(Path(sys.executable).read_bytes()).hexdigest(),
               'flint_version': flint.__version__, 'sympy_version': sp.__version__,
               'flint_file_sha256': {str(p.relative_to(flint_root)): sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(flint_root.rglob('*'))
                                    if p.is_file() and p.suffix in ['.py', '.pyd', '.dll']}}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(), 'mode': mode,
           'input_sha256_before': hashes, 'precision_bits':512, 'runtime':runtime,
           'frozen_sources_utf8': {p: (ROOT/p).read_text(encoding='utf-8') for p in OWN}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('all_frozen_inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs': []}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger, indent=2)+'\n').encode())
    log = ROOT/'evidence'/(STEM+'.log')
    text = json.dumps({'mode': mode, 'status': run['status'], 'checks': len(CHECKS),
                       'result': run.get('result'), 'traceback': run.get('traceback')}, indent=2)+'\n'
    with log.open('ab') as f:
        f.write(text.encode())
    print(text)
    if run['status'] != 'passed':
        raise SystemExit(1)
