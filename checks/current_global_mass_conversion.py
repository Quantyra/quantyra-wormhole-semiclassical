"""Full original mass-lift graph conversion; no physical coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import flint
import sympy

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-global-mass-conversion'
OWN = ['notes/'+STEM+'.md', 'checks/current_global_mass_conversion.py',
       'evidence/'+STEM+'-review.json']
INPUTS = OWN + [
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-lifted-boundary-graph-attempt.md',
    'notes/current-resolvent-corrected-angular-remainder.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-coupled-global-return-attempt.json',
    'evidence/current-lifted-boundary-graph-attempt.json',
    'evidence/'+STEM+'-exploration.json',
    'checks/current_physical_york_dual_assembly.py',
    'checks/current_all_angular_principal_operator.py',
    'checks/current_subtracted_thermal_action.py',
    'checks/current_coincident_dirac_action_assembly.py',
    'checks/current_endpoint_thermal_restoration.py',
    'checks/current_seed_radial_green_enclosure.py',
    'checks/current_seed_separated_quantum_response.py',
]
CHECKS = []

def check(name, value, detail=None):
    CHECKS.append({'name': name, 'passed': bool(value)})
    if not value:
        raise AssertionError(name)

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(32), 'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

sys.path.insert(0,str(ROOT/'checks'))
saved=sys.argv[:]
sys.argv=[v for v in sys.argv if v!='--acceptance']
import current_physical_york_dual_assembly as york
sys.argv=saved
york.check=check
york.prior.check=check

def poly_mul(p,q):
    r=[arb(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return r

def horner(p,x):
    v=arb(0)
    for a in reversed(p):v=v*x+a
    return v

def interval(lo,hi):
    return (lo+hi)/2+arb(0,((hi-lo)/2).abs_upper())

def sqrange(x):
    lo=arb(0) if x.contains(0) else arb(abs(x).lower())
    hi=arb(abs(x).upper())
    return interval(lo*lo,hi*hi)

def execute():
    ctx.prec=512
    p,h,zh,env=york.actual_parameters()
    N=80; cells=2048; sig=arb('4.005'); a=arb('.05'); b=arb('.1')
    R=arb('.2'); ln=arb(2).log()
    be=p['beta'];u=p['u'];ep=p['eps'];ta=p['tau']
    q=[arb(1),arb(0),u/be,(-4*be-2*u-ta)/be,
       (3*be+u+ep)/be,-2*ep/be,ep/be]
    dq=sum((abs(q[j])*R**j for j in range(1,7)),arb(0))
    dq1=sum((j*abs(q[j])*R**j for j in range(1,7)),arb(0))
    profile_guard=dq/(1-dq)+dq1/(4*(1-dq)**2)
    check('zero_free_disc_and_Cauchy_profile_guard',dq<1 and profile_guard<2*R**2)
    inv=[arb(1)]
    for n in range(1,N+1):
        inv.append(-sum((q[j]*inv[n-j] for j in range(1,min(6,n)+1)),arb(0)))
    isq=[sum((inv[j]*inv[n-j] for j in range(n+1)),arb(0)) for n in range(N+1)]
    bc={n:inv[n-3]+sum((j*q[j]*isq[n-j-3]/(n-j+1)
                       for j in range(1,min(6,n-3)+1)),arb(0))
        for n in range(3,N+1)}
    bc[3]-=1
    check('exact_removed_cubic_and_absent_quartic',bc[3].is_zero() and bc[4].is_zero())
    check('actual_reference_series_leading_terms',abs(bc[5])<arb('1e-9') and
          bc[6]>arb('.99') and bc[6]<arb('1.01'))
    KL=2*be*(p['Fl']+23*p['g']/(180*arb.pi()))/3
    pp=poly_mul(poly_mul([sig-4,arb(1)],[sig,arb(1)]),[sig-3,arb(1)])
    chi=[arb(1)]+[arb(0)]*21
    beta_norm=factorial(21)//factorial(10)**2
    for k in range(11):chi[11+k]=-arb(beta_norm*(-1)**k*comb(10,k))/(11+k)
    chider=[chi]
    for j in range(1,22):chider.append([i*v for i,v in enumerate(chider[-1])][1:])
    check('cutoff_endpoint_values',horner(chi,arb(0))==1 and horner(chi,arb(1)).is_zero())
    check('cutoff_C10_endpoint_jets',all(horner(chider[j],arb(0)).is_zero()
          and horner(chider[j],arb(1)).is_zero() for j in range(1,11)))
    polys=[[arb(0)]*5+[bc[n]*(arb(n)-sig)**k for n in range(5,N+1)] for k in range(7)]
    cs={n:KL*bc[n]*(n-4)*n*(n-3) for n in range(5,N+1)}
    inner=sum((cs[n]*cs[m]*a**(n+m-2*sig)*(1+(n-sig)*(m-sig))**3/(n+m-2*sig)
               for n in range(5,N+1) for m in range(5,N+1)),arb(0))
    check('inner_exact_integral_positive',inner>0)
    trans=[arb(0)]*4
    half=arb(1)/(2*cells)
    radius_powers=[half**k/factorial(k) for k in range(22)]
    for i in range(cells):
        mid=arb(2*i+1)/(2*cells)
        tt=mid+arb(0,half.abs_upper())
        zz=a*(ln*tt).exp()
        vals=[horner(c,mid) for c in chider]
        cut=[]
        for l in range(7):
            rad=sum((abs(vals[l+k])*radius_powers[k] for k in range(1,22-l)),arb(0))
            cut.append((vals[l]+arb(0,rad.abs_upper()))/ln**l)
        ek=[zz**(-sig)*horner(c,zz) for c in polys]
        for j in range(4):
            gj=KL*sum((pp[r]*sum((comb(j+r,l)*cut[l]*ek[j+r-l]
                                  for l in range(j+r+1)),arb(0)) for r in range(4)),arb(0))
            trans[j]+=sqrange(gj)*ln/cells
    normsq=inner+sum((comb(3,j)*trans[j] for j in range(4)),arb(0))
    check('whole_classical_norm_square_positive',normsq>0)
    classical=normsq.sqrt()
    B=[arb(1)]+[arb(factorial(21)//factorial(21-l)*2**(l-1))/ln**l for l in range(1,7)]
    ratios=[b/R*(1+1/(N+1+sig))**k for k in range(7)]
    check('entire_analytic_tail_ratios',all(v<1 for v in ratios))
    tail=[2*R**5*b**(-sig)*(N+1+sig)**k*(b/R)**(N+1)/
          ((2*(N+1-sig)).sqrt()*(1-ratios[k])) for k in range(7)]
    Aj=[]
    for j in range(4):
        Aj.append(abs(KL)*sum((abs(pp[r])*sum((comb(j+r,l)*B[l]*tail[j+r-l]
                                             for l in range(j+r+1)),arb(0)) for r in range(4)),arb(0)))
    Eseries=sum((comb(3,j)*Aj[j]**2 for j in range(4)),arb(0)).sqrt()
    previous=json.loads((ROOT/'evidence/current-coupled-global-return-attempt.json').read_bytes())['runs'][-1]['result']
    Qexact=dec(previous['global_cutoff_complete_normal_quantum_upper'])
    Qguard=arb('2.260e-10')
    check('inherited_complete_quantum_bound_applies_to_same_cutoff',Qexact<Qguard)
    check('complete_analytic_series_error_small',Eseries<arb('1e-8'))
    result=classical+arb(0,(Eseries+Qguard).abs_upper())
    check('full_conversion_lower_exceeds_2000',result>2000)
    check('full_conversion_upper_below_2400',result<2400)
    Dlo=arb(result.lower());Dhi=arb(result.upper());r=arb('1e-11')
    mass_needed=2*r/Dlo
    mass_split=r/(2*Dhi)
    m0=arb('4.2953e-8');Cg=arb('49.462')
    oldlower=dec(previous['global_cutoff_conversion_lower_per_mass'])
    check('conversion_improvement_over_previous_lower_exceeds_five_orders',Dlo/oldlower>arb('2e5'))
    check('both_graph_ceiling_necessary_mass_below_1e_minus14',mass_needed<arb('1e-14'))
    check('specified_positive_affine_envelope_fails',Dlo*m0>arb('8e-5'))
    check('specified_positive_graph_envelope_fails',1+Dlo*Cg>arb('90000'))
    check('explicit_split_conversion_sufficient',r/2+Dhi*mass_split<=r+arb('1e-150'))
    check('conditional_mass_column_debits_unchanged',(arb('5e-5')*m0).overlaps(arb('2.14765e-12')))
    check('conditional_fast_screen_preserved',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'N_reference_series':N,'transition_cells':cells,'precision_bits':ctx.prec,
       'parameters':{k:enc(p[k]) for k in ['beta','u','eps','tau','g','Fl']},
       'KL':enc(KL),'complex_disc_profile_guard':enc(profile_guard),
       'profile_coefficients':{str(n):enc(v) for n,v in bc.items()},
       'classical_inner_squared_H3':enc(inner),
       'classical_transition_derivative_L2':{str(j):enc(trans[j].sqrt()) for j in range(4)},
       'classical_total_H3':enc(classical),'analytic_series_H3_error':enc(Eseries),
       'complete_normal_quantum_H3_error':enc(Qguard),'full_conversion_H3':enc(result),
       'necessary_actual_mass_for_both_graph_ceilings':enc(mass_needed),
       'sufficient_actual_mass_with_half_lifted_graph_budget':enc(mass_split),
       'positive_estimator_affine_envelope_lower':enc(Dlo*m0),
       'positive_estimator_graph_slope_lower':enc(1+Dlo*Cg),
       'positive_estimator_affine_envelope_to_radius_lower':enc(Dlo*m0/r),
       'conditional_mass_allowance_to_necessary_mass_ratio':enc((m0+Cg*r)/mass_needed),
       'conversion_lower_improvement_factor':enc(Dlo/oldlower),
       'scope':{'complete_original_normal_graph_conversion':True,
                'actual_mass_computed':False,'actual_kappa_m_computed':False,
                'complete_physical_coupled_return':'unresolved',
                'common_regulated_core_computed':False,'nonlinear_admission':False,
                'construction_obstruction':False,'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    flint_root=Path(flint.__file__).resolve().parent
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
         'input_sha256_before':hashes,'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
         'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                    'flint_version':flint.__version__,'sympy_version':sympy.__version__,
                    'flint_file_sha256':{str(f.relative_to(flint_root)):sha256(f.read_bytes()).hexdigest()
                      for f in sorted(flint_root.rglob('*')) if f.is_file() and f.suffix in ['.py','.pyd','.dll']}}}
    try:
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('all_frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary={'mode':mode,'status':run['status'],'checks':len(CHECKS),
             'result':{k:v for k,v in run.get('result',{}).items() if k not in ['profile_coefficients','parameters']},
             'traceback':run.get('traceback')}
    message=json.dumps(summary,indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
