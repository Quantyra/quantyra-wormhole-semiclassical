"""Original-cutoff two-row positive estimator; not a coupled quantum solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-original-cutoff-coupled-estimator'
OWN=['notes/'+STEM+'.md','checks/current_original_cutoff_coupled_estimator.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-coupled-local-graph-bound.md',
 'checks/current_coupled_local_graph_bound.py',
 'evidence/current-coupled-local-graph-bound.json',
 'notes/current-global-mass-conversion.md',
 'evidence/current-global-mass-conversion.json',
 'notes/current-complete-weak-mass-residual.md',
 'evidence/current-complete-weak-mass-residual.json',
 'notes/current-actual-far-graph-energy-bound.md',
 'evidence/current-actual-far-graph-energy-bound.json',
 'notes/current-complete-strong-response-attempt.md',
 'notes/current-coupled-schur-graph-attempt.md',
 'notes/current-short-travel-strong-domain.md',
 'notes/current-strong-return-weight-repair.md']
CHECKS=[]
def check(name,ok):
    assert name not in {c['name'] for c in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(v):
    mid,rad,exp=v.mid_rad_10exp(50)
    return {'ball':str(v),'decimal_enclosure':{'mid':str(mid),'rad':str(rad),'exponent':str(exp)}}
def dec(v):
    d=v['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r['result']
def execute():
    ctx.prec=512
    previous=prior('current-coupled-local-graph-bound')
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    t=sp.Symbol('t')
    cut=1-sp.Rational(factorial(21),factorial(10)**2)*sum(
        ((-1)**k*comb(10,k)*t**(11+k)/sp.Integer(11+k) for k in range(11)))
    poly=sp.Poly(cut,t)
    check('original_cutoff_endpoints',cut.subs(t,0)==1 and cut.subs(t,1)==0)
    check('original_cutoff_endpoint_jets',all(sp.diff(cut,t,j).subs(t,0)==0 and sp.diff(cut,t,j).subs(t,1)==0 for j in range(1,11)))
    bern=[sum(poly.nth(j)*sp.Rational(comb(i,j),comb(21,j)) for j in range(i+1)) for i in range(22)]
    check('original_cutoff_bernstein_coefficients',bern==[1]*11+[0]*11)
    def rat(v):return arb(int(v.p))/int(v.q)
    c=[rat(poly.nth(k)) for k in range(22)]
    width=arb(2).log();cells=512;jet_cells=[]
    for cell in range(cells):
        mid=arb(2*cell+1)/(2*cells);rad=arb(1)/(2*cells)
        shifted=[sum((c[k]*comb(k,j)*mid**(k-j) for k in range(j,22)),arb(0)) for j in range(22)]
        jets=[]
        for j in range(4):
            cen=shifted[j]*factorial(j)
            err=sum((abs(shifted[k])*arb(factorial(k)//factorial(k-j))*rad**(k-j) for k in range(j+1,22)),arb(0))
            jets.append(cen+arb(0,err.abs_upper()))
        jet_cells.append(jets)
    cuts={};cut_rows={}
    for p in range(2,13):
        supremum=[arb(p)**j*arb('.5')**p for j in range(4)]
        integral=[supremum[j]**2/(2*p) for j in range(4)]
        for cell,jets in enumerate(jet_cells):
            exp=(p*width*(arb(cell+1)/cells-1)).exp()
            for j in range(4):
                rr=sum((comb(j,l)*arb(p)**(j-l)*jets[l]/width**l for l in range(j+1)),arb(0))
                bound=(exp*rr.abs_upper()).abs_upper()
                supremum[j]=max(supremum[j],bound)
                integral[j]+=width/cells*bound**2
        mm=sum((comb(3,j)*supremum[j] for j in range(4)),arb(0))
        ll=sum((comb(3,j)*integral[j].sqrt() for j in range(4)),arb(0))
        check('positive_whole_interval_cutoff_'+str(p),mm>0 and ll>0)
        cuts[p]=(mm,ll)
        cut_rows[str(p)]={'M_factoring_out_point_one_power':enc(mm),'L_factoring_out_point_one_power':enc(ll),
                         'derivative_supremum':[enc(v) for v in supremum],
                         'derivative_L2_squared':[enc(v) for v in integral]}
    terms={}
    for i in range(2):
        for j in range(2):
            entries=previous['matrix_bounds'][i][j]['terms']
            check('every_polynomial_power_retained_'+str(i)+str(j),{v['z_power'] for v in entries}==set(range(2,13)))
            terms[i,j]=[]
            for term in entries:
                p=term['z_power'];oldm,oldl=[dec(v) for v in previous['cutoff_constants'][str(p)]]
                assert oldm>0 and oldl>0
                lo=dec(term['low'])*10**p*cuts[p][1]/oldl
                hi=dec(term['high'])*10**p*cuts[p][0]/oldm
                terms[i,j].append((p,lo,hi))
    def matrix(removed):
        result=[];detail=[]
        for i in range(2):
            row=[];drow=[]
            for j in range(2):
                lo=sum((v[1] for v in terms[i,j] if v[0] not in removed),arb(0))
                hi=sum((v[2] for v in terms[i,j] if v[0] not in removed),arb(0))
                val=(lo**2+hi**2).sqrt();row.append(val)
                drow.append({'low':enc(lo),'high':enc(hi),'bound':enc(val),
                    'terms':[{'power':p,'low':enc(l),'high':enc(h)} for p,l,h in terms[i,j] if p not in removed]})
            result.append(row);detail.append(drow)
        return result,detail
    m,details=matrix(set());m3,d3=matrix({3});m34,d34=matrix({3,4})
    rho=(m[0][0]+m[1][1]+((m[0][0]-m[1][1])**2+4*m[0][1]*m[1][0]).sqrt())/2
    rows=[sum(row,arb(0)) for row in m]
    rem3=max(sum(row,arb(0)) for row in m3);rem34=max(sum(row,arb(0)) for row in m34)
    check('original_scalar_diagonal_majorant_exceeds_one',m[0][0]>1)
    check('original_shape_diagonal_majorant_exceeds_one',m[1][1]>1)
    check('full_majorant_spectral_radius_enclosure',rho>arb('360.4185') and rho<arb('360.4186'))
    check('original_max_row_enclosure',max(rows)>arb('365.9952') and max(rows)<arb('365.9953'))
    check('cubic_only_resummation_remainder_still_excessive',rem3>arb('19.4855') and rem3<arb('19.4856'))
    check('cubic_and_quartic_remainder_upper',rem34<arb('.002555496'))
    caps=[[arb('.0025554751'),arb('2.055945e-8')],[arb('4.696585e-8'),arb('1.046918e-7')]]
    check('all_four_residual_entries_enclosed',all(m34[i][j]<caps[i][j] for i in range(2) for j in range(2)))
    check('conditional_comparison_ceiling_rounded_down',arb('195.65673')*arb('.002555496')<arb('.5'))
    mass=prior('current-complete-weak-mass-residual')
    far=prior('current-actual-far-graph-energy-bound')
    check('actual_linear_mass_bound_retained',dec(mass['full_coupled_linear_mass_upper'])<arb('1.627532e-15'))
    check('actual_distant_bound_retained',dec(far['actual_linear_far_graph_upper'])<arb('1.24e-20'))
    mass_debit=2174*arb('1.627532e-15')
    shape_left=arb('1e-11')-mass_debit-arb('1.24e-20')
    check('shape_only_allowance_retained',shape_left>arb('6.4617454e-12'))
    k=arb('5e-5');affine=k*arb('4.2953e-8');feedback=k*arb('49.462')
    check('conditional_mass_affine_debit',(affine-arb('2.14765e-12')).contains(0))
    check('conditional_mass_feedback_debit',(feedback-arb('.0024731')).contains(0))
    energy_threshold=arb('.0012069')/(arb('.0195239')+arb('4.583e-13'))
    check('conditional_energy_threshold',energy_threshold<arb('.061816543'))
    screen=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('retained_fast_screen_arithmetic',screen<arb('.1599164') and (screen-arb('.135973')).contains(0))
    return {'precision_bits':ctx.prec,'cutoff_cells':cells,'cutoff_norms':cut_rows,
      'original_matrix':details,'spectral_radius_of_upper_bound_matrix':enc(rho),
      'row_sums_of_upper_bound_matrix':[enc(v) for v in rows],
      'remove_cubic_only_matrix':d3,'remove_cubic_only_max_row':enc(rem3),
      'remove_cubic_and_quartic_matrix':d34,'remove_cubic_and_quartic_max_row':enc(rem34),
      'conditional_fast_comparison_ceiling':enc(arb('.5')/arb('.002555496')),
      'original_mass_shape_conversion_debit':enc(mass_debit),'shape_only_remaining_allowance':enc(shape_left),
      'conditional_mass_column_affine_debit':enc(affine),'conditional_mass_column_feedback_debit':enc(feedback),
      'conditional_energy_comparison_ceiling':enc(energy_threshold),'conditional_two_fast_screen':enc(screen),
      'scope':{'specified_positive_estimator_failed':True,'lower_bound_on_actual_operator':False,
       'full_near_quantum_sums_evaluated':False,'quantum_contacts_assigned_zero':False,
       'finite_geometry_inverse_certified':False,'actual_coupled_graph_certified':False,
       'actual_kappa_m_computed':False,'independent_phase_columns_bounded':False,
       'nonlinear_admission':False,'nonlinear_contraction':False,'exact_compatibility':False,
       'physical_matching':False,'validity_noise':False,'all_sector_stability':False,
       'actual_controlled_solution':False,'construction_obstruction':False,'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
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
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'spectral_radius':run.get('result',{}).get('spectral_radius_of_upper_bound_matrix'),
       'remainder_max_row':run.get('result',{}).get('remove_cubic_and_quartic_max_row'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
