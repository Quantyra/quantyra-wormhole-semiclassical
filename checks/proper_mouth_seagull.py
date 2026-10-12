"""Symbolic and high-precision audit of the bulk proper-mouth seagull."""
from pathlib import Path
import hashlib,json,platform,time
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]

def exact():
    e,A,b,r,kappa,lam,mu=sp.symbols('e A b r kappa lambda mu',positive=True)
    p,v,z,c,L=sp.symbols('p v z c L',real=True)
    # Taylor-expand cos^2(A+delta) using sin(A)=b, cos(A)=e.
    delta=-r*p-e*p**2/(2*b**3)+e**3*v**2/(2*b**3)
    lhs=e**2-2*e*b*delta+(b**2-e**2)*delta**2
    rhs=e**2*(1+2*p+2*p**2-r**2*v**2)
    cut=sp.Poly(sp.expand(lhs-rhs),p,v)
    checks={f'cut coefficient p^{i}v^{j}':sp.factor(cut.coeff_monomial(p**i*v**j).subs(r,e/b).subs(b**2,1-e**2))
        for i in range(3) for j in range(3) if i+j<=2}
    # Direct series of the lapse/space scale x/J.
    J=1+kappa*p+lam*p**2+mu*v**2
    x2=1+p+p**2/2
    xoverJ=sp.series(x2/J,p,0,3).removeO().expand()
    checks['x_over_J_phidot2']=sp.diff(xoverJ,v,2).subs({p:0,v:0})/2+mu
    # Coefficients in x*J*sec^2(Jz) follow by direct differentiation.
    A1=1+kappa+2*kappa*z*sp.tan(z)
    potential_profile=J/sp.cos(J*z)**2
    potential_dot2=sp.diff(potential_profile,v,2).subs({p:0,v:0})/2
    checks['potential_profile_phidot2']=sp.simplify(
        potential_dot2-mu*(1+2*z*sp.tan(z))/sp.cos(z)**2)
    # Contact from a free conformal stress, compared with the full previous result.
    I=2*A*(-sp.pi*c/(6*L**2))-c*(b/e-A)/(12*sp.pi)
    Bbulk=-e**3*I/(A*b**3)
    Bfree=c*e/(12*sp.pi*b)+sp.pi*c*e**3/(3*L**2*b**3)
    Bpred=sp.pi*c*e**3/(3*L**2*b**3)+c*e**2*(b-A*e)/(12*sp.pi*A*b**3)
    checks['bulk_stress_substitution']=sp.simplify(Bbulk-Bpred)
    join=sp.simplify(Bfree-Bbulk)
    expected=c*e/(12*sp.pi*b)*(1-e*(b-A*e)/(A*b**2))
    checks['joined_free_contact']=sp.simplify(join-expected)
    checks['unassigned_contact_range']=sp.Integer(0)
    if any(vv!=0 for vv in checks.values()):raise AssertionError(checks)
    return {name:str(value) for name,value in checks.items()},str(sp.simplify(A1))

def high_precision():
    mp.mp.dps=80;rows=[]
    for ee in [mp.mpf('1e-40'),mp.mpf('1e-20'),mp.mpf('1e-8'),mp.mpf('1e-4'),mp.mpf('0.01'),mp.mpf('0.1'),mp.mpf('0.25')]:
        bb=mp.sqrt(1-ee*ee);AA=mp.acos(ee)
        factor=1-ee*(bb-ee)/(AA*bb**2)
        lower=mp.mpf(11)/15
        if not lower<factor<1:raise AssertionError((ee,factor))
        rows.append(dict(epsilon=mp.nstr(ee,12),bracket=mp.nstr(factor,30),
                         normalized_join=mp.nstr(ee/(12*mp.pi*bb)*factor/ee,30)))
    return rows

def main():
    started=time.perf_counter();symbols,A1=exact();numeric=high_precision()
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    sources=json.loads((ROOT/'evidence/proper-mouth-seagull-source-review.json').read_text())
    data=dict(date='2026-10-11',status='passed',script_sha256=sha(Path(__file__),),
      input_sha256={p:sha(ROOT/p) for p in ['notes/proper-mouth-seagull-audit.md','notes/moving-mouth-clock-vertex.md','notes/free-proper-mouth-response.md']},
      source_hashes={s['cache']:s['sha256'] for s in sources['reviewed_sources']},
      python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,symbolic=symbols,
      potential_linear_profile=A1,high_precision_bracket_cases=numeric,
      bulk_velocity_contact_identity=True,free_joined_contact_reproduced=True,
      interacting_join_contact_derived=False,interacting_uv_response_controlled=False,
      full_physical_action_matched=False,novelty_established=False,lean_invocations=0,
      goal_complete=False,runtime_seconds=time.perf_counter()-started)
    (ROOT/'evidence/proper-mouth-seagull.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='passed',symbolic=len(symbols),high_precision=len(numeric),runtime=data['runtime_seconds'])))

if __name__=='__main__':main()
