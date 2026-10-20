"""S150: averaged canonical pressure matching on a joint chart contract.

Manufactured Ward data test the map, not a physical quantum state.
The global coordinate/source contract is an unresolved hypothesis.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
STEM="averaged-ward-matching"
INPUTS=[
    "checks/averaged_ward_matching.py","notes/averaged-ward-matching.md",
    "notes/interior-pressure-and-bolt-primitive.md",
    "checks/interior_pressure_and_bolt_primitive.py",
    "notes/general-interior-dirac-response.md",
    "notes/outer-finite-background-response.md",
    "notes/outer-ward-and-constraint-transfer.md",
    "notes/spatial-response-norm-audit.md",
]
ROWS,NUMERICS=[],[]
PROTOCOL={"quadrature_nodes":[32,64],"digits":55,"steps":[".001",".0005"],
          "mixed_tolerance":"2e-6","identity_tolerance":"1e-10","refinement_ratio":".35"}


def check(name,passed,detail=""):
    assert not any(row["name"]==name for row in ROWS)
    row={"name":name,"passed":bool(passed),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def exact_checks():
    z,Y=s.symbols("z Y")
    poly=140*z**3*(1-z)**3
    Z=35*z**4-84*z**5+70*z**6-20*z**7
    check("overlap_density_normalization",s.integrate(poly,(z,0,1))==1)
    check("overlap_cumulative_polynomial",s.diff(Z,z)==s.expand(poly) and Z.subs(z,1)==1)
    check("overlap_c2_zero_extension",all(s.diff(poly,z,j).subs(z,endpoint)==0 for j in range(3) for endpoint in [0,1]))
    check("overlap_first_derivative_factor",s.simplify(s.diff(poly,z)-420*z*z*(1-z)**2*(1-2*z))==0)
    check("overlap_second_derivative_lobes",s.simplify(s.diff(poly,z,2)-840*z*(1-z)*(1-5*z+5*z*z))==0)
    check("overlap_peak_value",poly.subs(z,s.Rational(1,2))==s.Rational(35,16))
    # (Z a^2 P)'=zeta a^2 P+Z a a' S; integrate across the overlap.
    aa,ap,PP,SS,zz,zp=s.symbols("a ap P S Z zp")
    ward=ap/aa*(SS-2*PP)
    check("averaged_ward_recovery",s.simplify(zp*aa**2*PP+zz*(2*aa*ap*PP+aa**2*ward)-
                                            (zp*aa**2*PP+zz*aa*ap*SS))==0)
    ss,tt,h,k,yh,yk,yhk=s.symbols("ss tt h k yh yk yhk")
    zet=s.Function("zeta")
    field=s.exp(3*(ss*h+tt*k))*zet(Y+ss*yh+tt*yk+ss*tt*yhk)
    zero={ss:0,tt:0}
    first=3*h*zet(Y)+s.diff(zet(Y),Y)*yh
    mixed=9*h*k*zet(Y)+3*s.diff(zet(Y),Y)*(h*yk+k*yh)+s.diff(zet(Y),Y,2)*yh*yk+s.diff(zet(Y),Y)*yhk
    check("fixed_disk_overlap_first_kernel",s.simplify(s.diff(field,ss).subs(zero)-first)==0)
    # Evaluate derivative substitutions before comparing their dummy variables.
    check("fixed_disk_overlap_mixed_kernel",s.simplify((s.diff(field,ss,tt).subs(zero)-mixed).doit())==0)
    j,hp,kp,jp,w=s.symbols("j hp kp jp w")
    qa=s.exp(-2*(ss*h+tt*k+ss*tt*j))*(w-ss*hp-tt*kp-ss*tt*jp)
    check("proper_ward_factor_first",s.simplify(s.diff(qa,ss).subs(zero)+2*h*w+hp)==0)
    check("proper_ward_factor_mixed_acceleration",s.simplify(s.diff(qa,ss,tt).subs(zero)-
                (4*h*k*w+2*h*kp+2*k*hp-2*j*w-jp))==0)
    p0,ph,pk,phk=s.symbols("p0 ph pk phk")
    ff=s.exp(-2*(ss*h+tt*k+ss*tt*j))*(p0+ss*ph+tt*pk+ss*tt*phk)
    f_hk=s.diff(ff,ss,tt).subs(zero)
    check("pressure_recovery_second_response",s.simplify(
          f_hk+2*h*pk+2*k*ph+(2*j-4*h*k)*p0-phk)==0)
    G,B,C,x=s.symbols("G B C x",positive=True)
    scale=B*B*C/(2*G)
    check("normalized_mass_memory_prefactor",s.simplify(G/(B*B*x)*scale-C/(2*x))==0)
    check("normalized_lapse_memory_prefactor",s.simplify(G/(C*B*B)*scale-s.Rational(1,2))==0)
    eps=s.symbols("eps",positive=True)
    # Local Ward fixture W=1+x, trace x_+^(1/2), P=x_+^(3/2)/(3(1+x)).
    pp=eps**s.Rational(3,2)/(3*(1+eps))
    check("rough_trace_local_ward",s.simplify(s.diff((1+eps)*pp,eps)-s.sqrt(eps)/2)==0)
    mixedquot=(2*eps)**s.Rational(3,2)/(12*eps**2*(1+2*eps))
    check("rough_trace_moving_cut_hessian_diverges",s.limit(mixedquot,eps,0,dir="+")==s.oo)


def budget_checks():
    mp.mp.dps=60
    ya,yb,y0=Q(41,1000),Q(49,1000),Q(1,20)
    width=yb-ya
    tv=Q(35,8)/width
    tv2=Q(105)/width**2
    check("support_inside_interior_cover",mp.mpf(".041")*mp.exp(-mp.mpf(".001"))>mp.mpf(".04"))
    check("density_variation_constants",tv==Q(4375,8) and tv2==1640625)
    K0=Q(103,100)
    K1=K0*(3+yb*tv)
    K2=K0*(9+7*yb*tv+yb*yb*tv2)
    check("overlap_kernel_response_budgets",K1<31 and K2<4300,{"first":str(K1),"mixed":str(K2)})
    p0,p1,p2=Q(6,10**10),Q(36,10**8),Q(44,10**7)
    q0=K0*p0
    q1=K0*p1+31*p0
    q2=K0*p2+62*p1+4300*p0
    check("overlap_pressure_moment_budgets",q0==Q(618,10**12) and q1<Q(4,10**7) and q2<Q(3,10**5))
    eta=Q(85,10**34)
    cubic0=2*eta/(15*ya)
    check("full_outer_trace_value_budget",2*eta*4*10**25+cubic0<Q(7,10**7))
    check("full_outer_trace_first_budget",2*eta*2*10**25+cubic0*Q(75,1000)<Q(35,10**8))
    check("full_outer_trace_mixed_acceleration_budget",2*eta*(10**30+2*10**25)+cubic0*Q(75,1000)<Q(18,1000))
    check("proper_ward_factor_budgets",K0*Q(1,5)==Q(206,1000) and K0*140<145 and K0*(48000+140)<50000)
    Hbound=mp.mpf(".009")*(1+mp.log(mp.mpf(".041"))**2)
    check("overlap_log_weight_integral_budget",Hbound<mp.mpf(".102"),mp.nstr(Hbound,25))
    imix=Q(102,1000)*(Q(206,1000)*Q(18,1000)+290*Q(35,10**8)+50000*Q(7,10**7))
    check("outer_correction_mixed_budget",imix<Q(4,1000),str(imix))
    pmix=K0*(Q(3,10**5)+Q(4,1000))+400*Q(37,10**8)+40200*p0
    check("matching_pressure_mixed_budget",pmix<Q(44,10000),str(pmix))
    memM=Q(25,8)*(Q(5,100)*Q(44,10000)+2*Q(51,100)*Q(37,10**8)+Q(611,100)*p0)
    memA=Q(1,2)*(Q(1,10**5)*Q(44,10000)+2*Q(14,10**5)*Q(37,10**8)+Q(214,10**5)*p0)
    check("matching_mass_memory_budget",memM<Q(7,10000),str(memM))
    check("matching_lapse_memory_budget",memA<Q(23,10**9),str(memA))
    # First response estimates of S147's regular trace-forcing component.
    r0,r1=Q(4*10**25),Q(2*10**25)
    pp0,pp1=Q(284,100)*10**24,Q(400142,100000)*10**27
    firstM=2*Q(42,100)*Q(5,10**32)*(y0*(pp1+pp0)+Q(899,1000)*(r1+r0))
    firstA=16*Q(42,100)**2*eta*(Q(1,640000)*(2*pp1+8*pp0)+Q(1813,10**8)*(r1+4*r0))
    check("regular_mass_first_for_acceleration",firstM<Q(11,10**6),str(firstM))
    check("regular_lapse_first_for_acceleration",firstA<Q(4,10**10),str(firstA))
    check("combined_outer_component_budgets",
          Q(39,1000)+Q(11,10**6)+Q(7,10000)<Q(4,100) and
          Q(46,10**8)+Q(4,10**10)+Q(23,10**9)<Q(49,10**8))
    check("conditional_component_taylor_budgets",
          Q(4,100)/2*Q(1,10**14)==Q(2,10**16) and
          Q(49,10**8)/2*Q(1,10**14)==Q(245,10**23))
    NUMERICS.append({"moment_mixed_upper":str(q2),"correction_mixed_upper":str(imix),
                    "pressure_mixed_upper":str(pmix),"mass_memory_upper":str(memM),"lapse_memory_upper":str(memA),
                    "joint_global_contract_proved":False})


def unary(x,f,fp,fpp):
    if not isinstance(x,Jet):
        return f(x)
    a,b,c,d=x.v
    return Jet(f(a),fp(a)*b,fp(a)*c,fpp(a)*b*c+fp(a)*d)


def ex(x):
    return x.exp() if isinstance(x,Jet) else mp.exp(x)


def sh(x):
    return unary(x,mp.sinh,mp.cosh,mp.sinh)


def ch(x):
    return unary(x,mp.cosh,mp.sinh,mp.cosh)


def shape(x):
    x2=x*x
    return mp.mpf(".03")+x2/100, -mp.mpf(".02")+mp.mpf(".002")*x2*x2


def fields(x,ss,tt):
    h,k=shape(x)
    b=ss*h+tt*k
    bp=ss*x/50+tt*mp.mpf(".008")*x*x*x
    trace=1+x*x/5+ss*(mp.mpf(".3")+x*x/10)+tt*(-mp.mpf(".2")+mp.mpf(".07")*x*x*x*x)+ss*tt*(mp.mpf(".11")+mp.mpf(".03")*x*x)
    return b,bp,trace


def numerical_checks():
    mp.mp.dps=55
    ya,yb,y0,D=map(mp.mpf,[".041",".049",".05","2"])
    width=yb-ya
    snapshots=[]
    for n in PROTOCOL["quadrature_nodes"]:
        rx,rw=leggauss(n)
        nodes=[(mp.mpf(float(x))+1)/2 for x in rx]
        weights=[mp.mpf(float(w))/2 for w in rw]
        def quad(fun,lo,hi):
            return sum((hi-lo)*wt*fun(lo+(hi-lo)*z) for z,wt in zip(nodes,weights))
        # Small, separately chosen proper-map quadrature for the inverse solve.
        cx,cw=leggauss(16)
        cn=[(mp.mpf(float(x))+1)/2 for x in cx]
        cwt=[mp.mpf(float(w))/2 for w in cw]
        def properY(x,ss,tt):
            return sum((D-x)*wt*mp.exp(fields(x+(D-x)*z,ss,tt)[0]) for z,wt in zip(cn,cwt))
        def cut(Y,ss=0,tt=0):
            return mp.findroot(lambda x:properY(x,ss,tt)-Y,(D-Y-mp.mpf(".001"),D-Y+mp.mpf(".001")))
        def cutjet(Y,noacc=False):
            x=D-Y
            yh=quad(lambda v:shape(v)[0],x,D)
            yk=quad(lambda v:shape(v)[1],x,D)
            yhk=quad(lambda v:shape(v)[0]*shape(v)[1],x,D)
            h,k=shape(x)
            return Jet(x,yh,yk,0 if noacc else yhk-h*yk-k*yh)
        def pressure(c,ss,tt):
            bc,_,_=fields(c,ss,tt)
            denom=ex(2*bc)*sh(c)*sh(c)
            def kernel(z):
                x=c*z
                b,bp,T=fields(x,ss,tt)
                return c*ex(2*b)*sh(x)*(ch(x)+bp*sh(x))*T/denom
            return sum(wt*kernel(z) for z,wt in zip(nodes,weights))
        def zeta_jet(Y):
            z=(Y.v[0]-ya)/width
            v=140*z**3*(1-z)**3/width
            vp=420*z*z*(1-z)**2*(1-2*z)/width**2
            vpp=840*z*(1-z)*(1-5*z+5*z*z)/width**3
            return Jet(v,vp*Y.v[1],vp*Y.v[2],vpp*Y.v[1]*Y.v[2]+vp*Y.v[3])
        def zeta(Y):
            z=(Y-ya)/width
            return 140*z**3*(1-z)**3/width
        def cumulative(Y):
            if Y>=yb:
                return mp.mpf(1)
            z=(Y-ya)/width
            return 35*z**4-84*z**5+70*z**6-20*z**7
        sj,tj=Jet(0,1,0),Jet(0,0,1)
        def momentjet(freeze_transport=False):
            def integrand(x):
                yh=quad(lambda v:shape(v)[0],x,D)
                yk=quad(lambda v:shape(v)[1],x,D)
                yhk=quad(lambda v:shape(v)[0]*shape(v)[1],x,D)
                Yj=Jet(D-x) if freeze_transport else Jet(D-x,yh,yk,yhk)
                bj,_,_=fields(x,sj,tj)
                weight=zeta_jet(Yj)*ex(3*bj)*mp.sinh(x)**2/mp.sinh(D)**2
                return weight*pressure(x,sj,tj)
            return quad(integrand,D-yb,D-ya)
        def correctionjet(noacc=False):
            def integrand(Y):
                x=cutjet(Y,noacc)
                b,bp,T=fields(x,sj,tj)
                a=ex(b)*sh(x)/mp.sinh(D)
                ay=-(bp*sh(x)+ch(x))/mp.sinh(D)
                return cumulative(Y)*a*ay*T
            return quad(integrand,ya,yb)+quad(integrand,yb,y0)
        Qj,Ij=momentjet(),correctionjet()
        c0=cutjet(y0)
        b0,_,_=fields(c0,sj,tj)
        a0=ex(b0)*sh(c0)/mp.sinh(D)
        recovered=(Qj+Ij)/(a0*a0)
        direct=pressure(c0,sj,tj)
        tol=mp.mpf(PROTOCOL["identity_tolerance"])
        check(f"averaged_matching_all_response_components_{n}",
              all(abs(a-b)<tol*(1+abs(b)) for a,b in zip(recovered.v,direct.v)),
              {"recovered":[mp.nstr(v,25) for v in recovered.v],"direct":[mp.nstr(v,25) for v in direct.v]})
        check(f"overlap_transport_omission_detected_{n}",abs(Qj.v[3]-momentjet(True).v[3])>mp.mpf("1e-8"))
        check(f"outer_correction_omission_detected_{n}",abs(Ij.v[0])>mp.mpf("1e-5"))
        check(f"proper_coordinate_acceleration_omission_detected_{n}",
              abs(Ij.v[3]-correctionjet(True).v[3])>mp.mpf("1e-10"))
        def scalar_maps(ss,tt):
            def moment(Y):
                x=cut(Y,ss,tt)
                b,_,_=fields(x,ss,tt)
                a=mp.exp(b)*mp.sinh(x)/mp.sinh(D)
                return zeta(Y)*a*a*pressure(x,ss,tt)
            def correction(Y):
                x=cut(Y,ss,tt)
                b,bp,T=fields(x,ss,tt)
                a=mp.exp(b)*mp.sinh(x)/mp.sinh(D)
                ay=-(bp*mp.sinh(x)+mp.cosh(x))/mp.sinh(D)
                return cumulative(Y)*a*ay*T
            qv=quad(moment,ya,yb)
            iv=quad(correction,ya,yb)+quad(correction,yb,y0)
            x=cut(y0,ss,tt)
            b,_,_=fields(x,ss,tt)
            a=mp.exp(b)*mp.sinh(x)/mp.sinh(D)
            return qv,iv,(qv+iv)/(a*a),pressure(x,ss,tt)
        errors=[[],[],[],[]]
        targets=[Qj.v[3],Ij.v[3],recovered.v[3],direct.v[3]]
        for step_text in PROTOCOL["steps"]:
            step=mp.mpf(step_text)
            vals=[scalar_maps(ss,tt) for ss,tt in [(step,step),(step,-step),(-step,step),(-step,-step)]]
            for j in range(4):
                observed=(vals[0][j]-vals[1][j]-vals[2][j]+vals[3][j])/(4*step**2)
                errors[j].append(abs(observed-targets[j]))
        for label,errs,target in zip(["overlap_moment","outer_correction","recovered_pressure","direct_pressure"],errors,targets):
            check(f"{label}_mixed_difference_{n}",max(errs)<mp.mpf(PROTOCOL["mixed_tolerance"])*(1+abs(target)),
                  {"target":mp.nstr(target,25),"errors":[mp.nstr(e,16) for e in errs]})
            check(f"{label}_step_refinement_{n}",errs[1]<mp.mpf(PROTOCOL["refinement_ratio"])*errs[0] or max(errs)<mp.mpf("1e-35"))
        # Independently adaptive base pressure integral.
        c=D-y0
        adaptive=mp.quad(lambda x:mp.sinh(x)*mp.cosh(x)*(1+x*x/5),[0,c/2,c])/mp.sinh(c)**2
        check(f"independent_adaptive_base_pressure_{n}",abs(adaptive-direct.v[0])<tol)
        snapshots.append(targets)
        NUMERICS.append({"nodes":n,"mixed_targets":[mp.nstr(v,30) for v in targets],
                         "manufactured_ward_data_are_physical_states":False})
    check("quadrature_refinement_all_mixed_maps",
          all(abs(a-b)<mp.mpf(PROTOCOL["identity_tolerance"])*(1+abs(b)) for a,b in zip(*snapshots)))
    quotients=[mp.sqrt(2)/(6*mp.sqrt(eps)*(1+2*eps)) for eps in map(mp.mpf,[".001",".0001",".00001"])]
    check("missing_proper_source_contract_detected",all(b>3*a for a,b in zip(quotients,quotients[1:])),
          [mp.nstr(v,20) for v in quotients])


def hashes():
    return {name:sha256((ROOT/name).read_bytes()).hexdigest() for name in INPUTS}


def main():
    global Jet
    before=hashes()
    error=None
    try:
        spec=importlib.util.spec_from_file_location("s150_pressure_input",ROOT/"checks"/"interior_pressure_and_bolt_primitive.py")
        dependency=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(dependency)
        Jet=dependency.Jet
        exact_checks()
        budget_checks()
        numerical_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    after=hashes()
    ok=error is None and before==after and all(row["passed"] for row in ROWS)
    record={"story":"S150","created_utc":datetime.now(timezone.utc).isoformat(),"status":"passed" if ok else "failed",
            "python":sys.version,"numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
            "protocol":PROTOCOL,"input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
            "checks":ROWS,"numerics":NUMERICS,"exception":error,
            "scope":"averaged HH matching on a stated joint coordinate/source contract",
            "unresolved":["global preservation of the joint contract","cold weighted response","other source sectors",
                          "actual global inverse and exterior/state solve","all physical completion gates"]}
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"]=(ROOT/INPUTS[1]).read_text(encoding="utf-8")
    output=ROOT/"evidence"/f"{STEM}.json"
    history=json.loads(output.read_text(encoding="utf-8")) if output.exists() else {"runs":[]}
    history["runs"].append(record)
    output.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),"unchanged_inputs":before==after}),flush=True)
    return 0 if ok else 1


if __name__=="__main__":
    raise SystemExit(main())
