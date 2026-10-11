"""Checks of the Gaussian singlet's full subtracted stress comparison.

Analytic arguments establish uniformity. Finite matrices and meshes do not
certify a four-dimensional state, static renormalization or physical damping.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
import scipy
from scipy.integrate import quad
from scipy.linalg import eigh, null_space
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('graph',ROOT/'checks/schwinger_return_feedback.py')
graph=importlib.util.module_from_spec(spec);spec.loader.exec_module(graph)
PI=np.pi
BP=Poly([0,0,0,0,0,126,-420,540,-315,70])


def beta(t,order=0):
    t=np.asarray(t)
    middle=(t>.5)&(t<1)
    return np.where(middle,2**order*BP.deriv(order)(np.clip(2*t-1,0,1)),
                    np.where(t>=1,1.,0.) if order==0 else 0.)


def theta(t,order=0):
    return beta(t,order)-beta(np.asarray(t)/2,order)/2**order


S,R,W=sy.symbols('s r w')
BASE=R*W**2/(S*(W**2-R**2*S**2))
DERIV=[sy.lambdify((S,R,W),sy.diff(BASE,S,n),'numpy') for n in range(5)]


def g_derivatives(s,r,w):
    raw=[fn(s,r,w) for fn in DERIV]
    return np.array([sum(comb(k,j)*raw[j]*theta(s,k-j) for j in range(k+1)) for k in range(5)])


def s_quadrature(r,w,n):
    cuts=[.5,1.,2.]
    center=w.real/r
    for offset in [-16,-4,-1,0,1,4,16]:
        t=center+offset*w.imag/r
        if .5<t<2:cuts.append(t)
    cuts=sorted(set(cuts));nodes,weights=leggauss(n)
    xs=[];ws=[]
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        xs.extend((lo+hi)/2+(hi-lo)*nodes/2);ws.extend((hi-lo)*weights/2)
    return np.array(xs),np.array(ws)


def one_dimensional_norms(r,w,n=80):
    s,weights=s_quadrature(r,w,n)
    return abs(g_derivatives(s,r,w))@weights


def f_derivatives(a,eps):
    t=a/eps
    return [a*beta(t),beta(t)+t*beta(t,1),
            (2*beta(t,1)+t*beta(t,2))/eps]


def annular_norms(r,w,eps,n=40):
    ss,sw=s_quadrature(r,w,n)
    gd=g_derivatives(ss,r,w)
    nodes,weights=leggauss(n)
    result=np.zeros((3,3))
    for index,(s,weight) in enumerate(zip(ss,sw)):
        cuts=sorted(set([0.,s]+[v for v in [eps/2,eps,s-eps,s-eps/2] if 0<v<s]))
        aa=[];aw=[]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            aa.extend((lo+hi)/2+(hi-lo)*nodes/2);aw.extend((hi-lo)*weights/2)
        aa=np.array(aa);aw=np.array(aw);bb=s-aa
        fa,fb=f_derivatives(aa,eps),f_derivatives(bb,eps)
        for i in range(3):
            for j in range(3):
                value=sum(comb(i,k)*comb(j,l)*fa[i-k]*fb[j-l]*gd[k+l,index]
                          for k in range(i+1) for l in range(j+1))
                result[i,j]+=weight*np.dot(aw,abs(value))
    return result


def pair_kernel(a,b,w):
    total=a+b
    return w*w*a*b/(total*(w*w-total*total))


def exact_checks():
    a,b,w,r=sy.symbols('a b w r',positive=True)
    u=a*b/2
    kernel=w*w*a*b/((a+b)*(w*w-(a+b)**2))
    values={
      'pair_weight_and_subtraction':sy.simplify(kernel-2*w*w*u/((a+b)*(w*w-(a+b)**2))),
      'annular_frequency_rescaling':sy.simplify(kernel.subs({a:r*a,b:r*b})-r*w*w*a*b/((a+b)*(w*w-r*r*(a+b)**2))),
      'imaginary_frequency_positivity':sy.simplify(kernel.subs(w,sy.I)-a*b/((a+b)*(1+(a+b)**2))),
      'resonance_partial_fraction':sy.simplify(BASE-(R/S+R**2/2*(1/(W-R*S)-1/(W+R*S))))}
    assert all(v==0 for v in values.values())
    ts=np.geomspace(1e-5,1e5,151)
    partition=sum(theta(ts/(2.**k)) for k in range(-30,31))
    assert np.max(abs(partition-1))<2e-12
    return {k:str(v) for k,v in values.items()}|{'dyadic_partition_max_error':float(np.max(abs(partition-1)))}


def axis_checks():
    rows=[]
    reference=quad(lambda t:abs(2*beta(t,1)+t*beta(t,2)),.5,1,
                   points=[.75],epsabs=1e-10,epsrel=1e-10)[0]
    for eps in [.2,.02,.002,.0002]:
        value=quad(lambda a:abs(f_derivatives(np.array(a),eps)[2]),eps/2,eps,
                   points=[.75*eps],epsabs=1e-10,epsrel=1e-10)[0]
        assert abs(value-reference)<2e-9
        rows.append(dict(epsilon=eps,second_derivative_l1=value,scale_invariant_reference=reference))
    return rows


def annular_checks():
    rows=[]
    for r,w in [(1.,1.+.08j),(.03125,1.+.2j),(16.,1.+.2j)]:
        group=[]
        for eps in [.1,.01,.001]:
            coarse=annular_norms(r,w,eps,32)
            fine=annular_norms(r,w,eps,64)
            change=np.max(abs(fine-coarse))/(1+np.max(fine))
            assert change<.006
            row=dict(r=r,w=[w.real,w.imag],axis_epsilon=eps,l1_norms=fine.tolist(),
                     refinement=float(change),sum_norms=float(np.sum(fine)),
                     one_dimensional_derivative_norms=one_dimensional_norms(r,w).tolist())
            group.append(row);rows.append(row)
        assert group[-1]['sum_norms']<1.3*group[0]['sum_norms']
    tails=[]
    for r in [.015625,.03125,.0625,8.,16.,32.]:
        norms=one_dimensional_norms(r,1.+.2j)
        scaled=np.sum(norms)/(r if r<1 else 1/r)
        tails.append(dict(r=r,norm_sum=float(np.sum(norms)),normalized_tail_coefficient=float(scaled)))
    assert max(x['normalized_tail_coefficient'] for x in tails)<2*min(x['normalized_tail_coefficient'] for x in tails)
    resonance=[]
    for eta in [.2,.1,.05,.025]:
        n1=one_dimensional_norms(1.,1.+1j*eta,48)
        n2=one_dimensional_norms(1.,1.+1j*eta,96)
        assert np.max(abs(n1-n2))/(1+max(n2))<.001
        resonance.append(dict(eta=eta,derivative_norms=n2.tolist(),scaled_fourth_derivative=float(eta**4*n2[4])))
    return rows,tails,resonance


def finite_rank_checks():
    rng=np.random.default_rng(531011)
    H=np.diag(np.r_[-np.geomspace(.08,20,12),np.geomspace(.11,18,13)])
    vectors=rng.normal(size=(25,2))*.2
    J=H+vectors@vectors.T
    orth,_=np.linalg.qr(rng.normal(size=(25,25)));X=orth[:,:9]@orth[:,:9].T
    def statistic(matrix,w,r=1):
        ev,U=eigh(matrix/r);overlap=U.T@X@U
        left=ev<0;right=ev>0
        a=-ev[left,None];b=ev[None,right]
        return np.sum(pair_kernel(r*a,r*b,w)*abs(overlap[np.ix_(left,right)])**2)
    rows=[]
    for w in [.9+.3j,1.1+.08j]:
        first=statistic(H,w);second=statistic(J,w)
        for r in [.0625,1,16]:
            error=max(abs(statistic(H,w,r)-first),abs(statistic(J,w,r)-second))
            assert error<1e-10
            diff=np.linalg.inv(H/r-w*np.eye(25))-np.linalg.inv(J/r-w*np.eye(25))
            sv=np.linalg.svd(diff,compute_uv=False)
            assert np.linalg.norm(sv[2:])<1e-10 and sum(sv)<=4/w.imag+1e-10
            rows.append(dict(w=[w.real,w.imag],frequency_unit=r,statistic_difference=float(abs(first-second)),
                             rescaling_error=float(error),resolvent_rank_tail=float(np.linalg.norm(sv[2:])),
                             resolvent_trace_norm=float(sum(sv)),trace_bound=4/w.imag))
    return rows


def mass_matrix(T,q,n,mass,variable):
    if not variable:
        mt=graph.assemble_circles(T,np.ones(q),n)[5][:2*n,:2*n]
    else:
        mt=np.zeros((2*n,2*n));nodes,weights=leggauss(5)
        for edge in range(n):
            s=(edge+.5)*T/n+nodes*T/(2*n)
            profile=(1+.2*np.cos(2*PI*s/T))**2
            shape=np.array([(1-nodes)/2,(1+nodes)/2])
            elem=(shape*(weights*profile*T/(2*n)))@shape.T
            ids=[edge,edge+1];mt[np.ix_(ids,ids)]+=elem
    return mass*mass*np.kron(np.ones((q,q))/q,mt)


def stress_response(K,M,MT,KT,potential,w,cut=None):
    Q=np.eye(len(K)) if cut is None else null_space(cut)
    ev,U=eigh(Q.T@(K+potential)@Q,Q.T@M@Q)
    keep=ev>1e-8;freq=np.sqrt(ev[keep]);u=Q@U[:,keep]
    A=u.T@MT@u;B=u.T@(KT+potential)@u
    weights=graph.pair_weights(freq,A,B)
    total=freq[:,None]+freq[None,:]
    value=np.sum(2*w*w*weights/(total*(w*w-total*total)))
    imaginary=np.sum(2*weights/(total*(1+total*total)))
    assert imaginary>=0
    # Independent energy-representation formula for the entire complex sum.
    projection=(A-B/np.outer(freq,freq))/2
    alt=np.sum(pair_kernel(freq[:,None],freq[None,:],w)*abs(projection)**2)
    assert abs(value-alt)<2e-10
    return value,imaginary,len(freq)


def finite_element_checks():
    rows=[];T=PI;lengths=np.array([.7,2.1,12.8]);w=1.+.15j
    for n in [16,32,64]:
        K,M,C,CI,KT,MT=graph.assemble_circles(T,lengths,n)
        zero=np.zeros_like(K)
        free=stress_response(K,M,MT,KT,zero,w)
        freecut=stress_response(K,M,MT,KT,zero,w,C)
        for mass,variable in [(.2,False),(2.,False),(12.,False),(1.7,True)]:
            potential=mass_matrix(T,len(lengths),n,mass,variable)
            value=stress_response(K,M,MT,KT,potential,w)
            cut=stress_response(K,M,MT,KT,potential,w,C)
            cut_error=abs(cut[0]-freecut[0]);assert cut_error<2e-8
            difference=value[0]-free[0]
            rows.append(dict(cells_per_edge=n,mass=mass,variable_profile=variable,
                massive_response=[value[0].real,value[0].imag],free_response=[free[0].real,free[0].imag],
                difference=[difference.real,difference.imag],cut_mass_independence_error=cut_error,
                positive_imaginary_axis_value=value[1],retained_nonzero_scalar_modes=value[2]))
        print('Finite-element complex stress checks passed at',n,'cells per edge.',flush=True)
    refinements=[]
    for mass,variable in [(.2,False),(2.,False),(12.,False),(1.7,True)]:
        cases=[a for a in rows if a['mass']==mass and a['variable_profile']==variable]
        vals=[complex(*a['difference']) for a in cases]
        changes=[abs(vals[i+1]-vals[i]) for i in range(2)]
        assert changes[-1]<.02 and changes[-1]<changes[0]
        refinements.append(dict(mass=mass,variable_profile=variable,successive_mesh_changes=changes))
    return rows,refinements


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    exact=exact_checks();axes=axis_checks();annuli,tails,resonance=annular_checks();rank=finite_rank_checks()
    print('Kernel, axes, annuli and finite-rank checks passed.',flush=True)
    fem,refinements=finite_element_checks()
    r=dict(status='passed_singlet_complex_response_checks',date='2026-10-11',
      note_sha256=sha(ROOT/'notes/singlet-complex-response.md'),script_sha256=sha(Path(__file__)),
      dependency_hashes={f:sha(ROOT/f) for f in ['notes/singlet-stress-continuum.md','notes/finite-channel-decay.md','checks/schwinger_return_feedback.py']},
      versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
      exact_checks=exact,axis_cases=axes,annular_cases=annuli,tail_cases=tails,resonance_cases=resonance,
      finite_rank_cases=rank,finite_element_cases=fem,mesh_refinements=refinements,
      full_subtracted_gaussian_singlet_response_controlled=True,mass_uniform_bound='C(1+eta^-4)',
      damping_resolution_condition='q delta^5 -> infinity',
      static_curvature_bound_proved=False,nonspherical_interaction_bound_proved=False,
      physical_damping_proved=False,physical_state_constructed=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/singlet-complex-response.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=r['status'],counts={k:len(r[k]) for k in ['exact_checks','axis_cases','annular_cases','tail_cases','resonance_cases','finite_rank_cases','finite_element_cases','mesh_refinements']})))


if __name__=='__main__':main()
