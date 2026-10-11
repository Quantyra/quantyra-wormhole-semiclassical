"""Exact conditional energy structure for the spherical drift-response model.

No expansion in the background is made. The residual background stress is
defined by Einstein's equation, and its perturbation is prescribed by the
local algebraic drift lift. This checks that model, not its derivation from
four-dimensional quantum field theory or the existence of an exact vacuum.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp

spec=importlib.util.spec_from_file_location('throat_geometry',Path(__file__).with_name('throat_geometry.py'))
geometry=importlib.util.module_from_spec(spec);spec.loader.exec_module(geometry)


def verify():
    ((r,u,L,s,p),(h0,h1,w0,w1,X),(S,C),(Y,Yp),
     (simp,derivative),(metric,inverse,h,dinverse),
     (connection,dconnection),(field,field_up),
     (tq,einstein,maxwell))=geometry.derive_geometry()
    k=L-2
    results={}
    def zero(name,v):
        value=simp(v)
        assert value==0,(name,value)
        results[name]='exact_zero'
        print(name,'exact_zero',flush=True)

    pressure=sp.factor(tq[2,2]/metric[2,2])
    kval=[sp.factor(inverse[a,a]*tq[a,a]-pressure) for a in [0,1]]
    at=sp.factor(1-C*C*kval[0]/2)
    ar=sp.factor(1-C*C*kval[1]/2)
    zero('temporal_response_geometry',at-C*(C*sp.diff(S,r,2)-S*sp.diff(C,r,2)+2)/4)
    zero('radial_response_geometry',ar-(C*C*sp.diff(S,r,2)+C*S*sp.diff(C,r,2)+2*C-S*sp.diff(C,r)**2)/4)
    radial_null_stress=tq[0,0]/S+S*tq[1,1]
    zero('response_difference_null_energy',ar-at+C*C*radial_null_stress/2)
    zero('radial_null_focusing',radial_null_stress+S*(2*C*sp.diff(C,r,2)-sp.diff(C,r)**2)/(2*C*C))
    zero('local_characteristic_cone',
         (S/at)/(1/(S*ar))/S**2-1+C*C*radial_null_stress/(2*at))
    dtq=sp.zeros(4)
    for a,ha,wa in [(0,h0,w0),(1,h1,w1)]:
        tilt=-(ha+C*wa/2)*(1-u*u)*Yp
        dtq[a,3]=dtq[3,a]=pressure*h[a,3]+kval[a]*tilt
    coeff=[inverse[0,0]*kval[0]*(h0+C*w0/2),inverse[1,1]*kval[1]*(h1+C*w1/2)]
    jpol=sp.Matrix([-C*L*coeff[0]*Y,-C*L*coeff[1]*Y,
                    -(1-u*u)*(p*C*coeff[0]+sp.diff(C*coeff[1],r))*Yp,0])
    jsch=sp.Matrix([s*sp.diff(X,r)*Y,-s*p*X*Y,0,0])
    zero('polarisation_current_conservation',sum(derivative(jpol[a],a,True) for a in range(4)))
    zero('Schwinger_current_conservation',sum(derivative(jsch[a],a,True) for a in range(4)))
    mixed=inverse*tq
    dmixed=inverse*dtq+dinverse*tq
    for nu in range(4):
        div=sum(derivative(dmixed[mu,nu],mu,True) for mu in range(4))
        for mu in range(4):
            for lam in range(4):
                div+=connection[mu,mu,lam]*dmixed[lam,nu]-connection[lam,mu,nu]*dmixed[mu,lam]
                div+=dconnection[mu,mu,lam]*mixed[lam,nu]-dconnection[lam,mu,nu]*mixed[mu,lam]
        zero('quantum_force_balance_'+str(nu),div+sum(field[nu,lam]*jpol[lam]/(2*C) for lam in range(4)))

    electric=p*w1-sp.diff(w0,r)
    displacement=C*electric-s*X
    grav=sp.diff(h0,r)-p*h1-sp.diff(C,r)*h0/C
    b0,b1=w0+2*h0/C,w1+2*h1/C
    predicted_e0=-S*sp.diff(C*grav,r)/(2*C)+k*h0/(2*C)+at*b0/C
    predicted_e1=-p*grav/(2*S)+k*h1/(2*C)+ar*b1/C
    angular=-(1-u*u)*Yp
    zero('Einstein_tau_equation',(einstein[0,3]-dtq[0,3])-angular*predicted_e0)
    zero('Einstein_rho_equation',(einstein[1,3]-dtq[1,3])-angular*predicted_e1)
    zero('Maxwell_tau_equation',maxwell[0]-jsch[0]-jpol[0]-Y*(sp.diff(displacement,r)+L*at*b0/S))
    zero('Maxwell_rho_equation',maxwell[1]-jsch[1]-jpol[1]-Y*(-p*displacement-L*S*ar*b1))
    residual=maxwell-jsch-jpol
    zero('Maxwell_angular_constraint',p*residual[0]/Y+sp.diff(residual[1]/Y,r)
         -L*residual[2]/((1-u*u)*Yp))
    zero('Einstein_angular_equation',einstein[2,3]
         -(-L*Y+2*u*Yp)/(2*S)*(p*h0-S*sp.diff(S*h1,r)))
    zero('Maxwell_phi',residual[3])
    for i in range(4):
        for j in range(i,4):
            if (i,j) not in [(0,3),(1,3),(2,3)]:
                zero('other_Einstein_'+str(i)+str(j),einstein[i,j]-dtq[i,j])

    # Eliminate constraints using independent symbolic response coefficients.
    # Their specific geometric formulas are needed for the tensor match
    # above, but not for the following action and positivity identities.
    aa,bb=sp.Function('A_t')(r),sp.Function('A_r')(r)
    F,Phi,chi=[sp.Function(n)(r) for n in ['F','Phi','chi']]
    h0r=S*sp.diff(F,r)/k
    h1r=p*F/(k*S)
    w0r=-2*h0r/C-S*sp.diff(Phi,r)/(L*aa)
    w1r=-2*h1r/C-p*Phi/(L*S*bb)
    eom_f=p*p*F/(S*C)-sp.diff(S*sp.diff(F,r)/C,r)+k*(F-2*Phi/L)/C**2
    eom_phi=p*p*Phi/(S*bb)-sp.diff(S*sp.diff(Phi,r)/aa,r)+(L/C+4/C**2)*Phi-2*L*F/C**2+L*s*chi/C
    eom_chi=p*p*chi/S-sp.diff(S*sp.diff(chi,r),r)+(Phi+s*chi)/C
    gr=sp.diff(h0r,r)-p*h1r-sp.diff(C,r)*h0r/C
    er=p*w1r-sp.diff(w0r,r)
    zero('gravitational_definition_residual',C*gr+2*Phi/L-F+C*C*eom_f/k)
    zero('electric_definition_residual',er-(Phi+s*chi)/C+eom_phi/L+2*eom_f/k)
    # Check both remaining Einstein equations, including their off-shell
    # residuals, rather than relying only on their constraint derivation.
    et=-S*sp.diff(C*gr,r)/(2*C)+k*h0r/(2*C)+aa*(w0r+2*h0r/C)/C
    erho=-p*gr/(2*S)+k*h1r/(2*C)+bb*(w1r+2*h1r/C)/C
    zero('Einstein_tau_reconstruction',et-S*sp.diff(C*C*eom_f,r)/(2*k*C))
    zero('Einstein_rho_reconstruction',erho-p*C*eom_f/(2*k*S))
    zero('charge_anomaly_reconstruction',
         -p*p*chi/S+sp.diff(S*sp.diff(chi,r),r)-er
         +eom_chi-eom_phi/L-2*eom_f/k)
    zero('Einstein_angular_reconstruction',p*h0r-S*sp.diff(S*h1r,r))
    zero('Maxwell_tau_reconstruction',sp.diff(Phi,r)+L*aa*(w0r+2*h0r/C)/S)
    zero('Maxwell_rho_reconstruction',-p*Phi-L*S*bb*(w1r+2*h1r/C))

    # Euler-Lagrange equations and the local conservation law in the time domain.
    # Positive energy is one half of kinetic + gradient + potential squares.
    v=sp.Matrix(sp.symbols('F0 Phi0 X0',real=True))
    potential=(L*v[0]-2*v[1])**2/C**2+L*(v[1]+s*v[2])**2/C
    matrix=sp.hessian(potential,tuple(v))/2
    weight=sp.diag(L*L/k,1,L*s)
    couplings=sp.Matrix([[k/C**2,-2*k/(L*C**2),0],[-2*L/C**2,L/C+4/C**2,L*s/C],[0,1/C,s/C]])
    for i in range(3):
        for j in range(3):zero('potential_hessian_'+str(i)+str(j),matrix[i,j]-(weight*couplings)[i,j])
    tau=sp.Symbol('tau',real=True)
    fields=sp.Matrix([sp.Function(n)(tau,r) for n in ['F_t','Phi_t','X_t']])
    kinetic=sp.diag(L*L/(k*S*C),1/(S*bb),L*s/S)
    gradient=sp.diag(L*L*S/(k*C),S/aa,L*s*S)
    dot=fields.diff(tau);prime=fields.diff(r)
    density=((dot.T*kinetic*dot)[0]+(prime.T*gradient*prime)[0]+(fields.T*matrix*fields)[0])/2
    flux=-(dot.T*gradient*prime)[0]
    equations=kinetic*fields.diff(tau,2)-(gradient*prime).diff(r)+matrix*fields
    zero('local_energy_conservation',sp.diff(density,tau)+sp.diff(flux,r)-(dot.T*equations)[0])
    null=sp.Matrix([-2*s/L,-s,1])
    for i in range(3):zero('constant_shift_zero_mode_'+str(i),(matrix*null)[i])

    # The perturbative MMP throat has positive response coefficients for
    # small alpha. These expansions are not a global background existence proof.
    alpha=sp.Symbol('alpha',real=True);f=1+r*r
    gamma=-4*r*r/3-r*(r*r+3)*sp.atan(r)+sp.log(f)
    phi=r*sp.atan(r)+sp.Rational(4,3)
    substitutions={sp.diff(old,r,n):sp.diff(new,r,n) for old,new in [(S,f+alpha*gamma),(C,1+alpha*phi)] for n in range(3)}
    for name,actual,expected in [('A_t',at,(2*r*r-1)/(6*f)),('A_r',ar,(2*r*r+5)/(6*f))]:
        raw=actual.xreplace(substitutions)
        zero(name+'_MMP_leading',raw.subs(alpha,0)-1)
        zero(name+'_MMP_first_order',sp.diff(raw,alpha).subs(alpha,0)-expected)
    mass=sp.Symbol('M',real=True)
    mouth=[(S,1-2*mass/r+1/r**2),(C,r**2)]
    subs_mouth={sp.diff(old,r,n):sp.diff(new,r,n) for old,new in mouth for n in range(3)}
    zero('mouth_A_t',at.xreplace(subs_mouth)-1)
    zero('mouth_A_r',ar.xreplace(subs_mouth)-1)
    # A smooth two-ended coefficient family makes the exact hypotheses
    # nonempty. Its residual source is not asserted to be a fermion vacuum.
    bridge=sp.Symbol('a',positive=True)
    radius=1+sp.sqrt(r*r+bridge*bridge)
    sc=(1-1/radius)**2;cc=radius**2
    subs_bridge={sp.diff(old,r,n):sp.diff(new,r,n) for old,new in [(S,sc),(C,cc)] for n in range(3)}
    zero('smooth_example_A_t',at.xreplace(subs_bridge)-1-bridge**2/sp.sqrt(r*r+bridge**2))
    zero('smooth_example_A_r',ar.xreplace(subs_bridge)-1-bridge**2-2*bridge**2/sp.sqrt(r*r+bridge**2))
    return {'identities':results,
            'model':'Exact prescribed drift lift of a static spherical residual stress, conserved polarisation completion, Schwinger scalar',
            'positive_energy_hypotheses':['S>0','C>0','A_t>0','A_r>0','L>2','s>0','static coefficients','admissible vanishing or outward energy flux'],
            'time_domain_reconstruction_requires_nonzero_frequency':False,
            'constant_kernel':'(F,Phi,X)=constant*(-2*s/L,-s,1); reconstructed metric and electromagnetic perturbations vanish',
            'relative_radial_characteristic_speed_squared':['1','A_r/A_t','1'],
            'cone_identity':'A_r/A_t-1=-C^2*(kappa T_Q(n,n))/(2 A_t); n=S^(-1/2) partial_t+S^(1/2) partial_r',
            'causal_scope':'At a strict flaring minimum of R=sqrt(C), positive A_t implies a super-metric characteristic in this exact local closure. No microscopic front-velocity or external shortcut claim.',
            'A_t_geometric':str(at),'A_r_geometric':str(ar),
            'smooth_example':'R=1+sqrt(rho^2+a^2), C=R^2, S=(1-1/R)^2, a>0; A_t>1 and A_r>1 everywhere',
            'smooth_example_source':'Einstein residual with the prescribed drift perturbation; not a derived quantum vacuum',
            'microscopic_quantum_closure_proved':False,'exact_semiclassical_background_existence_proved':False}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_exact_conditional_drift_energy_structure',
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'dependency_sha256':{'throat_geometry':hashlib.sha256(Path(__file__).with_name('throat_geometry.py').read_bytes()).hexdigest()},
          'sympy_version':sp.__version__,'checks':verify(),
          'physical_MMP_stability_proved':False,'novelty':'unresolved','lean_certified':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
