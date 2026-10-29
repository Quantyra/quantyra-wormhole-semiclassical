"""Complete local/scalar mass residual, weighted original-domain estimates."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import argparse
import json
import sys
import traceback

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_mass_pairing_derivative_return as prev
m,fr,J=prev.m,prev.fr,prev.J
arb,ctx,dec,enc,sp=m.arb,m.ctx,m.dec,m.enc,m.sp
ROOT,BASELINE,runtime=m.ROOT,m.BASELINE,m.runtime
sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
STEM='current-mass-pairing-scalar-return'
digest=lambda p:sha256(p.read_bytes()).hexdigest()
up=prev.up
interval=prev.interval
poly=prev.polynomial


def bernstein(coeff):
    n=len(coeff)-1
    return [sum((coeff[j]*comb(k,j)/comb(n,j) for j in range(k+1)),arb(0))
            for k in range(n+1)]


def actual_reconstruction(p,bs,s0):
    L,ga,a=[p[k] for k in ['L','gamma','a']]
    muY=4*a/3;gap=muY-s0
    assert gap>0 and 1-s0>0
    R1=sum(((j-1)*abs(bs[j])*L**(j-3)/a for j in range(3,7)),arb(0))
    M1=(2*ga*L+a*L*L).sqrt();Dmax=arb('6.3')+arb('8.1')*L
    Mmax=arb('2.1')*L+arb('1.2');Gr=1/gap.sqrt()+s0.sqrt()/gap
    Knu=(2*L+(2*L*muY+Mmax)/gap+(14+L*R1)*M1*Gr+L*Dmax/gap)/(3*L**3)
    KT=(2*L+(2*Mmax/3+2*L*muY)/gap+(10+L*R1)*M1*Gr+2*L*Dmax/(3*gap))/L**3
    cy=arb(2)/3*(1+(s0.sqrt()+2/arb(5).sqrt())/(1-s0).sqrt())
    Kinv=1+cy/(5*(1-s0)).sqrt()+cy*(1+(s0/(1-s0)).sqrt())/2
    return Knu,KT,Kinv*KT


def local_coefficients(p,mp):
    core=SignedLocalCore();z=core.z
    values=dict(C=p['C'],B2=p['B2'],lambda0=arb(2000)**-10,
                cR=arb('1e22'),qflux=arb(2000)**4,gR=arb(2000)**-4,
                beta=p['beta'],eps=mp['eps'],tau=mp['tau'])
    subs=core.seed_substitution()
    arrays={name:core.original_A_coefficients(name) for name in core.densities}
    table=[];hashes={}
    for name in arrays:
        hashes[name]=sha256(str(arrays[name]).encode()).hexdigest()
    for i in range(2):
        row=[]
        for b in range(2):
            slot=[]
            for k in range(5):
                expr=sp.expand(-sum(arrays[name][i][b][k] for name in arrays).subs(subs))
                Lg=sp.Symbol('Lg')
                lp=sp.Poly(expr.subs(sp.log(z),Lg),Lg)
                assert lp.degree()<=1
                pieces=[]
                for n in range(2):
                    zp=sp.Poly(lp.nth(n),z)
                    degree=max(0,int(zp.degree())) if not zp.is_zero else 0
                    pieces.append([m.evaluate(zp.nth(j),values) for j in range(degree+1)])
                slot.append(pieces)
            row.append(slot)
        table.append(row)
    # Exact D^k(z^3 psi)/z^3, retaining both directions before norms.
    euler=[[1]]
    for k in range(4):
        old=euler[-1]
        euler.append([(j+3)*(old[j] if j<len(old) else 0)
                      +(old[j-1] if j else 0) for j in range(k+2)])
    return table,euler,hashes


def log_moments(z):
    if z.is_zero():return [arb(0)]*3
    l=z.log()
    return [z**3/3,z**3*(-l/3+arb(1)/9),
            z**3*(l*l/3-2*l/9+arb(2)/27)]


def integrate_sources(coeff,beta,L,ga,weightp,table,euler):
    cuts,Z,cutpolys=prev.cutoff_cells()
    cells=[(arb(i)/5120,arb(i+1)/5120,J([arb(1)]+[arb(0)]*4)) for i in range(512)]+cuts
    Wmax=((L+4*ga)/(L-arb('.2')+4*ga))**weightp
    totals=[arb(0)]*6
    maxrows=[arb(0),arb(0)]
    for count,(z0,z1,chi) in enumerate(cells):
        z=J([interval(z0,z1),arb(1)]+[arb(0)]*4)
        Q=J.cast(arb(0),6)
        z6=J(z.c+[arb(0)])
        for c in reversed(coeff):Q=Q*z6+c
        qp=Q.d();q=J(Q.c[:5])
        Zwide=J([interval(arb(0),z1),arb(1)]+[arb(0)]*4)
        Qwide=J.cast(arb(0),5)
        for c in reversed(coeff):Qwide=Qwide*Zwide+c
        assert Qwide.c[0]>arb('.15')
        gj=Qwide**-2
        I=J([gj.c[j]/(j+4) for j in range(5)])
        pn=3*beta*chi*q*I
        pt=-arb('1.5')*beta*chi*(1/q+z*qp*I)
        pu=pn+pt/2
        dirs=[]
        for psi in [pn,pt]:
            dirs.append([sum((euler[k][j]*z.c[0]**j*factorial(j)*psi.c[j]
                              for j in range(k+1)),arb(0)) for k in range(5)])
        logs=[b-a for a,b in zip(log_moments(z0),log_moments(z1))]
        assert all(x>0 for x in logs)
        for i in range(2):
            R=[sum((poly(table[i][b][k][n],z.c[0])*dirs[b][k]
                    for b in range(2) for k in range(5)),arb(0)) for n in range(2)]
            A,B=up(R[0]),up(R[1])
            totals[i]+=Wmax**2*(A*A*logs[0]+2*A*B*logs[1]+B*B*logs[2])
            maxrows[i]=prev.maxup(maxrows[i],A+B)
        zz=z.c[0];qq=q.c[0];qprime=qp.c[0]
        S=arb(2)/3*(zz*pt.c[1]-pt.c[0])
        localg=3*pu.c[0]+zz*pu.c[1]-S/4
        locall=zz**2*(3*qprime*pu.c[0]+(4*qq+zz*qprime)*pu.c[1]
                     +2*zz*qq*pu.c[2]-(qprime*S+arb(4)/3*zz*qq*pt.c[2])/4)
        dz=z1-z0
        totals[2]+=Wmax**2*up(pt.c[0])**2*logs[0]
        totals[3]+=Wmax**2*up(pu.c[0])**2*logs[0]
        totals[4]+=Wmax**2*up(qq)*up(localg)**2*logs[0]
        totals[5]+=Wmax**2*up(locall)**2*dz
        if count%512==511: print('complete support cells: '+str(count+1),flush=True)
    return [up(x).sqrt() for x in totals],Z,cutpolys,maxrows


def execute():
    ctx.prec=768
    p,geo,bs=fr.load()
    mp={k:dec(x) for k,x in m.result('current-global-mass-conversion')['parameters'].items()}
    v={k:prev.old.ball(prev.old.rational(x)) for k,x in
       m.result('current-correlated-signed-trial-test-attempt3')['exact_nominal_parameters'].items()}
    shape=m.constants('current-mass-pairing-derivative-return')
    bridge=m.constants('current-filtered-test-mass-trace-attempt2')
    adj=m.constants('current-filtered-mass-adjoint-bound-attempt2')
    mixed=m.constants('current-scalar-green-source-residual')
    green=m.constants('current-filtered-mixed-green-residual')
    L,ga,a,beta,d,C,B2=[p[k] for k in ['L','gamma','a','beta','d','C','B2']]
    coeff=[beta,arb(0),mp['u'],-4*beta-2*mp['u']-mp['tau'],
           3*beta+mp['u']+mp['eps'],-2*mp['eps'],mp['eps']]
    assert (96*arb.pi()*arb('1e22')*C/B2-C/(60*arb.pi()*B2)).overlaps(d)
    wp=arb(3)/4;sig=arb('1.001')*wp*wp
    assert bs[3]+sum((abs(bs[j])*L**(j-3) for j in range(4,7)),arb(0))<0
    # Exact polynomial guards on the WHOLE original interval [0,L] subset [0,1].
    guards={}
    for order,limit in [(0,arb('.17')),(1,arb('.3')),(2,arb('2.1')),(3,arb('8.1'))]:
        co=[coeff[j]*factorial(j)/factorial(j-order) for j in range(order,7)]
        for i in range(64):
            left=arb(i)/64;step=arb(1)/64
            shifted=[sum((co[j]*comb(j,k)*left**(j-k)*step**k for j in range(k,len(co))),arb(0))
                     for k in range(len(co))]
            assert all(abs(x)<limit for x in bernstein(shifted)),('global_geometry',order,i)
    for label,c1,c2,ceiling in [('QP_y',arb(2)/3,-arb(8)/3,arb(1)),
                               ('QP_xi',arb(1),-arb(4),arb('1.9'))]:
        co=[(6+c1*j*(j-1)+c2*j)*q for j,q in enumerate(coeff)]
        co[2]-=ceiling
        bern=bernstein(co);assert all(x>0 for x in bern),label
        guards[label]=[enc(x) for x in bern]
    # Accepted Pxi>=4.5/z2 guard and Py follows by Py=2Pxi/3+2/z2.
    assert arb('1.5')*beta-arb('4.5')*abs(mp['u'])-20*(abs(mp['tau'])+abs(mp['eps']))>0
    epsbar=C/(60*arb.pi()*B2);Cv=1+11*epsbar*beta/3
    w0=arb('1.1e-27')*beta;tau=arb('3.01')*beta
    KA=(2*d+11*epsbar/6)*w0+22*epsbar*tau/3
    cp=4*beta-sig;cd=Cv-d*sig
    assert cp>0 and cd>0
    GP=1/cp.sqrt()+sig.sqrt()/cp
    GD=1/(d*cd).sqrt()+sig.sqrt()/cd
    rA=(KA*GD*GP+2*w0/(cd*cp))/6
    assert rA<1
    Yc=1/(6*cd*(1-rA));Nk0=Yc/cp;NkG=Yc*GP;NkL=Yc*(1+4*beta/cp)
    Qratio=arb(1)/2+sum((j*abs(bs[j])*L**(j-2) for j in range(2,7)),arb(0))
    Do=wp*(L*L*arb('1.001')+L*L*Qratio+2*L*arb('1.001')*(L+4*ga))
    Gu=NkG+sig.sqrt()*Nk0
    Lu=NkL+2*sig.sqrt()*NkG+(Do+sig)*Nk0
    Dorb=Lu+2*arb('.17').sqrt()*Gu
    width=ga*Nk0/(4*Dorb);radius=arb('1e-110')
    assert 0<width<L
    Kobs=((4*ga+radius)/(L+4*ga))**wp*(4*(Nk0*Dorb/ga).sqrt()+4*radius.sqrt()*Dorb/ga)
    print('weighted scalar inverse guards established',flush=True)

    # Complete weighted rational trial on the entire auxiliary spectrum.
    degree=((1+4*v['muT']/a).sqrt()-1)/2
    assert 2*degree+1-2*wp>0 and (v['muT']*radius/(2*ga)).exp()<2
    BQ2=(4+1/a)**(2*wp)*(arb(3).log()**2+2*arb(3).log()+2)/(4*a)
    BQ2+=arb.pi()**2*(1+4*a)**(2*wp)/(a**(2*wp+1)*(2*degree+1-2*wp))
    A0=3*v['alpha']/(4*v['gamma']*v['A'])
    t=v['twice_a'];assert t>sig and v['muS']>sig
    Kfw=2*BQ2.sqrt()*A0*(1+abs(v['muY']-t)/(t-sig))*(1+2*v['scalar_c']/(v['muS']-sig))
    h=v['d1']/2;uW=abs(h*h-v['d0'])/(h-sig)**2
    assert h>sig and uW<1 and v['Cv']-v['delta']*h>0
    LFw=v['kappa']*((v['Cv']-v['delta']*h)/(h-sig)**2+v['delta']/(h-sig))/(1-uW)
    Nfw0=ga**(wp+arb(1)/2)*Kfw/(L+4*ga)**wp;Nfw1=LFw*Nfw0
    Knu,KT,Kxi=actual_reconstruction(p,bs,sig)
    Knu0,KT0,Kxi0=actual_reconstruction(p,bs,arb(0))
    Vw0=Kxi*Nfw0;Vw1=Kxi*Nfw1
    V00=Kxi0*shape['original_trial0_profile_L2_upper']
    V01=Kxi0*shape['original_trial1_profile_L2_upper']
    kx=1-sig/arb('1.9')
    cY=1+((sig/arb('1.9')).sqrt()+2/arb('4.5').sqrt())/kx.sqrt()
    yz=cY/(arb('4.5')*kx).sqrt()
    yp=cY*(1+(sig/(arb('1.9')*kx)).sqrt())
    CUw=arb(1)/2+yz;CVw=1+yz+yp/2
    print('complete weighted rational and York comparisons established',flush=True)
    table,euler,hashes=local_coefficients(p,mp)
    print('both complete physical local rows assembled',flush=True)
    (Fn,Ft,TH,U0H,Gloc,Lloc),Z,cutpolys,maxrows=integrate_sources(coeff,beta,L,ga,wp,table,euler)
    ky=1-sig
    Ey=arb(2)/3*TH*(1+(2+arb('.2')*wp/(L-arb('.2')+4*ga))/(5*ky).sqrt())
    Yz=Ey/(5*ky).sqrt();Yp=(1+(sig/ky).sqrt())*Ey;YP=Ey/ky.sqrt()
    Ns0=U0H+Yp/4+Yz/2
    NsG=Gloc+arb('3.22').sqrt()*YP/4+arb('.17').sqrt()*Yz/2
    NsL=Lloc+arb('3.59')*Yz+arb('.72')*Yp
    Esm=((6*d*NsL**2+(6*Cv+24*beta*d+KA)*NsG**2+(24*beta*Cv+2*w0)*Ns0**2)/beta).sqrt()
    Ekm=Fn/(beta*arb('23.99').sqrt())
    Kbrem=mixed['complete_mixed_remainder_constant']
    Bs_main=4*tau/beta*(Cv*Ns0+d*NsL)*Vw0
    Bs_rem=Kbrem/beta*Esm*V00
    Bk_main=4*tau/beta*(Cv*Nk0+d*NkL)*Fn*(Vw0+Vw1)
    Bk_rem=Kbrem/beta*Ekm*(V00+V01)
    q=arb('.338');Kg=adj['actual_filtered_mass_functional_norm_upper']
    alpha=green['alpha_original_radius_cost']
    Hloc=(Fn*Knu+Ft*KT)/beta*(q*Nfw0+Nfw1)
    Obs=Kg*q*alpha*Kobs*Fn
    Smain=Kg*(q+1)*Bs_main;Srem=Kg*(q+1)*Bs_rem
    Kmain=Kg*Bk_main;Krem=Kg*Bk_rem;Hdebit=Kg*Hloc
    Elocal_generic=Obs+Smain+Srem+Kmain+Krem+Hdebit
    # SAME full mixed estimate, using the equation's actual derivatives
    # before introducing the generic energy loss d^(-1/2).
    Rjets=arb('3e6')*(4400+2*4*2*arb(2).sqrt()*20+8)
    Cjets=100*Rjets;w2=arb('4.3e-25')*beta**2
    def brem(N0,NG,NL):
        NH=(NL*NL+13*beta*NG*NG).sqrt()
        curvature=arb(22).sqrt()*(100*Cjets*(NH+2*NG+N0)+arb('1e6')*N0)
        return (4*tau*(11*epsbar*w0/36)*N0+arb(22).sqrt()*d*w2*N0+epsbar*curvature
                +20/(60*arb.pi())*arb('4e-27')*beta*N0
                +arb('9e-28')*beta*(5/beta.sqrt())*NG)
    Bsrem_direct=brem(Ns0,NsG,NsL)*V00/beta
    Bkrem_direct=brem(Fn*Nk0,Fn*NkG,Fn*NkL)*(V00+V01)/beta
    Sdirect=Kg*(q+1)*Bsrem_direct;Kdirect=Kg*Bkrem_direct
    Elocal=Obs+Smain+Sdirect+Kmain+Kdirect+Hdebit
    Eshape=shape['complete_massive_and_LLL_mass_pairing_upper'];Emass=Eshape+Elocal
    conv=bridge['full_E_to_retained_Et_norm_comparison'];Eb=bridge['full_to_retained_j_energy_error']
    target=arb('1.49e34');remaining=(target-Eb)/conv-Emass
    # Absolute a0/a1 intervals are valid; no uncomputed sign is claimed.
    local0=alpha*Kobs*Fn+(Fn*Knu+Ft*KT)*Nfw0/beta+Bs_main+Bsrem_direct
    local1=Bs_main+Bsrem_direct+Bk_main+Bkrem_direct+(Fn*Knu+Ft*KT)*Nfw1/beta
    vals=dict(weight_exponent=wp,weight_gradient_square=sig,scalar_P_gap=cp,scalar_D_gap=cd,
        complete_scalar_remainder_inverse_ratio=rA,scalar_weighted_L2_coefficient=Nk0,
        scalar_weighted_gradient_coefficient=NkG,scalar_weighted_Laplacian_coefficient=NkL,
        weight_Laplacian_upper=Do,observed_scalar_return_coefficient=Kobs,
        scalar_bolt_average_upper=Kobs*Fn,scalar_observation_interval_width=width,
        auxiliary_Q_weighted_norm_upper=BQ2.sqrt(),complete_weighted_Fref_ratio=uW,
        complete_weighted_Fref_norm=LFw,original_weighted_trial0_Hmu_upper=Vw0,
        original_weighted_trial1_Hmu_upper=Vw1,original_unweighted_trial0_Hmu_upper=V00,
        original_unweighted_trial1_Hmu_upper=V01,weighted_York_U_constant=CUw,
        weighted_York_V_constant=CVw,cutoff_normalization=Z,
        complete_weighted_canonical_nu_reconstruction=Knu,
        complete_weighted_canonical_T_reconstruction=KT,
        complete_weighted_TT_reconstruction=Kxi,complete_unweighted_TT_reconstruction=Kxi0,
        complete_local_scalar_source_weighted_Hmu_upper=Fn,
        complete_local_shape_source_weighted_Hmu_upper=Ft,mass_T_weighted_Hmu_upper=TH,
        mass_scalar_weighted_Hmu_upper=Ns0,mass_scalar_weighted_gradient_upper=NsG,
        mass_scalar_weighted_Laplacian_upper=NsL,mass_scalar_energy_upper=Esm,
        scalar_return_energy_upper=Ekm,complete_mixed_remainder_constant=Kbrem,
        scalar_observation_mass_debit=Obs,scalar_mass_mixed_main_debit=Smain,
        scalar_mass_mixed_remainder_debit=Srem,scalar_return_mixed_main_debit=Kmain,
        scalar_return_mixed_remainder_debit=Krem,complete_local_canonical_pairing_debit=Hdebit,
        generic_energy_local_scalar_mass_block_upper=Elocal_generic,
        generic_energy_complete_mass_estimator=Eshape+Elocal_generic,
        generic_energy_estimator_over_whole_Dj_target=conv*(Eshape+Elocal_generic)/target,
        direct_scalar_mass_mixed_remainder_debit=Sdirect,
        direct_scalar_return_mixed_remainder_debit=Kdirect,
        direct_mixed_remainder_gain=(Srem+Krem)/(Sdirect+Kdirect),
        complete_local_scalar_mass_block_upper=Elocal,
        local_residual0_absolute_upper=local0,local_residual1_absolute_upper=local1,
        unchanged_massive_and_LLL_mass_block_upper=Eshape,complete_mass_vector_norm_upper=Emass,
        complete_mass_vector_bolt_units=green['pref_t']*conv*Emass,
        mass_fraction_of_Dj_diagnostic=conv*Emass/target,
        conditional_full_E_account_for_unknown_mass_free_and_errors=remaining,
        full_E_to_retained_Et_comparison=conv,unchanged_full_response_bridge=Eb)
    assert all(x.is_finite() for x in vals.values())
    return dict(constants={k:enc(x) for k,x in vals.items()},whole_cell_count=1536,
        exact_local_operator_sha256=hashes,exact_euler_coefficients=euler,
        global_polynomial_Bernstein_guards=guards,cutoff_derivative_polynomials=cutpolys,
        numerical_contract_pass=bool(remaining>0),actual_mass_free_residual=None,
        remaining_application_errors=None,actual_Dj=None,actual_Gt=None,
        scope=dict(complete_local_scalar_mass_contribution_bounded=True,
            complete_mass_contribution_bounded=True,actual_paired_fixed_seed_domain=True,
            complete_signed_response=False,complete_Dj_diagnostic=False,
            actual_signed_Gt=False,full_bolt_allowance=False,available_physical_margin=False,
            original_P32_defect=False,actual_kappa_m=False,all_original_coupled_outputs=False,
            nonlinear_admission_contraction=False,nonlinear_neighborhood_transfer=False,
            nonlinear_existence=False,exact_compatibility=False,physical_matching=False,
            semiclassical_validity_noise=False,all_sector_stability=False,
            construction_obstruction=False,full_goal_complete=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
    out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable receipts'
    prior=json.loads((ROOT/'evidence/current-mass-pairing-derivative-return.json').read_bytes())
    scalar=json.loads((ROOT/'evidence/current-scalar-green-source-residual.json').read_bytes())
    scalar=scalar.get('runs',[scalar])[-1]
    inputs=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_mass_pairing_scalar_return.py',
        'operators/signed_local_core.py','evidence/current-mass-pairing-derivative-return.json',
        'evidence/current-mass-pairing-derivative-return-verification.json',
        'evidence/current-scalar-green-source-residual.json',*prior['input_sha256_after'],
        *scalar['input_sha256_after']]))
    before={p:digest(ROOT/p) for p in inputs};protected=json.loads(BASELINE.read_bytes())
    rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        assert all(digest(ROOT/p)==h for p,h in prior['input_sha256_after'].items())
        assert all(digest(ROOT/p)==h for p,h in scalar['input_sha256_after'].items())
        assert all(digest(Path(p))==h for p,h in protected.items())
        assert digest(Path(sys.executable))==runtime['python_sha256']
        for pkg in runtime['packages'].values():
            assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
        rec['result']=execute()
        after={p:digest(ROOT/p) for p in inputs};assert before==after
        assert all(digest(Path(p))==h for p,h in protected.items())
        rec.update(status='complete_mass_contribution_bounded_full_diagnostic_open',
            input_sha256_after=after,protected_paths_unchanged=len(protected),
            python_sha256=runtime['python_sha256'],
            pinned_packages={k:len(v['files']) for k,v in runtime['packages'].items()})
    except Exception:rec.update(status='failed',traceback=traceback.format_exc())
    rec['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
        numerical_contract_pass=rec.get('result',{}).get('numerical_contract_pass'),
        constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
    if rec['status']=='failed':raise SystemExit(1)
