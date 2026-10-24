"""Whole-domain error accounting for the signed resolvent enclosure.

The ordinary argument is in current-signed-resolvent-projection.md.
This module never selects physical thermal modes or changes the boundary.
"""
from math import factorial,comb

def certify(c,p,P,V,N,PB,VB,R,RV,J,JV,AJ,VJ,Fmax):
 a=c.arb; enc=c.encoded; dec=c.decode; check=c.check
 pi=a.pi(); M1=((a(p['q'])+1)/p['C']).sqrt()
 H=a('1.2e-8');fm=a('.409');AA=a('1.01');cc=AA/4;cv=(a('.5')+AA)/4
 B=a('.008');BV=a('.02');Gnorm=1/(pi*p['B2']*p['Fl']);yc=a(2)**-49
 tag=lambda s:s+'_'+str(N)
 geo=next(r for r in reversed(c.NUMERICS) if r['kind']=='complete_joint_radial_coefficient_bounds' and r['N']==N)
 L=dec(geo['proper_bulk_length_upper']);Hr=H*L
 check(tag('actual_geometry_lipschitz_guards'),AJ[0]<AA and VJ[1]<a('.5') and Fmax<AA)
 check(tag('physical_f_guard'),p['beta']+p['u']/16<fm**2 and p['eps']<0 and p['tau']>0)
 check(tag('approximation_denominator_guard'),J<a('.1') and JV<a('.1'))
 eta=cc*max(1+J,2*J)/(2*M1);b0=2*(1-eta)
 Ca=(1+(cv*(1+J)+5*cc)/(2*M1))/(1-eta)
 buffer=7*a(2).log()/AA;half=buffer/2
 xe=cc/(2*M1-cc);xbound=xe+(1+xe)*(-(2*M1-cc)*half).exp()
 ve=(B+cv*(1+B**2)/(2*M1))/(1-cc*B/M1)
 vbound=ve+(a('2.5')+ve)*(-(2*M1-2*cc*B)*half).exp()
 check(tag('true_regular_buffer_smallness'),xbound<B and vbound<BV)
 masses=[(a(n)*(n+p['q'])/p['C']).sqrt() for n in range(1,81)]
 ds=[a(2)*(p['q']+2*n) for n in range(1,81)]
 sm=lambda k:sum((d/M**k for d,M in zip(ds,masses)),a(0))
 je=lambda e:pi.sqrt()*(e-a('.5')).gamma()/(2*e.gamma())
 def profile(power,exponent):
  # gamma sum of a monotone profile, with its entire tail.
  return fm*je(a(exponent))*sm(power-1)+p['gamma']*sm(power)/2
 Kerr=((9*Ca+9+36*BV)*R+6*RV)/(2*(1-eta))
 Erad=Gnorm*Hr*Kerr*profile(N,a(N+2)/2)

 # Nonnegative rational majorant, exact to ALL orders.
 lim=6
 def mul(x,y):return [sum((x[j]*y[k-j] for j in range(k+1)),a(0)) for k in range(lim+1)]
 plus=lambda x,y:[xx+yy for xx,yy in zip(x,y)]
 scale=lambda x,b:[xx*b for xx in x]
 one=[a(1)]+[a(0)]*lim
 qa=mul(PB[:lim+1],PB[:lim+1]);qva=scale(mul(PB[:lim+1],VB[:lim+1]),2)
 for j in range(1,lim+1,2):qa[j]=a(0);qva[j]=a(0)
 sa=[2*PB[j] if j%2==0 else a(0) for j in range(lim+1)]
 sva=[2*VB[j] if j%2==0 else a(0) for j in range(lim+1)]
 na=plus(scale(qa,4),scale(sa,2))
 nva=plus(plus(scale(qva,4),scale(qa,8)),plus(scale(sva,2),scale(sa,2)))
 q=scale(qa,a(1)/4);qv=scale(qva,a(1)/4);qq=mul(q,q);qqq=mul(qq,q)
 inv=plus(plus(plus(one,q),qq),qqq)
 inv2=plus(plus(plus(one,scale(q,2)),scale(qq,3)),scale(qqq,4))
 hm=plus(mul(nva,inv),mul(mul(na,qv),inv2))
 # Even parts without cancellation: sum over all i+j even.
 qval=sum((PB[i]*PB[j]/M1**(i+j) for i in range(1,N+1) for j in range(1,N+1) if (i+j)%2==0),a(0))
 qvval=sum((2*PB[i]*VB[j]/M1**(i+j) for i in range(1,N+1) for j in range(1,N+1) if (i+j)%2==0),a(0))
 sval=2*sum((PB[j]/M1**j for j in range(2,N+1,2)),a(0))
 svval=2*sum((VB[j]/M1**j for j in range(2,N+1,2)),a(0))
 qn=qval/4;check(tag('rational_majorant_denominator'),qn<a('.01'))
 hval=(4*qvval+8*qval+2*svval+2*sval)/(1-qn)+(4*qval+2*sval)*qvval/4/(1-qn)**2
 Ktail=M1**8*(hval-sum((hm[j]/M1**j for j in [2,4,6]),a(0)))
 check(tag('positive_entire_rational_tail'),Ktail>0)
 Eapp=Gnorm*Ktail*(H*a(16)/315*sm(6)+p['gamma']*Hr*sm(7))

 # Complete midpoint rule error for the retained local polynomials.
 # Each H_2k has degree <=12. Markov on [0,1]: |H'|<=288 sup|H|.
 # The weight r^2(1-r^2)^(k+.5) has TV<=2, so TV <=290 Hmajor.
 Etherm=Gnorm*p['gamma']*Hr*(sum((145*hm[2*k]*sm(2*k-1) for k in [1,2,3]),a(0))+3*sm(-1))

 # Inner data remain an enclosure of the true regular branch, no wall.
 amin=b0*M1*buffer
 check(tag('inner_transient_profile_monotonicity'),amin>13)
 e0=1+J;ev0=a('2.5')+JV
 Tinner=(3*ev0+(3*Ca*amin+3+12*BV)*e0)*(-amin).exp()
 # Right-inherited transient and its sensitivity receive a factor four.
 # x exp(-x) and exp(-x) retain the same lmax^6 envelope for x>=amin.
 Einner=4*Gnorm*Hr*Tinner*profile(-1,a(6))

 # Boundary jets computed directly from the original Q(z), not reflected jets.
 order=14
 qs=[p['beta'],a(0),p['u'],-4*p['beta']-2*p['u']-p['tau'],3*p['beta']+p['u']+p['eps'],-2*p['eps'],p['eps']]+[a(0)]*(order-7)
 ff=[qs[0].sqrt()]
 for k in range(1,order):ff.append((qs[k]-sum((ff[j]*ff[k-j] for j in range(1,k)),a(0)))/(2*ff[0]))
 avec=[-(j+1)*ff[j+1] for j in range(order-1)]
 vvec=[a(0)]*4+[a(-1)]+[a(0)]*(order-5)
 ajet=[];vjet=[]
 for j in range(N+1):
  ajet.append(avec[0]);vjet.append(vvec[0]);avec=c.jet_D(avec,ff);vvec=c.jet_D(vvec,ff)
 values=ajet+[a(0)]*(len(c.a)-len(ajet))+vjet+[a(0)]*(len(c.v)-len(vjet))
 SB=[a(0)]+[c.eval_poly(c.compile_poly(P[j]),values) for j in range(1,N+1)]
 SD=[a(0)]+[c.eval_poly(c.compile_poly(V[j]),values) for j in range(1,N+1)]
 check(tag('exact_seam_low_odd_jets'),SB[1].is_zero() and all(SD[j].is_zero() for j in [1,3,5]))
 Eb=a(0);bterms=[]
 # Rows (coefficient, power of distance, inverse power of M) for
 # 3 Q(d)+12 P(d)+3 d^4 in the note's boundary proof.
 rows=[(3*(B+J)/4,4,0),(3*(B+J)*AA/20,5,0),
       (3*cc,5,0),(6*cc,4,1),(12*cc,3,2),(18*cc,2,3),(18*cc,1,4),
       (a(33),4,0),(a(48),3,1),(a(72),2,2),(a(72),1,3),(a(36),0,4)]
 for j in range(1,N+1,2):
  var=3*SD[j]/(1-eta)*profile(j,a(j+2)/2)
  val=2*SB[j]*fm**4*(3/(1-eta))*144/b0**5*profile(j+4,a(j+6)/2)
  for coef,k,m in rows:
   val+=2*SB[j]*fm**4*coef*factorial(k)/b0**(k+1)*profile(j+k+m,a(j+k+2)/2)
  Eb+=val+var;bterms.append(dict(j=j,P=enc(SB[j]),V=enc(SD[j]),value_error=enc(Gnorm*H*val),variation_error=enc(Gnorm*H*var)))
 Eboundary=Gnorm*H*Eb

 # Actual cold observation, all thermal and angular modes majorize n<=80.
 FP=c.prior.FP;mon=c.prior.mon
 zz=mon(1);Q=(1-zz)**2*(p['beta']*(1+2*zz+3*zz**2)+p['u']*zz**2+p['eps']*zz**4)-p['tau']*zz**3
 vv=-zz**4;DD=lambda x:(Q*x.diff()).diff();RR=-Q.diff().diff();LL=lambda x:DD(x)+RR*x
 O4=DD(DD(vv))+RR*DD(vv)+DD(RR*vv)+RR**2*vv/2+2*vv
 O6=-(LL(DD(LL(vv)))+DD(RR)*LL(vv)+DD(vv*DD(RR)))/105-(16*RR**3*vv+24*RR**2*DD(vv)+24*DD(RR**2*vv)+24*DD(RR*DD(vv)))/1008
 def sup_poly(x):return sum((abs(co) for co in x.c.values()),a(0))
 sc=yc/(1-yc);Hc=a('1e-6')*sc*(1+(a('.01')/sc).log())
 Uc=Hc*yc*sup_poly(O6+4*vv/63)/(4*pi);Bc=1101*Hc;t0=a('.001');Ic=Hc*yc
 S1=Hc*c.prior.absolute_poly_integral(DD(vv)-2*vv)/(12*pi)
 S2=Hc*c.prior.absolute_poly_integral(O4)/(120*pi)
 Ch=2*a(2).sqrt()*a('2.25')*1521*Hc/a(1).exp();lam=M1**2
 theta=(-lam*t0).exp()*(2*(p['q']+2)+2*p['C']/t0)
 coldpos=2*theta/(p['B2']*p['Fl'])*(Uc*t0**2/4+Bc*t0**3/6+Ch*t0.exp()*lam/(lam-1)+Ic/(pi*t0)+Ic/(6*pi)+S1+S2*(t0+1/lam))
 colduv=4*p['C']/(p['B2']*p['Fl'])*(Uc*t0+Bc*t0**2/2)
 coldproduct=a('2.984e-37')*Ic
 # Raw A3 = product-subtracted A3 minus its known local product density.
 rawA3=a('2.3454197e-13')+a('6.704e-12')/(63*pi)
 E3=(rawA3+Hc*yc*sup_poly(O6)/(4*pi))*sm(4)/(p['B2']*p['Fl'])
 errors={'local_rest':a('1.875137e-42'),'higher_angular_UV':a('1.5483364e-43'),
  'higher_angular_positive_time':a('4.624e-45'),'low_angular_product_contact':a('2e-48'),
  'radial_residual':Erad,'entire_rational_remainder':Eapp,'local_A3':E3,
  'paired_boundary_reflection':Eboundary,'regular_bolt_buffer':Einner,
  'complete_thermal_midpoint_difference':Etherm,'cold_positive_time':coldpos,
  'cold_UV':colduv,'cold_product_conversion':coldproduct}
 total=sum(errors.values(),a(0));target=a('1.287e-40')
 c.NUMERICS.append(dict(kind='whole_signed_projection_certificate',N=N,
  P_hat_rest='0',signed_value_measured=False,complete_signed_interval=True,
  errors={k:enc(x) for k,x in errors.items()},E_rest=enc(total),target=enc(target),passes=bool(total<target),
  guards={k:enc(x) for k,x in {'eta':eta,'Ca':Ca,'Kerr':Kerr,'Ktail':Ktail,'proper_bulk_length':L,
   'Hrho':Hr,'buffer_length_lower':buffer,'actual_X_upper':xbound,'actual_Xv_upper':vbound,
   'inner_exponent_lower':amin,'cold_H':Hc,'raw_A3':rawA3}.items()},
  complete_angular_sums={str(k):enc(sm(k)) for k in [-2,-1,2,4,6,7,8,9,10,11]},
  rational_majorant_coefficients=[enc(x) for x in hm],seam_terms=bterms,
  actual_controlled_solution=False,rough_nonlinear_C2_action_proved=False))
 print(__import__('json').dumps(dict(event='whole_signed_projection_certificate',N=N,
  errors={k:str(x) for k,x in errors.items()},E_rest=str(total),passes=bool(total<target))),flush=True)
 return total
