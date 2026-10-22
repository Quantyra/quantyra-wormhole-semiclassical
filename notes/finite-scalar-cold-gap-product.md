# Gap Hessian and cold-factor scalar products

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, complete source allocation,
physical state/domain, magnetic flux, mass and return. Improve the
two dominant terms of [the refined scalar bound](finite-scalar-residual-six-refined.md).
The new exact rational ledger defines C_tau in the same estimate

    ||f_tau||_(H_W6)/sqrt(gamma0)<=C_tau ell^2/g.

No transfer, nonlinear existence or stability conclusion is assumed.
The source remains the bookkeeping action gradient, with its exact
quadratic path zero and the full noncompact physical volume.

## 1. A small gap Hessian through the bolt

Use the seed canonical proper radius rho, b=calB/C, p=sqrt(b),
tau=1-z, nu=tau V, beta=N'/N, alpha=p/tau and L=p/z. Primes in
this section denote rho derivatives. The original radius equation
from [the cold reconstruction](cold-metric-ten-jets.md) is

    nu''+(beta+2L)nu'=2S u+c_nu nu-F+gamma z^2 T+gamma E_nu,
    S=-b_z/z+2b/z^2, E_nu=z^2 E_delta.

Here T is the complete canonical product-subtracted trace. Use its
retained small VALUE bound, not the enlarged C8 reserve. The
[repaired source dependency](repaired-source-metric-transfer.md)
retains the cold 1e20 a(s) trace bound; including the explicit LLL
difference fits |T|<2e20 a. The local real value bound is
|E_delta|<1e22 a, hence |E_nu|<1e22 a since z<=1. The larger
1e24 allowance for differentiated complex local coefficients is
unnecessary for this zeroth-order equation.

The other inputs are gamma<6.8e-28, |F|<1e-10 a,
a/tau<82, |u|,|V|<1.2e-8 and |V'|<3.4e-7 on the cold patch.
The last bound follows from |V_y|<8e-7 and
d rho/dy=sqrt(C/x)/t, C/x>5.9, t<1.001.

The centered seed polynomial gives

    .99<b_zz/2<1.1, |b_z|/tau<2.1,
    alpha<1.01, L<.0011, |S|/tau<3, |c_nu|<10,
    |tau''+(beta+2L)tau'|/tau<3.

For example b=tau^2-B3 tau^3+B4 tau^4-4e tau^5+e tau^6
-delta_bar(1-tau)^3, with |B3|<1.334, |B4|<.501,
|e|<1e-24 and delta_bar<1.001h^2. The lower bound on b_zz/2
is greater than 1-3*1.334*.001 minus the tiny remaining terms.
The same polynomial bounds |b_z|/tau. Also p/tau<1.01 follows
globally from b<=tau^2 g0 and
g0<=6/5.9+1e-10+1e-24<1.01^2. No inverse thermal-gap constant
is introduced.

Dividing the exact nu equation by tau gives

    V''+beta V'=G,
    G=-(2L+2alpha)V'
      +[c_nu-(tau''+(beta+2L)tau')/tau]V
      +2S u/tau-F/tau+gamma z^2 T/tau+gamma E_nu/tau.

Its absolute cold bound is the following rational sum:

    |G| <=(2*.0011+2*1.01)*3.4e-7+13*1.2e-8
           +6*1.2e-8+82e-10
           +6.8e-28*82*(2e20+1e22) < .0006.

Every product subtraction occurs before this division. The source
and local terms retain a(s); a(s)/tau<82 removes the apparent loss.

Regularity fixes V'(0)=0 and therefore

    V'(rho)=N(rho)^-1 integral_0^rho N(s)G(s)ds.

The comparison for N''=kN, .99<k<1.1, N(0)=0,N'(0)=1,
gives .99<beta<1.1+1/rho. For the lower bound, differentiate
N'^2-.99N^2 and use N'>0. For the upper bound, compare the
Riccati equation with sqrt(1.1)coth(sqrt(1.1)rho), whose value is
below 1.1+1/rho. If rho<=1, monotonicity of N bounds
beta integral N/N by 2.1. If rho>=1, the lower beta bound gives
integral N/N<1/.99 and beta<2.1. Thus the product is below 2.13
on the entire cold interval. Consequently

    |beta V'|<2.13*.0006, |V''|<3.13*.0006,
    |Hess_gamma V|<5.26*.0006<.0032.

These are the full two-dimensional Hessian slots, including their
regular limits at the bolt. With D=z0 nabla_gamma and W=t0^2,
W|D^2V|=|Hess_gamma V|. A safe cold reserve is .004.

## 2. Global low values, gradients and Hessian

Keep the [corrected low physical equations](finite-trace-gradient-and-projection-transversality.md).
For y=sqrt(lambda)rho_physical they give
tau|V_yy|<1.1e-5, |A0 nu_y|<3.432e-6, A0 R0<=11,
R0<2.24 and |V_y|<8e-7. Since tau_y=z R0,

    tau|A0 V_y|<3.432e-6+11*1.2e-8=3.564e-6.

The canonical Hessian slots are exactly

    D^2V_radial=(V_yy+R0 V_y)/(C/x),
    D^2V_angular=(A0-R0)V_y/(C/x).

On 1.001<=t<=50.1, the maximum of t^2/tau=t^3/(t-1)
occurs at an endpoint and is below 2600. Thus

    W|D^2V| < [2600*(1.1e-5+3.564e-6)
                        +50.1^2*2*2.24*8e-7]/5.9 < .008.

Without W, tau>=.001/1.001 gives |D^2V|<.0025.
On t>=40, retain the independent bounds
t^2|V_y|<4.3e-8, t^2|V_yy|<2.1e-6 and A0<2;
the same formula gives W|D^2V|<1e-6. The covers overlap. Hence

    |D^2V|<.004, W|D^2V|<.01.                 (1)

The independent value bounds on t<=4 and t>=4 also give
W|u|,W|nu|,W|V|<2e-7. For gradients use |u_y|<3.2e-7/t
up to t=40 and the old outer tensor bound W|Du|<3e-6 beyond it.
For nu,V use the core bound up to four and
nu_y=tau V_y+tau_y V with the outer |V|<2e-9/t^2 thereafter.
Together these imply

    W|Du|,W|Dnu|,W|DV|<6e-6.                  (2)

No derivative is inferred by differentiating a value inequality.
In the previous divided-jet arrays, one may therefore replace
U0,Vbar0 by 1.2e-8, K0,KV0 by 2e-7, K1,KV1 by 6e-6,
and Vbar2,KV2 by .002,.005 respectively, taking minima with
any stronger entries. The other established entries remain valid.

There is a further useful scalar value bound

    W|L4nu|/tau<.005.                         (3)

In the cold region use the undifferentiated nu equation above:
the RHS divided by tau is at most
6*1.2e-8+10*1.2e-8+82e-10+6.8e-28*82*(2e20+1e22).
Since L4nu=(C/x)z^2[nu''+(beta+2L)nu'], multiplying by W
costs at most 6.1. The result is below .004. In the middle use
the corrected |L4nu|<6.3e-7 and t^2/tau<2600. In the outer
region, differentiate nu=tau V twice, using the retained
|tau_y|<=2.24/t, |tau_yy|/tau<21 and the actual outer derivatives
above. Then W|nu_yy|<2.2e-6, W|nu_y|<5e-8, so
W|L4nu|/tau<[2.2e-6+6.48*5e-8]/.975<3e-6.

## 3. Retain tau on high gap derivatives

Let v_j be the improved W-weighted divided D-jet array KV.
Define h_j bounding W tau|D^jV|/j!, separately from v_j.
The relative inverse-tau disk bound is 1.2/tau at each real
observation, both in the middle and cold charts. The middle bound
1200 used previously was its uniform enlargement.

Let r_j=1.2 sum_(i=0)^j (.2)^i N_(j-i) in normalized radial
units j!100^j. Apply the same full tensor recurrence as before:

* Middle: use r_j, then multiply by 50.1^2.
* Cold away: use r_j and the away recurrence, then multiply by
  .00331*1.001^2, retaining a(s) until multiplication by tau.
* Cold near: use 82r_j and the Cartesian/connection recurrence,
  then multiply by .001*1.001^2. Here the radial Taylor integral
  uses a(s)/tau<82 along the segment, while the observation has
  tau<.001. It is not assumed that a(s)/tau increases along it.
* Outer: the previous weighted V array also bounds tau times it.

Take the maximum over the cover, multiply by 100^j, and take the
minimum with v_j. This proves h_j without dividing by a lower
thermal gap. Also |D tau|<=1.01 tau, as shown in section 1.

## 4. Two weighted product bounds

Let q=1+Zeta majorize the ordinary tau jets, with q0=1, and
define a_i=(i+1)v_(i+1), b_i=(i+1)h_(i+1). For n<=6 put

    G_n=sum_(k+i+j=n) g_(k,i,j),
    g_(0,i,j)=min(a_i b_j,b_i a_j),
    g_(1,i,j)=1.01 min(a_i b_j,b_i a_j),
    g_(k,i,j)=q_k a_i a_j for k>=2.

Then G_n bounds W^2|D^n[tau |DV|^2]|/n!. This is the full
ordered tensor Leibniz rule: contraction with gamma is parallel,
and divided derivatives absorb the multinomial coefficients.
For k=0 place tau on whichever factor has the better bound; for
k=1 retain it using |D tau|/tau<=1.01. Terms k>=2 contain only
V derivatives through five and keep the previous ordinary tau jets.
No tau derivative or tensor slot is discarded.

Let LN be the previous W L4nu array, with its zeroth entry reduced
to .005 using (3). Define

    H_n=sum_(i=0)^n v_i LN_(n-i),

except that the i=n term is bounded by
min(v_n LN0,.005 h_n). This bounds W^2 D^n[V L4nu]/n! and
uses (3) only on its undifferentiated L4nu factor.

With O,E and the other arrays defined in the preceding exact
remainder proof, replace its third and fifth terms by

    N3=6*1.001 O q E H,
    N5=6*6.1*1.001 O E q^2 G.

Recompute the other four terms with the stronger low entries from
(1)-(2). Keep the old finite-field S array, anomaly, reference and
Lie terms. Take coefficientwise minima with the previous valid
curvature remainder. The same actual-connection conversion and
physical-volume integral then produce the new scalar H6 constant.
Report its size and the scalar-only transfer budget before deciding
whether this direct estimate has become useful.

## Verification scope

Freeze this argument, checker and review before execution. Check
the exact divided radius equation, rational source/seed margins,
global Hessian cover, retained regional h_j arrays and every product
index. Compare the regular bolt integral with an independent closed
solution for nonconstant forcing at 60/90 digits. Independently
differentiate a two-dimensional Cartesian polynomial gradient
square through six, including tau derivatives and all tensor slots.
These checks support an ordinary proof, not formal certification or
samples of the actual charged quantum source. Preserve all runs.

Remaining to-do list: assess useful scalar and quantum/TT source
bounds; higher lift/inverse and positive transfer; nonlinear mapping,
preserved neighborhood, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
