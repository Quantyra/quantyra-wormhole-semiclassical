# Low tangent reconstruction and the regular divided gap equation

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, state/domain, full angular tower,
source subtraction, mass and fixed return. Improve the derivative
allowances feeding the complete scalar residual, using the existing
source orders and exact equations. No nonlinear solution or physical
stability conclusion is assumed.

## 1. Middle reconstruction starts at its actual low data

Use the complete middle source array Tm and the P,F arrays in the
[per-order reconstruction](finite-source-order-scalar.md). Their
normalization is j!100^j in canonical proper radius rho. The full
equations and analytic coefficient bounds are those of the
[middle proof](middle-metric-ten-jets.md).

The retained global bounds on both u,nu are

    |u|,|nu|<1.2e-8, |D(u,nu)|<4e-7,
    W|D(u,nu)|<6e-6, D=z0 nabla_gamma, W=t0^2.

At a real observation, canonical radial first derivatives therefore
are below min(4e-7 t,6e-6/t)<1.6e-6. Hence start the divided
radial arrays at U0=N0=1.2e-8 and U1=N1=1.6e-8. This uses
the proved pointwise gradient bounds, not an H1 embedding.

Apply the exact positive recurrences (7) of the middle proof for
n=0,...,8, instead of starting at n=4. Keep every P,F,Tm and
local term. Take coefficientwise minima with the old radial bounds
at each newly obtained order. Each recurrence uses only previous
metric indices. Thus it gives original tangent derivatives through
ten with a smaller low-order start and the same complete source.

## 2. Divide the cold equations before estimating derivatives

Use the [cold proof](cold-metric-ten-jets.md), with

    tau=1-z, p=sqrt(b), alpha=p/tau, L=p/z,
    beta=N'/N, S=-b_z/z+2b/z^2, nu=tau V.

Here tau>0 for the retained positive temperature, including at the
bolt. The following bounds introduce no inverse power of that
temperature. The exact identities beta p=-b_z/2 and tau''=-b_z/2
give

    tau''+(beta+2L)tau'=z S.

Writing f=F/tau, the exact equations become

    f'=-(L+alpha)f-2gamma alpha z(Q-T+E_rho),
    u''+beta u'=b_zz u-(z b_zzz tau/2)V
                                  +3F-3gamma z^2 T+gamma E_u,
    V''+(beta+2L+2alpha)V'=(2S/tau)u
                   +(c_nu-zS/tau)V-f+gamma z^2(T/tau)+gamma E_nu/tau.

The product-subtracted local quantities remain exactly
E_u=E_A-3z^2E_delta and E_nu=z^2E_delta. The canonical source
T and regular Ward pressure Q are the full retained ones. The
same cancellation in the radius equation is used; no physical
source or contact is discarded.

## 3. Uniform analytic coefficients and complete cold sources

On each .005 canonical proper disk, the old seed estimates give
|tau|<.004, |z| in (.996,1.004), |alpha|<3 and |L|<1. The
relative inverse-tau proof on the larger .05 disk gives coefficients
at most 1.2(.2)^j/tau0 in factorial/100 normalization. The
centered polynomial and delta_bar/tau0^2<1.002 imply

    |b_z|/|tau|<2.1, |S|/|tau|<4,
    |L+alpha|<4, |alpha z|<4,
    |2S/tau|<8, |c_nu-zS/tau|<16,
    |z b_zzz tau/2|<.04, |b_zz|<12.

For example delta_bar/|tau|<1.002*1.2*.001, and the explicit
b_z polynomial gives the first estimate. Then use |p/tau|<3,
|tau|<.004 and |z|>.996 in S/tau. Each displayed constant
has normalized derivative majorant equal to that constant times
2^j, by Cauchy on the .005 disk.

The cold proof supplies |E_rho|,|E_delta|,|E_A|<1e22 s0 on
each such disk, before replacing s0 by its larger logarithmic
weight a(s0)=s0[1+log(.01/s0)]. Therefore

    |E_rho|<1e19,
    |E_u|<4.1e19,
    |E_nu/tau|<1.3e22.

The last two inequalities use 1+3*1.004^2<4.1 and
1.004^2*1.2*(s0/tau0)<1.3, with s0/tau0=1+s0<=1.001.
Their normalized derivative reserves are these constants times
2^j. This differentiates actual holomorphic local formulas, not
an unsigned value inequality. The nonlocal T is differentiated
only using its existing real actual-domain source theorem.

Let Tc_j,Q_j and Fa_j be the frozen complete source, Ward and
weighted mass arrays. Set A=.00331, g0=6.8e-28 and

    T_j^u=A Tc_j, Q_j^u=A Q_j,
    T_j^v=1.2*82 sum_(i=0)^j (.2)^i Tc_(j-i).

These bound T,Q and T/tau respectively. The weight bounds
a<=A and a/tau<82 are retained. Define C_n(c,X)=sum c2^iX_(n-i)
and D_n(c,X)=sum c2^i100(n-i+1)X_(n-i+1). With f0=82e-10,

    f_(n+1)={C_n(4,f)+2g0 C_n(4,Q^u+T^u+E)}/[100(n+1)],
    E_j=1e19*2^j.

Use n=0,...,7. Bound F either by A Fa_j or by the product tau f,
taking the minimum. For the latter use tau0<=.001 and normalized
positive derivatives .004*2^j. Every mass and Ward term remains.

## 4. Regular coupled radial recurrence through the bolt

On rho<=1, beta=1/rho+c has regular analytic c; on rho>=.5
the full beta is regular. The old common coefficient allowance
12*2^j covers either drift. The V drift has the larger allowance
20*2^j after adding 2L+2alpha. The regular flat radial inverse is

    X^(n+2)=G^(n)-integral_0^1 t^(n+1)G^(n)(t rho)dt

for X''+X'/rho=G. Its factor-two bound uses only the nth
source derivative, extends through rho=0, and also covers the
smaller away estimate. All new cold arrays are uniform, so the
interior segment is covered by the same bounds.

Start unweighted radial arrays at

    U0=V0=1.2e-8, U1=4.1e-9, V1=3.4e-9.

The u first derivative follows from |Du|<4e-7 and t<1.001;
the V derivative is the retained canonical cold bound 3.4e-7.
For n=0,...,8 use

    U_(n+2)=2{D_n(12,U)+C_n(12,U)+C_n(.04,V)+3F_n
                         +3g0 C_n(2,T^u)+g0*4.1e19*2^n}/den,
    V_(n+2)=2{D_n(20,V)+C_n(8,U)+C_n(16,V)+f_n
                         +g0 C_n(2,T^v)+g0*1.3e22*2^n}/den,
    den=1e4(n+1)(n+2).

Take minima with all earlier applicable radial estimates, including
the cold Hessian bound for V at order two. Recover nu by the
ordinary product tau V and take its minimum with the old weighted
radial estimate enlarged by a<=A. These are the original u,nu,V;
there is no new metric or boundary choice.

## 5. Full tensors retain the fixed-observation powers of z

Apply the existing middle, cold-Cartesian and cold-away tensor
recurrences to the nonconstant radial arrays through ten. Keep
the frozen outer arrays, including analytic mass remainder and
correlated cubic. The cold scalar arrays use the regular Taylor
argument through the center and every ordered connection slot.

The tensor recurrence gives canonical nabla_gamma jets. At a fixed
observation D^j=z0^j nabla_gamma^j, exactly: z0 is not differentiated.
Consequently W D^j has factor t0^(2-j). On the middle region a
safe factor is 50.1^(2-j) for j=0,1 and 1 for j>=2. Retaining
this factor avoids multiplying every high derivative by 50.1^2.

For middle V=nu/tau use the .05-disk inverse coefficients. Its
unweighted radial allowance remains 1200 sum(.2)^i N_(j-i),
as proved in the middle note. For the weighted tensor bound first
retain the relative coefficient 1.2/tau0, then use

    W z0^j/tau0=t0^(3-j)/(t0-1)
       <=2600 (j=0), <=1003 (j=1), <=1001 (j>=2).

For j=0,1 the only interior critical point is a minimum; the
endpoint bounds suffice. For j>=2 the function decreases with t.
For W tau0 D^jV the factor is just 1.2 t0^(2-j).
On the cold region use W<1.001^2 and tau0<.001; this preserves
all endpoint weights without an inverse temperature. Take minima
with the earlier global U,K,Vbar,KV,h arrays, retaining their
independent small values, gradients and gap Hessian.

## 6. Complete scalar consequence and independent checks

Feed these new tangent arrays into the unchanged six-step finite
geometry refinement and complete scalar/anomaly formulas of
[the preceding audit](finite-anomaly-tangent-geometry.md). Every
curvature product, metric contraction, magnetic term, seed monomial,
reference/Lie term, actual covariant conversion and noncompact
volume factor remains. Report the full C_tau, its upward enclosure
and the scalar-only transfer budget before reassessing the next step.

Freeze note, checker and review before execution. Verify coefficient
and local bounds by exact rational inequalities, raw Leibniz forms
of the radial recurrences, all full tensor conversions and retained
inputs. Independently verify the divided equations on nonconstant
regular manufactured profiles with tau=h cosh(rho), including the
exact bolt and two very small positive h values. Check derivatives
through eight at 60/90 digits. Also check the regular radial inverse
through ten on nonconstant even polynomial profiles, retaining its
interior integral. Manufactured tests do not sample the charged
quantum source, and supporting checks are not formal certification.

Remaining to-do list: assess useful complete source and transfer
bounds; quantum/TT H4 and higher lift/inverse; nonlinear mapping,
preserved neighborhood, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
