# Fixed-chart global constraint transfer

2026-10-20, E003/S152 supporting S139/S074. Ordinary argument and
verification protocol before execution.

The complete mass/lapse integration map can be differentiated in the
fixed common disk. It does not require a separately differentiated
pressure at a moving proper-coordinate cut. This avoids S151's
coordinate regularity cost for these integrated outputs. The actual
quantum-source map must still be proved in the weighted input space
below, and the full higher-order equation inverse remains open.

This is a global SOURCE-TO-CONSTRAINT theorem with quantitative
finite-background geometry factors. It is not a new existence or
stability theorem. It extends the particular endpoint calculation
in [the complete residual](stationary-centered-complete-residual.md)
to arbitrary families satisfying explicit low geometry and source
norms. That earlier endpoint defect is not promoted to a solution.

## 1. Fixed geometry and exact densities

Retain B,A,C=B^2/A^2,x=lambda_e B^2, h, the physical period,
boundary clock, paired domain and subtraction scheme. Set

    t=r0/B, s=t-1, v=1+h, tau=s/t,
    eta_G=GC/B^2<8.5e-33, C/x<6.1,
    h<1e-30, L_h=1+log(.01/h)<81.

The seed r0(chi) increases from B(1+h) to infinity on the fixed
disk 0<=chi<D0. Use arbitrary real fields b,V with

    sigma=sigma0+b, r=B+(r0-B)exp(V),
    z=s partial_t V, T=r/r0=1+tau(exp(V)-1),
    max(|b|,|V|,|z|)<=.001.                         (1)

The fields satisfy the chosen bolt, clock and asymptotic boundary
conditions. The estimate below does not replace those conditions.
In particular z is a derivative of V, not an independent radius
degree of freedom. At the bolt it is its regular limiting value.
The radius is monotone since 1+z>0.

Let rho,p be the combined sphere-integrated physical source whose
mass and lapse constraints are being evaluated. Write

    H_rho=r^2 rho, H_N=r^2(rho+p).

The retained canonical, Weyl, LLL and reference/radius contributions
must be combined in the prescribed common scheme before imposing
the outer norm. A single unremoved canonical mass cubic does not
belong to that norm. This theorem also applies to an already
integrable source difference; it does not invent a new subtraction.

The exact relation F=C exp(-2sigma)(r_chi/r)^2 gives

    r_chi/r0,chi=exp(V)(1+z),
    F/F0=exp(2V-2b)(1+z)^2/T^2,
    G rho r_chi=G J_M H_rho r0,chi/r0^2,
    G(rho+p)r_chi/(rF)=G J_A H_N r0,chi/(r0^3 F0),
    J_M=exp(V)(1+z)/T^2,
    J_A=exp(2b-V)/[T(1+z)].                       (2)

These identities are geometric; the metric need not solve the
constraint. The changing radius/frame factor is counted exactly
once in (2). Both sources are scalar profiles pulled to fixed chi.

Use affine common-chart parameters for b,V and therefore z.
For a direction h=(h_b,h_V,h_z), let
U=max(||h_b||inf,||h_V||inf,||h_z||inf), and similarly K for k.
Put a=tau exp(V)/T=(r-B)/r, so 0<=a<1. Exact logarithmic
derivatives of the multipliers are

    l_M,h=(1-2a)h_V+h_z/(1+z),
    l_A,h=2h_b-(1+a)h_V-h_z/(1+z),
    l_M,hk=-2a(1-a)h_V k_V-h_z k_z/(1+z)^2,
    l_A,hk=-a(1-a)h_V k_V+h_z k_z/(1+z)^2.        (3)

Then J_h=J l_h and J_hk=J(l_h l_k+l_hk).
The following bounds hold everywhere, independently of h's small
thermal-gap value (h in that phrase is the fixed gap, not a direction):

| Multiplier | Value | First / U | Mixed / UK |
| --- | ---: | ---: | ---: |
| J_M | 1.003 | 2.01 | 5.54 |
| J_A | 1.006 | 5.04 | 26.5 |

For the value bounds use
|log J_M|<=.001-log(.999) and
|log J_A|<=.004-log(.999).
Here |1-2a|<=1, a(1-a)<=1/4,
1/(1+z)<=1/.999 and 1/(1+z)^2<=1/.999^2.
Inserting these in (3) proves every entry of the table. In
particular a nonzero background radius does not freeze the
denominator or its second variation.

## 2. Global input space with four regional weights

All regions and weights are FIXED BY THE SEED. There is no
parameter-dependent division into charts. Define

    s_c=.001, t_c=1.001, t_e=40,
    L_c(s)=1+log(.01/s), L_e(t)=1+log t.

For a profile pair H=(H_rho,H_N), use four bounds:

| Region | Required profile bound, after division by C^2 |
| --- | --- |
| cold, h<=s<=s_c | |H_rho|/C^2 <= Z_c s L_c(s) |
| cold away from bolt, chi>=1 and s<=s_c | |H_N|/C^2 <= Z_c s L_c(s) |
| bolt neighborhood, 0<=chi<=1 | |H_N|/C^2 <= Z_b h L_h sqrt(chi) |
| middle, t_c<=t<=40 | max(|H_rho|,|H_N|)/C^2 <= Z_m |
| outer, t>=40 | max(|H_rho|,|H_N|)/C^2 <= Z_e L_e(t)^2 |

The bolt neighborhood lies in the cold region: the retained seed
has s(chi)<=2h for chi<=1. Values at the two interfaces may be
assigned both bounds; their measure-zero overlap is immaterial.
These weighted sup norms, or their measurable completions, define
the input space. The lapse condition includes H_N(0)=0.

For j=0,1,2 write Z_{i,j}, i=c,b,m,e, for bounds on the profile,
its first parameter response divided by U, and mixed response
divided by UK. The hypothesis is C2 dependence INTO THIS SPACE,
not merely pointwise C2 dependence with an unspecified majorant.
No spatial derivative of a varying source is assumed.

The four constants must be retained separately in applications.
An outer logarithmic source bound must not be substituted into
the cold column; doing so can destroy a useful quantitative result.

## 3. Fixed linear integration kernels

For a mass density f_M and lapse density f_A on the common disk,
the output norms are

    Q_M=sup_(chi>0) |integral_0^chi f_M|/[r0 F0],
    Q_inf=(BC)^(-1) integral_0^D0 |f_M|,
    Q_A=integral_0^D0 |f_A|.                       (4)

Here f_M=G H_rho r0,chi/r0^2 and
f_A=G H_N r0,chi/(r0^3 F0) for the linear kernels.
The first output is a zero-horizon-datum mass primitive. The
second controls absolute geometric mass density; it is not a
renormalized Hamiltonian charge. The third is an L1 density norm.
We differentiate the maps into these Banach spaces, NOT the
absolute-value or supremum norm as a scalar function.

The following matrix bounds the linear map, with overall factor
eta_G and column order (Z_c,Z_b,Z_m,Z_e):

                    cold       bolt    middle     outer
    Q_M          :    46         0       279       .00016
    Q_inf        : .000002       0         1       .835
    Q_A          : 13122        486      1300      .00003.        (5)

Here is a direct derivation using retained seed guards. They follow
from the exact seed in
[the centered comparison](stationary-centered-covariant-seed.md)
and the regional estimates in the complete residual:

    F0>=.9C(s^2-h^2)/t^2                  (cold),
    F0>=.9C s^2/t^2                       (s>=.001),
    F0>=.5x t^2                          (t>=40),
    C dt/F0<=3 dchi/(h chi)              (chi<=1),
    C ds/F0<=4 ds/s^2                    (chi>=1, cold),
    C integral_tc^infinity dt/F0<1300.                (6)

The second lower bound also follows directly by factoring the
seed F_ext: its square-bracket divided by C is at least
1-O(1e-24) for t>=1. The subtracted thermal term divided by
C s^2/t^2 is below 2h^2 t/s^2<1e-53 on s>=.001.
No extrapolation of the cold quantum estimate is involved.

For Q_M in the cold region,

    integral_h^s y L_c(y)/(1+y)^2 dy
        <=L_h(s^2-h^2)/2,
    C integral_h^s .../[t F0(t)] <=1.001 L_h/1.8<46.

The denominator t F0(t) increases with t, so the completed cold
contribution cannot grow farther outward. The middle part costs

    (s-s_c)/[.9 t_c s^2] <=1/[3.6 t_c s_c]<279.

For an outer-only source, use

    I_e=integral_40^infinity L_e(t)^2/t^2 dt
       =[(log40)^2+4log40+5]/40<.835,
    (2C/x) I_e/40^3 <.00016.

This proves the first row, keeping the cold anchor at its own
denominator rather than dividing total mass by a tiny cold gap.
For Q_inf the cold integral is below

    integral_0^.001 s L_c(s)ds
       =(.001)^2[1.5+log10]/2 <.000002,

the middle integral of 1/t^2 is below one, and the outer is I_e.
For Q_A the bolt part is at most
3L_h integral_0^1 chi^(-1/2)dchi=6L_h<486.
The remaining cold part is at most

    4 integral_h^.001 L_c(s) ds/s <2L_h^2<13122.

The middle is below 1300. The outer is below

    (2C/x) integral_40^infinity L_e(t)^2/t^5 dt
    <12.2/40^4 [(1+log40)^2/4+(1+log40)/8+1/32]
    <.00003.

This proves (5). The uniform constants use the retained finite
preparation range L_h<81. They are not a zero-temperature limit:
the unrounded cold bounds still grow logarithmically as h->0.

If all four input bounds are replaced by a common Z, the three
rows are below (330,2,15000) eta_G Z. The regional matrix is
the useful result; this common-Z envelope is only a convenience.

## 4. Complete finite-background first and mixed maps

Insert J_M H_rho and J_A H_N into the fixed kernels.
In each region i define

    B_M,i,0=1.003 Z_i,0,
    B_M,i,1=1.003 Z_i,1+2.01 Z_i,0,
    B_M,i,2=1.003 Z_i,2+4.02 Z_i,1+5.54 Z_i,0,
    B_A,i,0=1.006 Z_i,0,
    B_A,i,1=1.006 Z_i,1+5.04 Z_i,0,
    B_A,i,2=1.006 Z_i,2+10.08 Z_i,1+26.5 Z_i,0.    (7)

For j=0,1,2, multiply the first two rows of (5) by the
B_M,i,j and the last by B_A,i,j, with overall eta_G. These
are value, first and mixed output bounds, with the usual
direction factors 1,U,UK. Every source-Jacobian cross term
is included.

The map is C2 because the scalar multipliers and their first
two derivatives are uniformly bounded in L-infinity on (1),
the source is C2 in the weighted space, and (5) gives bounded
linear maps to the output spaces. This is also a dominated
convergence proof at the cold bolt and infinity. Only low
geometry jets appear in (2)-(3); no moving C5-to-C5 coordinate
map is differentiated.

For a one-parameter affine common-chart path preserving all
hypotheses, one half of the j=2 bound times its parameter
length squared bounds the component Taylor remainder. The
numerical source constants Z_i,j have NOT been proved for
the unknown full quantum family, so (7) is not assigned an
unconditional physical numerical error.

## 5. Why an artificial matching cut disappears

Let c(b,V) be any moving interior cut, and split an actual
fixed-chart density f as

    integral_0^D0 f = integral_0^c f + integral_c^D0 f.        (8)

For smooth data each mixed derivative of the first piece has
boundary terms

    f_h(c)c_k+f_k(c)c_h+f_chi(c)c_h c_k+f(c)c_hk.

The second piece has their negatives. Thus all cut transport
and acceleration terms cancel for the complete integral.
The right proof for nonsmooth source values is even simpler:
combine (8) BEFORE differentiation and apply the bounded
fixed-domain map. The individual pieces need not be C2.

For example f(x)=sqrt(max(x,0)) on [-1,1] and c=s+t gives
the lower piece (2/3)(s+t)_+^(3/2), whose mixed centered
quotient diverges as the step tends to zero. The upper piece
has the opposite divergence; their sum is exactly 2/3.
This example is a regularity diagnostic, not a physical state.

This cancellation requires the SAME physical source on both
sides. In the canonical sector S148's HH Ward prescription
fixes a single pressure functional, hence its boundary value
is not an independently adjustable constant. If inconsistent
pressures or subtraction allocations are supplied on the two
sides, (8) is not their combined source integral and no such
cancellation follows.

S150 remains a valid estimate for a separately exposed proper-cut
pressure under its joint contract. S151 remains a valid stronger
coordinate theorem. Neither is logically required merely to
differentiate the complete integrals once the source estimates
in this fixed-chart space are available. This observation does
NOT transfer S146's proper-coordinate source theorem into this
space automatically, nor does it provide the missing cold
oscillation response. Those source tasks remain substantive.

## 6. Limits and verification protocol

The norm needs both cold gap scaling and a bolt modulus. In
the elementary model r0=B(1+h cosh chi), sigma0=0,

    F0=C[h sinh chi/(1+h cosh chi)]^2.

For H_N=C^2 h phi(chi), its lapse density is exactly

    eta_G phi(chi)/[(1+h cosh chi)sinh chi].

Thus phi=chi^2/(1+chi^2) is integrable uniformly in small h.
Dropping the factor h loses that uniformity; taking phi=1
leaves a logarithmic bolt divergence even with the factor h.
Source-value continuity alone must not be substituted for the
weighted hypotheses. This model tests the norm, not the actual
exterior or state.

Before execution, verify exact multiplier/variation identities,
all regional integrals and constants, and interface cancellation.
Use manufactured even b,V and independent profile jets in the
preceding model, on 0<=chi<=2.5. Compare direct physical
constraint integrals with the fixed-chart forms and independent
centered differences. Use 55 digits, 32/64 high-precision Gauss
nodes, steps .001,.0005, mixed tolerance 2e-6(1+|target|),
refinement ratio .35 unless both errors are below 1e-30,
and identity/refinement tolerance 1e-11(1+|target|).
Test h=1e-5 and 1e-30 and the analytic failure diagnostics.
Retain all runs, source snapshots for failures, and exact
input hashes before and after.

The retained physical model and ensemble are those discussed by
[Bintanja, Freivogel and Rolph](https://arxiv.org/abs/2308.00871v2).
That work does not supply the new profile-space bounds or a
full higher-order inverse here.

The next source proof should target the four Z_i,j directly in
the common disk, with canonical/Weyl recombination before the
outer absolute norm and a quantitative cold modulus. The full
inverse must additionally recover the required metric regularity
and enforce the horizon, clock and state conditions. Equation
(5) bounds only low constraint outputs; it is not that inverse.

Remaining to-do list: actual global source responses in this space
and a usable full inverse; remaining source sectors and coupled
exterior/state, charge and signed-source solve; opening, stability,
finite-energy transport, preparation and separate creation.
