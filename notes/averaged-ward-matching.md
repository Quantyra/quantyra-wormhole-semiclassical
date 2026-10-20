# Averaged Ward matching without a pointwise cut derivative

2026-10-19, E003/S150 supporting S139/S074. Ordinary argument and
verification protocol before execution.

A smooth overlap average determines the HH matching pressure at the
outer proper cut. Its second response can be bounded using scalar
source values in the two coordinate descriptions, without inserting
an estimate for S_chi into S148's pointwise chain rule. The result
requires an explicit common-family contract. It does not prove that
an arbitrary raw C5 metric family satisfies that contract.

## 1. Common family and units

Keep B,A,C,x,d, the physical period, paired domain and common scale
of S145-S149. Use S148's normalized canonical trace and pressure:

    S=2G Theta/(B^2 C), P=2G (B^2 p_hat)/(B^2 C).

On the fixed disk set E=exp(sigma), W=E^2 sinh^2 chi,
sigma=sigma0+b, and

    Y_b(chi)=integral_chi^D0 E(t)dt,
    a(Y_b(chi))=Omega E(chi)sinh chi,
    Omega=kappa A/sqrt(x).

Omega is fixed by the physical clock; it is not reset along a path.
Y increases inward. The same scalar pressure is used in both charts.

The following is the sufficient JOINT coordinate-family contract:

1. In the fixed disk b=b0+s h+t k is affine. The background and
   directions meet S145/S148 on the cover through seed Y0=.04,
   with background norm at most .001 and direction norms at most U,V.
   The global value bounds also hold. Thus
   |P|<6e-10, |P_h|<3.6e-7 U, |P_hk|<4.4e-6 UV.
2. In outer proper coordinates write a=a0 exp(-beta). The background
   meets S146/S147, including N(beta)<=1e-7 and .99<a<1.01.
   Its first coordinate directions obey N(beta_h)<=U and
   N(beta_k)<=V. Its mixed acceleration j=beta_hk obeys N(j)<=UV.
   All these are jets of the SAME geometry family and fixed domain.
   In particular |beta_h|,|beta_h,Y|<=100U, and similarly for k,j.
3. The actual canonical source at fixed Y meets the S146 source
   contract for those jets, including its actual-domain terms.
   Composition with j adds the first-source bound to the mixed
   bound. It is not omitted because b is affine in the other chart.
4. For outer constraint transfer, u=B/r obeys S147's value/first-jet
   guards. Its first relative value/derivative direction norms are
   at most U,V; the mixed acceleration has relative norm at most UV.
   The correlated mass/radius jets remain imposed.

These are conditional coefficient and source hypotheses. A C5 field
composed with its varying proper-coordinate map need not be a C2
path in C5. S149's flat gradient obstruction is not bypassed by
silently asserting otherwise. Compatibility and preservation of this
joint contract by the global solve remain unproved.

## 2. Exact smooth-average prescription

Choose fixed proper numbers

    ya=.041, yb=.049, y0=.05, L=yb-ya=.008.

Let z=(Y-ya)/L and define the C2, compactly supported density

    zeta(Y)=140 z^3(1-z)^3/L, ya<=Y<=yb,
    zeta(Y)=0 otherwise,
    Z(Y)=integral_ya^Y zeta(v)dv.

Use Z=0 below ya and Z=1 above yb. Its integral is one.
All its first two derivatives vanish at the support endpoints.
The moving support stays inside the S148 interior cover because
ya exp(-.001)>.04.

Define the pressure moment

    Q(b)=integral_ya^yb zeta(Y) a(Y)^2 P(Y)dY
        =integral_0^D0 G_b(chi) P_b(chi)dchi,
    G_b=zeta(Y_b) E_b a(Y_b)^2
       =Omega^2 zeta(Y_b) E_b^3 sinh^2 chi.          (1)

The canonical Ward equation is (a^2 P)_Y=a a_Y S.
Integrating it from Y to y0 and averaging proves exactly

    a(y0)^2 Pm=Q+I,
    I=integral_ya^y0 Z(Y) a(Y)a_Y(Y)S(Y)dY,
    Pm=P(y0).                                     (2)

The moment is calculated in fixed disk coordinates and I in fixed
proper coordinates. No pressure integration constant is assigned:
P_b inside Q is S148's bolt-regular HH pressure functional.
The nonzero pressure carried from the interior is retained.

The full canonical trace occurs in I. Its explicit cubic part is
not dropped. In the outer source budget write

    S=(2G/(B^2 C))R + 2eta_G c/(5pi Y).

The cubic and its first/mixed jets on [ya,y0] are negligible relative
to the following enlarged bounds, obtained from S147 and |c|<1,
|c_h|<.075U, |c_hk|<.075UV:

    |S|<7e-7 H(Y), |S_h|<3.5e-7 U H(Y),
    |S_hk|<.018 UV H(Y), H=1+log^2 Y.              (3)

For the mixed regular term the bound before rounding is
2eta_G(1e30+2e25)<.01700034. The second summand is the
coordinate acceleration contribution.

## 3. Fixed-disk transport kernel

At fixed chi,

    Y_h=integral_chi^D0 E h, Y_hk=integral_chi^D0 E hk,
    |Y_h|<=Y U, |Y_hk|<=Y UV.

With zeta and its derivatives evaluated at Y_b(chi), direct
differentiation of (1) gives

    G_h=E a^2[3h zeta+zeta' Y_h],
    G_hk=E a^2[9hk zeta+3zeta'(hY_k+kY_h)
                         +zeta''Y_hY_k+zeta'Y_hk]. (4)

The factor three includes both the orbit area factor and the
proper-distance Jacobian. No derivative of the quantum trace
appears. Because zeta is C2 and has compact support strictly
inside the cover, differentiation under this fixed integral is
justified by the uniform source bounds.

Unimodality gives integral |zeta'|=35/(8L)=546.875.
The derivative z^3(1-z)^3 is 3z^2(1-z)^2(1-2z),
of absolute value at most 3/16 and with one positive and one
negative lobe. Thus

    integral |zeta''|<=105/L^2=1640625.

Using a^2<1.03 and E dchi=-dY yields the kernel norms

    ||G||_1<1.03,
    ||G_h||_1<1.03[3+yb*546.875] U<31U,
    ||G_hk||_1
      <1.03[9+7yb*546.875+yb^2*1640625]UV<4300UV.

Consequently

    |Q|<6.18e-10, |Q_h|<4e-7 U, |Q_hk|<3e-5 UV. (5)

These are response bounds on the actual curved family under
section 1, not a flat multiplier substitution.

## 4. Proper-coordinate correction and pressure Hessian

Put w=a_Y/a and q_a=a a_Y=a^2 w. On [ya,y0],
|w|<.2. With h_o=beta_h, k_o=beta_k and j=beta_hk,

    (q_a)_h=a^2[-2h_o w-h_o,Y],
    (q_a)_hk=a^2[4h_o k_o w
                    +2(h_o k_o,Y+k_o h_o,Y)-2j w-j_Y].

The joint contract therefore gives

    |q_a|<.206, |(q_a)_h|<145U, |(q_a)_hk|<50000UV.

Also integral_ya^y0 H(Y)dY<.102. Equations (2),(3) give

    |I_hk|<.102[.206*.018+290*3.5e-7+50000*7e-7]UV
            <.004UV.                              (6)

S148 already bounds |Pm|<6e-10 and |Pm_h|<3.7e-7 U
using no spatial trace derivative. Differentiate (2) twice:

    Pm_hk=a(y0)^-2(Q_hk+I_hk)
            +2h_o Pm_k+2k_o Pm_h+(2j-4h_o k_o)Pm.

The right side is bounded by

    [1.03(.00003+.004)+400*3.7e-7+40200*6e-10]UV
       <.0044 UV.                                 (7)

This replaces the explicit missing S_chi in the matching-pressure
bound ON THE JOINT CONTRACT. It is not an unconditional C5
moving-cut theorem. If the source is only C2 in fixed disk values
but not in fixed proper values, hypothesis 3 can fail and (6)
cannot be invoked.

For example a continuous fixed-disk trace (chi-c0)_+^alpha,
0<alpha<1, with a cut crossing c0 can have a pressure whose first
derivative is only alpha-Hölder. Its moving-cut second response
need not exist. The proper-coordinate trace jets then fail the
contract. Smooth averaging alone does not fix that failure.

## 5. Outer constraint consequence, with acceleration retained

Allocate the outer canonical pressure as

    P=H_b Pm+L_b[R]+L_b[cubic],

with both forcing pieces chosen zero at y0. The first term now
carries the FULL interior pressure. This is an exact linear
allocation, not a second choice of state. The cubic forcing must
still be combined with its Weyl partner before the physical
absolute mass integral is bounded.

S147's memory coefficients, after converting from its unnormalized
pressure to P here, carry factors C/(2x)<25/8 for mass and
1/2 for lapse. Its affine/log mixed coefficient bounds acquire
at most their first bounds for the allowed coordinate acceleration:

    mass coefficients: (.05,.51,5.6+.51),
    lapse coefficients: (1e-5,1.4e-4,.002+.00014).

Using (7) and the sharper lower pressure jets gives

    |D^2 M_mem|<(25/8)[.05*.0044
                       +2*.51*3.7e-7+6.11*6e-10]UV<.0007UV,
    |D^2 A_mem|<(1/2)[1e-5*.0044
                       +2*1.4e-4*3.7e-7+.00214*6e-10]UV<2.3e-8UV.

For the regular trace-forcing component, direct first product
rules in S147 give |D M_reg|<1.1e-5 U and
|D A_reg|<4e-10 U. Thus its mixed bounds on the nonaffine proper
path are at most .039+1.1e-5 and 4.6e-7+4e-10.
Combining the regular forcing and full interior memory,

    |D^2(M_reg+M_mem)|<.04 UV,
    |D^2(A_reg+A_mem)|<4.9e-7 UV.                  (8)

The corresponding component Taylor remainders at segment norm
1e-7 are below 2e-16 and 2.45e-21, IF the joint contract
holds along the segment. These are not full physical error bounds:
the cubic/Weyl contribution, other source sectors, interior/cold
constraints and actual global inverse remain outside (8).

## 6. Verification and research decision

Before execution freeze exact checks of the averaging normalization,
Ward recovery, geometric transport, scalar normalization, coordinate
acceleration and all response budgets. Smooth manufactured traces
must satisfy the same Ward equation in both coordinates. They are
not quantum states or metrics satisfying the physical solution gates.

Use 32 and 64 Gauss nodes with 55-digit arithmetic for pressure and
overlap integrals; compare the average formula with independently
computed moving-cut pressure. Mixed finite differences at .001 and
.0005 must have error below 2e-6(1+|target|), and ratio below .35
unless both errors are below 1e-35. Quadrature refinement and direct
identity comparisons must agree to 1e-10(1+|target|). Omitting the
overlap transport, correction integral or coordinate acceleration
must be detected. Failed runs retain source snapshots.

The model and boundary-state distinctions remain those of
[Bintanja--Freivogel--Rolph](https://arxiv.org/html/2308.00871v2).
This averaging argument follows from the retained Ward identities;
the original construction is not treated as a proof of its norm
constants or of the global contract.

The smooth average closes matching memory as a conditional map
between the existing regional estimates. The next global proof
must either establish the joint coordinate/source contract or work
in a fixed-coordinate integrated norm that supplies equivalent
control. S149's endpoint result and S148's cold limitations remain
in force. No additional physical gate is completed.

Remaining to-do list: prove and preserve the global coordinate/source
contract and cold weighted response; complete cubic/Weyl, LLL,
reference and radius maps; solve the global exterior/state, charge
and signed source problem; opening, stability, finite-energy
transport, preparation and creation.
