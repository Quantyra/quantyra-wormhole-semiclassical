# Cold response from real endpoint estimates

2026-10-20, E003/S153 supporting S139/S074. Ordinary proof and
verification protocol before execution.

The actual canonical Dirac source has finite-background cold first
and mixed response bounds on the intersection of the established
real endpoint and fixed-domain coefficient contracts. A real
endpoint estimate is NOT simply subjected to complex Cauchy.
Instead an elementary smallness theorem transfers it to derivatives
with logarithmic losses. The resulting quarter-power bolt modulus
still controls absolute lapse integration.

This supplies a cold canonical component for
[S152](fixed-chart-global-constraint-transfer.md). The constants are
large. This is not a contraction theorem, an unconditional C5
theorem, or proof that the previously constructed endpoint and
unknown solver preserve the required metric contracts.

## 1. Admissible common-chart families and retained inputs

Use the fixed physical period, B,A,C,q, mu0, paired actual domain
and common scheme of [S145](general-interior-dirac-response.md).
Write h for the positive thermal gap, s=r0/B-1, and

    eta_G=GC/B^2<8.5e-33, q/C<400,
    exp(-80)<h<1e-30, L_h=1+log(.01/h)<81,
    a_x=s[1+log(.01/s)], h<=s<=.001,
    a_x<.00331, a_x<=2h L_h on chi<=1.              (1)

The increasing weight a_x and the regions are seed quantities.
Let N be a real norm on common conformal fields, with E=1e-6.
The family must have the following admissibility properties:

1. The actual global/cap value and gradient norms and fixed
   interior raw C5 chart norms of S145 are bounded by N.
   These are required uniformly over the cold observation cover.
   The same linear representation has the complex coefficient
   bounds needed for its parameter holomorphy.
2. For each REAL endpoint N(v)<=E, the seed-to-endpoint comparison
   satisfies the sufficient cold, cap, global-domain and compatible
   exterior hypotheses of
   [the real canonical response](stationary-centered-canonical-response.md),
   with cold size at most E. Its endpoint coordinate identifications
   are made once, as in that proof. They are not differentiated in
   the family parameter or evaluated at a complex spatial point.
3. In the common disk, directions obey
   |u|<=N(u)a_x and |u_chi|<=8N(u)a_x in the cold region.
   Boundary/bolt regularity, the same actual paired domain and the
   physical clock are preserved.

These are explicit intersected geometry hypotheses. In particular,
condition 2 includes more than a pointwise C5 value bound. We do
not prove here that an arbitrary raw C5 ball or the computed endpoint
satisfies it. It is a contract for applying the existing real
endpoint theorem to every real member, not an assumed bound on
the new quantum-source derivatives.

Take a real background b with N(b)<=E/2 and directions u,w
of norms U,V. Let

    S[b]=(2G/(B^2 C))Theta[b]

be the Lorentzian canonical scalar trace. Its sign is opposite
to S145's action source, with identical absolute bounds.
The real endpoint estimates in the cited proof, equations
(16),(18), imply for every admissible REAL v with N(v)<=E,

    |S[v](x)-S[0](x)|<=1.7e-15 a_x,
    |delta S[v](chi)-delta S[v](0)|
        <=1.7e-14 h L_h sqrt(chi), chi<=1,          (2)

where delta S[v]=S[v]-S[0]. The factors are
2 eta_G * 1e23 E and 2 eta_G * 1e24 E, respectively.
This uses the whole angular tower and the actual global domain.
No endpoint difference is relabeled as a second response.

There is also a DIFFERENT, unweighted complex estimate. S145
bounds the nonconstant normalized scalar on its parameter disks
by K<5e-12, after separating the fixed flat-vacuum scalar.
For N(d)<=1 consider b+(E/2) z d, |z|<=1. Its coefficient size
is at most E, well inside S145's .011 complex cone. Therefore

    f_x(z)=S[b+(E/2)zd](x)-S[0](x),
    g_chi(z)=f_chi(z)-f_0(z)

are holomorphic and each is bounded by

    B_c=2e-11.                                     (3)

For g there are four nonconstant scalar values, hence 4K.
For f two suffice. The fixed vacuum value cancels; it is not
set to zero. Equation (3) is unweighted. No weighted cold
complex extension of (2) is assumed.

## 2. Elementary analytic smallness lemma

Suppose f is holomorphic on the unit disk, |f|<=B there,
and |f(x)|<=a on the real diameter, where 0<a<=B. Set

    L=max(1,log(B/a)).

In the upper half disk the harmonic measure of the semicircle is

    omega(z)=(2/pi) Im log[(1+z)/(1-z)].

It is zero on the diameter and one on the semicircle. Applying
the maximum principle to log|f|, with the usual limiting
interpretation at zeros and corner points, gives

    |f(z)|<=a exp(omega(z)log(B/a)).

The same argument applies below the diameter. On |z|<=rho<=1/2,

    |omega(z)|<= (4/pi) artanh(rho)
                  <=16rho/(3pi)<2rho.

The logarithmic form and branch are the standard ones in
[DLMF 4.37](https://dlmf.nist.gov/4.37); the harmonic-measure
and maximum-principle argument is given here explicitly.
Choose rho=1/(2L). Then |f|<=e a on that smaller circle, so

    |f'(0)|<=2e a L,
    |f''(0)|<=8e a L^2.                            (4)

If a=0, the identity theorem gives f=0. Boundary continuity is
not needed: apply the argument on slightly smaller half disks
and take the limit. This proves (4) without treating a real
bound as a bound on the original complex disk.

The logarithmic price cannot in general be omitted. The
functions f_n(z)=exp(-n)cos(nz) are bounded by one on the
unit disk and by exp(-n) on the real diameter, but
|f_n''(0)|/exp(-n)=n^2.

## 3. Actual cold scalar first and mixed bounds

For f_x in (3), put a=1.7e-15 a_x. Since a_x>=h L_h>exp(-80),

    log(B_c/a)<80+log(2e-11/1.7e-15)<90.

For g_chi use a=1.7e-14 h L_h sqrt(chi). Then

    log(B_c/a)<90+.5 log(1/chi).

Neither a exceeds B_c. Differentiating the scaled parameter
costs 2/E for the first derivative and 4/E^2 for the second.
Equation (4) therefore gives the following bounds after rounding:

    |D S[b][u](x)|<2e-6 U a_x,
    |D^2 S[b][u,w](x)|<2500 U V a_x,               (5)
    |D S[b][u](chi)-D S[b][u](0)|
        <2e-5 U h L_h chi^(1/4),
    |D^2 S[b][u,w](chi)-D^2 S[b][u,w](0)|
        <25000 U V h L_h chi^(1/4).               (6)

For the bilinear estimate use polarization of the symmetric
second derivative on u/U +/- w/V. Both directions have norm
at most two. This costs a factor two compared with the
unit-direction quadratic bound, not a factor one.

For the spatial power in (6), set t=log(1/chi). For j=1,2,

    exp(-t/4)(90+t/2)^j <=90^j, t>=0.

Its logarithmic derivative is negative. Thus
sqrt(chi)[90+.5log(1/chi)]^j <=90^j chi^(1/4).
The unrounded quadratic constant in (5) is below 1200;
its polarized bound is below 2400. The enlarged 2500
and ten-times-larger oscillation bound have strict slack.

The response is C2 INTO the weighted space with quarter-power
bolt modulus. On chi>=delta>0 this follows from uniform
unweighted parameter holomorphy. Near zero the actual, sharper
sqrt(chi) log^j(1/chi) estimate divided by chi^(1/4)
tends uniformly to zero. Splitting at delta proves norm
continuity of the first and second derivatives. Integrating
the derivatives along real segments then proves the Frechet
difference quotient statements in this space. It does not
prove C2 into the earlier square-root weighted space.

That distinction is necessary from these inputs alone. With
0<q<=1, set a_q=a0 sqrt(q), n_q=log(B/a_q), and

    f_q(z)=a_q cos(n_q z), f_0(z)=0.

This is a bounded holomorphic C0-valued family on smaller
parameter disks, with real bound a_q. Its second derivative
divided by sqrt(q) grows like log^2(1/q). It obeys the
quarter-power estimate. This is an analytic diagnostic,
not a claim that this manufactured family is a Dirac source.

## 4. HH pressure and the canonical cold component

Let S_* be the frozen homogeneous canonical scalar at the
same B,A and mu0. Separate it exactly:

    Q=S-S_*, R=P-S_*/2,
    R=F_b[Q], F_b[f](c)=[2W_b(c)]^-1 integral_0^c W_b' f,
    W_b=exp(2sigma0+2b)sinh^2 chi.                  (7)

The homogeneous term remains in the full action/source
allocation; it is not discarded or retuned. The canonical
deviation profiles below must be combined with that fixed
term and the other sectors in the eventual full source.
The null combination 2P-S is unchanged by this separation.

S148 proves ||F_b||<=1.2, ||D_u F_b||<=700 U and
||D^2_uw F_b||<=5500 UV on the fixed interior cover.
Condition 3 and monotonicity of a_x sharpen these on [0,c] to

    ||D_u F_b||<=700 U a_c,
    ||D^2_uw F_b||<=5500 UV a_c^2, a_c=a_x(c).     (8)

Indeed its kernel differences are at most 2U a_c, and the
direction derivatives at most 8U a_c. Insert those bounds
in S148's exact kernel derivatives. All these operator
derivatives annihilate spatial constants.

The fixed-seed estimate (24) of the real canonical response is

    |nabla^j(Theta_seed-Theta_product)|/(j!100^j)
        <=1e20 q^2 a_x, j<=3.

After normalization, the value coefficient is below 2.72e-7.
For the first spatial derivative multiply by 100 and convert
proper distance to chi with exp(sigma0)<8. On chi<=1
also use a_x<=2hL_h. Combining with (2) yields

    |Q|<3e-7 a_x,
    |S(chi)-S(0)|<5e-4 h L_h chi^(1/4).           (9)

These are conservative finite-background VALUE bounds. In
particular the product term in (7) has not been assigned a
zero value.

In normalized units define the canonical deviation profiles

    Rho_c=R-Q, N_c=2R-Q.

They equal (2G/(B^2 C)) times the corresponding physical
hat profiles H_rho,c and H_N,c. The latter null profile is
the full canonical null profile. Equations (5)-(9), including
every changing-Ward-kernel term, give:

| Profile weight | Value | First / U | Mixed / UV |
| --- | ---: | ---: | ---: |
| Maximum absolute density/null deviation / a_x | 1.1e-6 | 9e-6 | 8600 |
| Absolute null profile / [h L_h chi^(1/4)], chi<=1 | .0017 | 7e-5 | 86000 |

For the first row, the pressure first bound is
1.2(2e-6)+700(.00331)(3e-7)<3.2e-6. Its mixed
bound is
1.2(2500)+1400(.00331)(2e-6)
              +5500(.00331)^2(3e-7)<3001.
Combining with the trace gives the table.

For the bolt row apply the operators to S-S(0).
The value coefficient is 3.4(5e-4)=.0017.
The first coefficient is
3.4(2e-5)+1400 a_c(5e-4)<7e-5.
The mixed coefficient is
3.4(25000)+2800 a_c(2e-5)+11000 a_c^2(5e-4)<86000.
Here a_c<=2hL_h<2e-28, so the kernel variations do not
lose the cold gap. No varying source gradient is used.

## 5. Absolute cold constraint consequence and its size

S152's integration proof applies to any integrable bolt
modulus. With chi^(1/4), only its bolt lapse coefficient
changes:

    3L_h integral_0^1 chi^(-3/4)dchi
       =12L_h<972,

in place of 486 for the square-root weight. Its other
regional coefficients and Jacobian bounds are unchanged.
This is an explicit variant, not a claim that (6) meets
S152's original square-root input condition.

To use the preceding normalized profile table in S152,
divide each entry by 2 eta_G to obtain Z_i,j. The overall
eta_G in the constraint table then cancels once, leaving
one half times the displayed normalized bounds.
Assume the separate low radius/shape geometry hypotheses
of S152, with their combined direction norm bounding U,V.
For the COLD canonical deviation component alone,

    ||D^2 mass_primitive|| <2e5 U V,
    ||D^2 absolute_mass_density_map|| <.0087 U V,
    ||D^2 lapse_density_map|| <1e8 U V.             (10)

The first is the fixed-seed zero-datum primitive norm,
the second the map into the absolute-density space in BC
units, and the third the L1 lapse-density norm. We do not
differentiate an absolute value as a scalar function.
The exact budgets use

    m2=1.003(8600)+4.02(9e-6)+5.54(1.1e-6),
    c2=1.006(8600)+10.08(9e-6)+26.5(1.1e-6),
    b2=1.006(86000)+10.08(7e-5)+26.5(.0017),
    mass: 46 m2/2,
    infinity: .000002 m2/2,
    lapse: [13122 c2+972 b2]/2.

These are finite but coarse. For example the lapse Taylor
bound at a path length 1e-6 is 5e-5, not a demonstrated
global error allowance. At length 1e-10 it would be
5e-13, PROVIDED the path remains admissible; no such
correction length or full inverse is supplied here.
The mere smallness of eta_G is not a contraction proof.

## 6. Focused verification and remaining proof

Before execution, verify the harmonic-measure boundary
identities and derivative budgets; the parameter rescaling,
polarization and logarithmic weight bounds; every Ward
and regional-constraint coefficient; and the quarter-power
integration. Independently compare analytic scalar families
against Cauchy contour integrals and real finite differences.
Use 60 digits, 64/128 trapezoidal contour nodes, real steps
1e-4 and 5e-5, relative-to-(1+|target|) tolerance 1e-7,
refinement ratio .35 unless both errors are below 1e-40,
and contour identity tolerance 1e-20(1+|target|).

Test manufactured families at thermal analogues exp(-20),
exp(-50), exp(-79) and bolt separations 1,exp(-20),exp(-80).
Normalize derivative diagnostics by their known amplitude so
tiny cold values do not make the tolerance vacuous. Use frequency
log(B/a)/16 for the stated finite-difference steps; use the full
frequency in the exact logarithmic-loss counterexample.
Also compare an exactly solvable Ward pressure and its absolute
bolt integral for a square-root/logarithmic analytic trace family.
Check that dropping logarithmic losses or asserting the
old square-root response space fails analytically. These
fixtures test the mathematical transfer, not the actual
continuum determinant. The actual-source inputs remain
the retained S145 and real endpoint proofs.

The result closes the canonical cold RESPONSE estimate on
the explicitly stated family contract, using an adequate
weaker bolt modulus. It does not prove family admissibility
for the constructed endpoint, solver preservation, the other
source sectors, the outer recombined response or the full
higher-order inverse. No physical completion gate is met.

Remaining to-do list: establish an admissible global family
and outer/other source responses; sharpen or use these bounds
with a full inverse; solve exterior/state, charge, return and
signed source; complete opening, stability, transport,
preparation and separate creation.
