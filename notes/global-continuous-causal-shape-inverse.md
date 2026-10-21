# Global continuous causal shape bound and variational inverse

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

For the retained charged seed there is a finite Gamma_G such that at
EVERY real Laplace frequency gamma>=Gamma_G the actual curvature-free
shape form has a closed auxiliary realization on H=L2((0,infinity),dy).
Its associated operator satisfies

    A_gamma,F>=gamma^2/(40G), ||A_gamma,F^-1||<=40G/gamma^2.       (1)

The canonical slot is the actual HH causal response, with its full
angular tower and unchanged contacts. The scalar model introduced below
is a comparison used to prove closure. The conclusion covers compact
shape inputs anywhere in the exterior and their global conformal tails.
It extends the [outer continuous theorem](outer-continuous-causal-shape-bound.md)
through the horizon. The threshold is not numerically calibrated.

This supplies a uniformly bounded strong feedback at all sufficiently
high real gamma, with the same curvature graph domain as before.
Continuity of that feedback in gamma, its physical boundary/state
interpretation, transfer to an exact nonlinear geometry and stability
are still unproved. In particular this note does not invoke the
conditional continuous resonance conclusion.

## 1. Global scalar comparison at the actual thermal anchors

Use physical optical time and the ACTUAL even doubled canonical orbit,

    gbar=dY^2+a(Y)^2 dt_E^2, a=sqrt(CF)/r, dY=a dy.

On the half orbit 0<=Y<=Y_h, the outer seam is Y=0 and the bolt is
Y=Y_h. The double is compact with two smooth bolts and its retained
C2,1 seam. At the seam a>0 and a'=0; at a bolt, in proper distance rho,
a=kappa_h rho+O(rho^3). Derivatives here are with respect to Y.
Set L=-Delta_gbar and Qsc[T]=<T,Kflat(sqrt(L))T>, using exactly the
multiplier and common fourth contact of the
[uniform comparison](uniform-logarithmic-frequency-comparison.md).
Divide even forms by two and by the actual time-mode norm, as before.

For every epsilon>0 there are finite K_epsilon,L_epsilon such that

    |Qgeo[T]-Qsc[T]|<=C[epsilon J2+K_epsilon J1+L_epsilon J0]     (2)

for all compact smooth thermal Fourier sums on the open half orbit,
uniformly in their support and frequency. Here Jj are the global scalar
norms, not the outer auxiliary-cylinder norms. The extension of the
previous kernel argument to the whole double is as follows.

Cover the compact double by finitely many normal neighborhoods: smooth
Cartesian neighborhoods at the bolts and in the interior, and the
already controlled C2,1 neighborhoods at the seam. Choose one positive
normal radius and finite curvature bounds for this fixed geometry.
The Jacobi estimate |Delta r-1/r|<=k r and |j-1|<=k r^2, and the
Gaussian barriers exp(+-2kt+-kr^2/4)(4pi t)^-1 exp(-r^2/(4t)),
give precisely the local comparison in the preceding note.

The full heat kernel has the necessary uniform return bound on this
compact double. A finite coordinate partition gives the two-dimensional
Nash inequality ||u||_2^2<=K||u||_1||grad u||_2+K||u||_1^2.
Applying it to the conservative heat semigroup and differentiating its
squared L2 norm gives ||P_t||_(1->infinity)<=K(1+t^-1).
Within each smaller coordinate ball the diffusion coefficients and
drift are bounded; the exponential martingale bound for exiting a
fixed ball is K exp(-c/t). The strong Markov property, a time split
and symmetry combine this with the preceding kernel bound to bound
returns by K t^-1 exp(-c/t). This also justifies the killed-kernel
boundary correction in the barrier argument. For large times the
kernel is bounded, with its constant zero mode retained. All these
terms are integrable against the unchanged angular heat weight
beta(t)<=.8 C pi^-1 t^-3 exp(-648t).

Consequently the separated scalar/geodesic kernel difference is bounded
by KC/r^4 plus a bounded tail, uniformly on the double. Its common
fourth restoration differs only by
C F_eff/(360pi) integral Rbar |grad T|^2. Both global models annihilate
the constant scalar. The finite second contacts are fixed by the same
cutoff coordinate quadratics q_x with grad q_x(x)=0. In the smooth caps
all required jets are bounded; in the seam caps the previous bound
|L^2 q_x(z)-L^2 q_x(x)|<=K distance(x,z) applies. The resulting heat
remainder is O(t^(5/2)), and the geodesic remainder is O(r^5).
Their summed contacts are finite because sum d_m/M_m^3<infinity.
Long geodesic returns, including conjugate points, have bounded row
mass by Liouville invariance and the retained Bessel moments.
The anchored gain-two difference estimate and the integrated Bochner
identity now give (2), with arbitrary epsilon. This proof includes
the bolt charts; it does not bound the singular polar time coefficient.

Combine (2), at epsilon=.001, with the retained GLOBAL actual/geodesic
error .001C J2+3e24C J1+3e26C J0. At gamma_n=n kappa_h, s_n=gamma_n^2,
the actual causal form D therefore obeys

    |D(s_n)-Qsc(s_n)|<=C[J2(s_n)/500+A0 J1(s_n)+B0 J0],          (3)

where A0,B0 are finite nonnegative constants. No time rescaling is
made. Equality of causal and Euclidean responses at these nonzero
anchors is the retained thermal spectral identity.

## 2. The unbounded bolt multiplier has a radial graph bound

For all real s>0 define the scalar Friedrichs operator

    H_s=A+s B_tau, A=-a^-1 partial_Y(a partial_Y), B_tau=a^-2,
    H_s>=s b0, b0=(max a)^-2>0,
    J0=||v||^2, J1=<v,H_s v>, J2=||H_s v||^2                  (4)

in L2(a dY), with even seam condition. It agrees with the Fourier
restriction of L at every anchor. The positive forms defining H_s have
a common domain for all s>0, even though B_tau is unbounded. Put

    K_B=sup(4a'^2-2aa'')_+, K_V=sup(a'^2/4-aa''/2)_+.

Both constants are finite, including at the bolts. Direct integration
by parts on compact inputs gives the exact identity

    ||H_s v||^2=||A v||^2+2s integral a B_tau |v'|^2
                +integral a[s^2-s(4a'^2-2aa'')]B_tau^2 |v|^2.

Here Delta_rad B_tau=(4a'^2-2aa'')B_tau^2. Thus for s>=2K_B,

    ||B_tau v||<=sqrt(2) ||H_s v||/s<=2||H_s v||/s.             (5)

Only this compact-input graph bound is needed below. Resolvent
differentiation uses the form bound B_tau<=H_s/s, not an assumed
operator norm for B_tau. Under the half-density unitary, H_s is

    -partial_Y^2+(s+aa''/2-a'^2/4)/a^2.                        (6)

For gamma/kappa_h>1 its inverse-square endpoints are limit point.
The even condition at the seam remains fixed. We will use gamma>2kappa_h.

## 3. Factor the fourth order before differentiating the logarithm

Write Kflat(sqrt(x))=x^2 g(x), c_alpha=alpha(1-alpha). The identity
R(z)/z^2=integral_0^1 u[1-(1+zu)^-1]du gives

    g(x)=-C F_eff/(120pi)
       +sum_m d_m/(2pi M_m^2)
          integral_0^1 c_alpha^2 integral_0^1
                 u[1-(1+c_alpha u x/M_m^2)^-1] du d alpha.     (7)

Here M_m^2=m(m+q)/C, d_m=2(q+2m), m>=1 and q>=1.
The counting function satisfies

    N(T)=2N(N+q)+2N<=3CT,
    sum_m d_m/(M_m^2+z)^2<=3C/z, z>0.

The first N on the right denotes the largest admitted integer m;
if none is admitted the count is zero. Integrating the counting bound
proves the second inequality. Using integral c_alpha^2=1/30 and
integral u=1/2 in (7) therefore gives

    0<=x g'(x)<=C/(40pi).                                    (8)

For one resolvent R=(I+tH_s)^-1 put
Y=t R^(1/2) B_tau R^(1/2). Its form bound gives 0<=Y<=1/s.
The j-th derivative of 1-R has sign (-1)^(j+1) and equals
j! R^(1/2)Y^j R^(1/2) with that sign removed. Thus Y^j<=s^(1-j)Y,
and B_tau<=H_s/s followed by (8) gives, for G_s=g(H_s),

    0<=(-1)^(j+1) G_s^(j)<=j! C/(40pi s^j) I, j>=1.           (9)

No B_tau factor is commuted through a resolvent. One may first truncate
the positive mass integrals and sum. Uniform bounds for derivatives
through order four and their strong monotone limits justify the first
three form derivatives on compact vectors by integration in s.
Operator-norm convergence of the angular tails is not required or
asserted. The undifferentiated g is logarithmically unbounded.

For fixed compact v set h_s=H_s v and b=B_tau v. Product differentiation
of Qsc(s)[v]=<h_s,G_s h_s>, with h_s'=b, gives

    Qsc'''=<h_s,G_s'''h_s>+6 Re<b,G_s''h_s>+6<b,G_s'b>.

Equations (5) and (9) imply

    |Qsc'''(s)[v]|<=C K3 J2(s)/s^3, K3=27/(20pi).              (10)

This is the global derivative estimate. Its right side is a graph
norm, allowing supports arbitrarily close to the bolt.

## 4. Four actual anchors give a global continuous comparison

Use a=s_(n-1), b=s_n, c=s_(n+1), d=s_(n+2), s in [b,c]. The retained
positive spectral measure gives the exact bracket

    I_(b,c,d)D(s)<=D(s)<=I_(a,b,c)D(s).

All fixed quadratic subtraction polynomials are reproduced. For n>=20
both sums of absolute interpolation weights are <=3/2, and throughout
[a,d], |t-s|/s<=1/4. Equations (4)-(5) give

    J2(t)<=(3/2)^2 J2(s), J1(t)<=(3/2)J1(s).

Require the whole interval above 2K_B. The two-sided anchor error (3)
therefore costs at most C[27J2/4000+(9/4)A0 J1+(3/2)B0 J0].
The usual degree-two interpolation remainder costs at most
sup|Qsc'''| (d-a)(c-b)^2/24. Since t>=3s/4, (10) bounds that
supremum by (16/3)C K3 J2(s)/s^3. Elementary squared-node arithmetic
gives (d-a)(c-b)^2/(24s^3)<=2/n^3. The remaining cost is at most
72C J2/(5pi n^3)<C J2/1000. Hence

    |D(s)-Qsc(s)|<=C[31J2/4000+A1 J1+B1 J0],
    A1=9A0/4, B1=3B0/2.                                    (11)

This is uniform in support throughout the exterior. It assumes no
endpoint-uniform bound on the actual continuous spectral measure.

## 5. Global logarithmic compensation by operator Jensen

Take s>=K_V. The half-density operator (6) dominates the Dirichlet
free Laplacian on I=(-Y_h,Y_h). The variational resolvent inequality
compares that Dirichlet inverse with the compressed full-line inverse
on zero extensions. Integrating the log resolvent representation gives

    <u,log sqrt(H_s)u>>=<u_ext,log|D_Y|u_ext>
       >= integral -log|Y| |u|^2+D_H||u||^2,
    D_H=psi(1/4)+log 2>-4.                                  (12)

The last inequality is the retained logarithmic uncertainty bound,
obtained by differentiating equation (3), page 1 of
[Beckner](https://arxiv.org/pdf/math/0612486v2), with its Fourier
normalization converted as in the preceding audit. Use a positive
spectral regulator before taking the free logarithm limit.

Use the unchanged even partition chi_o=cos(theta), chi_i=sin(theta),
transitioning between |Y|=.01 and .02. Put W_chi=sum|chi_j'|^2 and
W_*=||W_chi||_infinity. The exact second-order IMS identity is
sum chi_j H_s chi_j=H_s+W_chi. Operator concavity of log gives

    sum chi_j log H_s chi_j<=log(H_s+W_chi)
        <=log H_s+W_*/(s b0) I.                             (13)

For unbounded H_s the first inequality follows directly by applying
the variational inverse/compression inequality to H_s+t and integrating
its resolvents. Thus it is not an unsupported substitution into the
bounded-operator version of
[Hansen--Pedersen, Theorem 2.1](https://arxiv.org/pdf/math/0204049).

Apply (13) to u=H_s v, after half-density conjugation. On chi_o u,
(12) and omega<=-log|Y|-.5log(x/C)+.5log 2, together with
-log(mu0 ell_A)>9, leave the same coefficient >33/20 after subtracting
omega+log(mu0 ell_o)+38/15+19/60. On chi_i u, omega<6 and
log sqrt(H_s)>=.5log(s b0). Choose a finite threshold such that

    .5log(s b0)-6-log(mu0 ell_o)-38/15-19/60>=2,
    s b0>=10W_*.

Then the loss from (13) for log sqrt is <=1/20, leaving 8/5 globally.
The retained scalar multiplier lower bound therefore gives

    Qsc(s)[v]-C/(30pi) integral(omega+19/60)|H_s v|^2
                         >=4C J2/(75pi).                    (14)

Localization acts on H_s v inside the logarithm. It never assumes
that a fourth-order expression commutes with a cutoff.

## 6. Closure of a comparison form and restoration of the actual response

Let a_sc,gamma be the full original classical and Weyl curvature-free
form with Qsc in the canonical slot. Use the retained exact global
reconstruction p=Q_gamma w/r, w=r v. The
[classical estimate](global-classical-curvature-free-form.md) holds
for real gamma. The [Weyl square completion](global-thermal-curvature-free-form.md)
is a polynomial identity in gamma^2 with residual <1e-24 N_gamma/G.
For compact v its global p has outer O(y^4) jets and horizon
rho^nu times a smooth function of rho^2, nu=gamma/kappa_h.
For nu>2 all fourth-order quadratic integrals and boundary fluxes
converge by the power counting in the preceding continuous audit.
This does not call a noninteger mode a smooth Euclidean harmonic.
Equations (14) and those retained estimates give

    a_sc,gamma[w]>=N_gamma/(33G)+4C J2[v]/(75pi)
                                      +C J2[p+v/3]/(40pi).   (15)

We first prove this model form closable in physical H. For compact
v supported away from seam and bolts, v belongs to D(H_s^3).
Equation (8) shows Kflat(sqrt(x))=O(x^2(1+log(1+x))), so
z=Kflat(sqrt(H_s))v is in L2(a dY). Its physical w-row output is
f_sc=(CF/r^3)z, and

    ||f_sc||_H^2<=sup(CF/r^4) ||z||_(L2(a dY))^2<infinity.

The local Euler outputs fw,fq belong to H as in
[the thermal operator proof](global-thermal-shape-inverse.md).
The outer argument is unchanged: fw,fq=O(y). At the noninteger horizon
four proper derivatives give at worst rho^(nu-4); the optical volume
contributes rho^2, so fw,fq=O(rho^(nu-2)), square integrable in dy for
nu>2. The bounded reconstruction adjoint, ||Q_gamma^*||<4/5, gives
the full local output fw+Q_gamma^*fq in H. Endpoint integrations have
zero flux by the same jets. Consequently a_sc is the form of a densely
defined symmetric positive core operator D=C_c^infinity(0,infinity)
in w. Its Friedrichs form closure exists. This is an operator-to-H
argument; positivity alone has not been used to assert closability.

Now let e_gamma=D(s)-Qsc(s). In addition to all preceding thresholds
require s b0>=max(1,2000 A1,sqrt(2000 B1)). Since
J1<=J2/(s b0), J0<=J2/(s b0)^2, (11) gives

    |e_gamma[v]|<=7C J2[v]/800,
    |e_gamma[v]|<=eta a_sc,gamma[w], eta=525pi/3200<3/5.        (16)

Hermitian polarization extends e continuously to the model form
completion. The relative bound is strictly less than one, so the
[KLMN theorem, Teschl 6.24](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf)
gives a closed actual form a_gamma=a_sc,gamma+e_gamma on precisely that
completion, with equivalent norms and the original compact core dense.
Its restriction to the core is the unchanged actual causal form.
We do not need a pointwise nonthermal Dirac-kernel mapping to H.

Since 4/(75pi)-7/800>1/(200pi), the actual form retains

    a_gamma>=N_gamma/(33G),
    a_gamma>=Gv/(40G)+C J2[v]/(200pi)+C J2[p+v/3]/(40pi).      (17)

The second line uses Gv<=1.001 N_gamma. These bounds pass to closure
by the closed weak derivative maps and bounded reconstruction, proving
(1). The response bounds of the thermal inverse hold unchanged:

    ||w||<=40G||f||/gamma^2, ||r p||<=32G||f||/gamma^2,
    a_gamma[w]<=40G||f||^2/gamma^2,
    J2[v]<=8000pi G||f||^2/(C gamma^2).

Choose Gamma_G large enough that the preceding requirements hold for
all s>=Gamma_G^2, including n=floor(gamma/kappa_h)>=20 and
(n-1)^2 kappa_h^2>=2K_B. This is a finite definition, not an evaluated
benchmark frequency or a claim that these frequencies are physically
admitted within a controlled semiclassical approximation.

## 7. Uniform strong feedback and the next decisive task

The [mixed row](thermal-curvature-shape-schur.md) is local in frequency
apart from the shape inverse just constructed. Its exact residual
contacts and reconstruction identities hold at real gamma. From (17),

    ||D_gamma||_(E_gamma->H)<=D_*,
    D_*=epsilon K2 sqrt(200pi)+epsilon K1 sqrt(33G)
                                      +K0 sqrt(40G)/Gamma_*,
    ||K_cross^dagger||_(H->E_gamma)<=C_R D_*/(2G),
    Gamma_*=max(Gamma_G,Gamma_R), C_R=1/(18lambda).

The [strong remainder bound](high-frequency-zero-fast-injectivity.md)
therefore holds at EVERY real gamma>=Gamma_*, with its same finite
formula for M_* and unchanged curvature graph domain. The minimal
operator is injective with range of codimension one once
C_a M_*/gamma^2<1. On the maximal domain its homogeneous kernel has
dimension one and its reflecting phase satisfies the same pointwise
O(gamma^-2) comparison in the
[boundary-phase theorem](coupled-boundary-phase-winding.md).

This establishes the uniform-bound portion of that theorem's additional
continuous-family hypothesis. Parameter continuity of the eliminated
strong feedback V_gamma remains to be proved. Neither separate
closures nor analyticity on fixed compact probes identifies or controls
the varying form domains. Do not infer fixed-phase resonance sequences,
analytic poles, or Lorentzian instability from this pointwise result.

Verification checks the radial identity, angular counting and resolvent
budgets, Jensen localization, interpolation and relative-form margins.
Independent smooth-bolt quadrature and noncommuting matrix derivatives
and logarithms test the mechanisms at two precisions. These fixtures
are not actual charged spectral data; the continuum proof and finite
thresholds remain ordinary mathematics, with no Lean run selected.
No action term, source condition, boundary phase or physical state is
changed. The seed remains off shell; S167 and S168.3-S168.4 stay open.

Remaining to-do list: continuity of the full eliminated causal response
and physical boundary/constraint/state admissibility; quantitative
thresholds and signed static compatibility; nonlinear existence and
all-sector stability; exterior/state, charge/return, opening, transport,
preparation, separate creation and full audit.
