# Continuous retarded growth rates have the proved outer normal form

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The actual centered seed's coupled retarded metric form has the same
outer annular limit at continuous positive imaginary frequencies as
the previously proved thermal Euclidean form. The proof uses the
positive third divided difference of the actual spectral response and
four neighboring nonzero Matsubara values. It does not fit an analytic
continuation or replace the fixed local contact polynomial.

More precisely, for fixed compact radial profiles and annulus scale
d->0, every continuous growth-rate sequence gamma_d with
d gamma_d->kappa>0 has that limit. The result also covers kappa=0
when gamma_d is at least twice the fixed Matsubara spacing. The
low-frequency interval below that threshold is not settled here.
This is convergence of core forms, including uniform convergence on
each fixed finite-dimensional profile space. It is not convergence
of full inverses, a physical boundary choice or a mode-existence theorem.

## 1. Fixed causal response and thermal samples

Use the same centered seed, full Dirac tower, neutral HH restriction,
prepared clock, zero-source domain and common subtraction as
[the thermal spectral derivation](thermal-shape-response-and-subtractions.md)
and [the coupled outer limit](outer-coupled-scalar-separation.md).
The [neutral current result](neutral-spherical-current-decoupling.md)
permits delta A=0 in the selected spherical sector. The background
remains off shell; no residual contact is dropped.

Use the normalized canonical time tau of the outer-limit note and
write h=2pi/beta_tau>0 for its FIXED Matsubara spacing. It is not a
new temperature or a value to be varied with d. The retarded convention
is exp(-i omega tau); omega=i gamma, gamma>0, represents a growth rate.
Let sigma=gamma^2. For a real compact pair (p,v), define D_(p,v)(sigma)
as the full fixed-density retarded quadratic response in the same
Fourier-amplitude convention as the Euclidean form divided by its
cosine time norm beta_tau/2.

Only the canonical v,v row is nonlocal. All conformal, mixed and
remaining local terms are polynomials of degree at most two in sigma,
as derived in the local frequency-polynomial note. The thermal
spectral measure in the nonlocal part is nonnegative and depends on
v, not p. The complete response therefore has the exact remainder

    D(sigma)-I_(a,b,c)D(sigma)
      =(sigma-a)(sigma-b)(sigma-c)
         integral dmu_v(t)/[(t+sigma)(t+a)(t+b)(t+c)],   (1)

for distinct positive a,b,c. I is quadratic interpolation of the
ACTUAL three values. The third subtracted integral converges for
each fixed compact probe by the preceding Hadamard/spectral argument.
Adding the local quadratic polynomial changes neither side of (1).
Thus all chart, seagull, measure and off-shell local contacts remain.

At nonzero Matsubara points b_n=h^2 n^2, n>=1, the KMS spectral
denominators give equality of this retarded response and the
Euclidean Hessian with the same contacts. The equal-energy static
covariance term integrates to zero at these frequencies. No zero-
frequency equilibrium Hessian is used as a causal sample.

The need for consistent retarded/Euclidean contact definitions is
established in [Romatschke--Son, section II, (13)-(17)](https://arxiv.org/pdf/0903.3946v1).
The retained thermal note proves the actual optical-channel spectral
identity used here. No quark-gluon-plasma sum rule is imported.

## 2. A four-sample enclosure on the positive axis

Choose n>=2 and sigma in [b_n,b_(n+1)]. Set

    L_n(sigma)=I_(b_n,b_(n+1),b_(n+2))D(sigma),
    U_n(sigma)=I_(b_(n-1),b_n,b_(n+1))D(sigma).         (2)

The product in (1) is nonnegative for the first interpolation and
nonpositive for the second. Hence, exactly,

    L_n(sigma)<=D(sigma)<=U_n(sigma).                  (3)

At a sample endpoint this is an equality. The inequalities hold as
quadratic-form inequalities on the common compact pair core: for
each pair the only spectral measure is its nonnegative mu_v. No
simultaneous spatial diagonalization or positivity of D itself is
assumed. In particular local gravitational signs may be indefinite.

There is a uniform numerical bound on cancellation in the weights.
For the lower interpolation only the b_(n+2) weight is negative;
for the upper only the b_(n-1) weight is negative. For sigma between
the middle sample points, their maximal magnitudes are bounded by

    |ell_lower,negative| <= (2n+1)^2/[4(4n+4)(2n+3)] < 1/8,
    |ell_upper,negative| <= (2n+1)^2/[16n(2n-1)] <=25/96.

The weights sum to one, so

    sum |ell_lower|<=5/4,
    sum |ell_upper|<=73/48, for n>=2.                  (4)

For example the slack in the last negative-weight bound is
(n-2)(26n+3)/[96n(2n-1)], which is nonnegative. The estimates are
independent of the thermal spacing. Negative weights are essential;
substituting only lower bounds at all three anchors would not give
a justified lower bound.

For ANY real target number Q and E=max_(j=n-1,...,n+2)|D(b_j)-Q|,
(3)-(4) imply the useful finite-sample estimate

    |D(sigma)-Q| <= (73/48) E.                         (5)

For a symmetric matrix of finitely many profiles, the same statement
holds in operator norm if E is the maximum sample matrix error norm.
Indeed each sample lies between Q-E I and Q+E I; the signed weights
and their absolute sums bound both interpolants as forms, and (3)
then bounds D. This is an error transfer, not a claim that the actual
sample errors have already been calibrated numerically.

## 3. Actual seed annular limit at continuous growth rates

Let p,v be fixed smooth profiles supported in a compact subinterval
of (0,infinity), and set p_d(Y)=d^-1/2 p(Y/d),
v_d(Y)=d^-1/2 v(Y/d). Let gamma_d>=2h and d gamma_d->kappa>=0.
Choose n_d=floor(gamma_d/h), taking n_d at the lower endpoint when
gamma_d is an exact sample. Then n_d>=2, and for each
j=n_d-1,...,n_d+2,

    |d h j-d gamma_d|<=2dh ->0.                        (6)

Every one of these four nonzero thermal sequences therefore has
the same scaled frequency kappa. The preceding uniform thermal
comparison and outer-limit theorem give

    d^4 D_(p_d,v_d)(b_j)/C -> Q_kappa[p,v]              (7)

for each of the four samples, with the actual local contacts intact.
The arbitrary-epsilon gain-two comparison is the input behind (7);
a fixed .001 fourth-order error would not prove this exact limit.

Apply (3) to each FIXED d and its fixed compact probe, then multiply
by d^4/C. Apply (5) with target Q_kappa[p,v]. The four errors in (7)
tend to zero, and the weight bound is independent of n_d. It follows
that

    lim_(d->0) d^4 D_(p_d,v_d)(gamma_d^2)/C
        =Q_kappa[p,v].                                (8)

In particular d gamma_d->kappa>0 automatically makes gamma_d>=2h
eventually, so (8) covers every such continuous positive growth-rate
sequence. It also covers fixed gamma>=2h with kappa=0. No assertion
is made here about a fixed gamma in (0,2h), a zero-frequency limit,
or arbitrary off-axis complex frequencies.

No spectral-measure bound uniform on escaping supports is assumed in
this transfer. The sign identity (1) is applied before the annular
limit, for each compact probe separately. Uniformity needed to pass
the limit is supplied by the four established thermal limits and
the elementary weight bound (4). This is why the earlier lack of a
numerical endpoint-uniform Lehmann constant does not obstruct (8).

Polarization gives convergence of all entries on a fixed finite set
of compact pairs. In finite dimension, entry convergence gives norm
convergence, so (8) is uniform on the unit sphere of any fixed finite-
dimensional profile space. The dimension and profiles must stay fixed.
This argument does not prove uniformity on an infinite-dimensional
unit ball or on families with growing radial oscillation.

## 4. The retained coupled retarded normal form

For clarity the limit in (8) is precisely the prior full coupled form.
In dmu=eta^-4 d eta, define

    D_k=eta^2(partial_eta^2-kappa^2)-2eta partial_eta,
    A_k=D_k-4, B_k=eta^2(partial_eta^2-kappa^2)-12,
    e=60pi/g+11/3, g=G lambda.

With F_k[v] the explicit direct-shape form of the outer-limit note,

    Q_k[p,v]=<p,(D_k+e)A_k p>_mu/(20pi)
             +<p,(D_k+e)B_k v>_mu/(30pi)+F_k[v].        (9)

All frequency signs here refer to omega=i gamma, not real oscillatory
Lorentzian frequency. The stable inverse A_k^-1 and the change
s=p+(1/3)A_k^-1 B_k v separate this NORMAL form as already proved:

    Q_k=<s,(D_k+e)A_k s>_mu/(20pi)+R_k[v],
    R_k[v]>=19||(-partial_eta^2+kappa^2)v||^2/(150pi)
                  +(1/(27g)+23/(4860pi))G_k[v].         (10)

The fast factor D_k+e remains. The positive closed realization and
inverse estimate of R_k for kappa>0 retain their NORMAL comparison
scope. The change of variables can have noncompact auxiliary tails;
it is applied after the core limit, in its already established normal
domain. No unproved interchange of an actual exterior inverse with
the limit is used.

Thus the coupled normal object is now identified with the actual
seed's retarded positive-imaginary-axis asymptotics, not only with
thermal Fourier samples. This closes that causal-continuation step.
It does not promote the seed to an exact solution. Pointwise or
finite-dimensional form convergence alone does not give an actual
infinite-dimensional inverse, a metric pole or a physical unstable
mode. Neither does the positivity of R_k settle the remaining factor.

## 5. Verification and next required transfer

The focused checker verifies (1)-(5), both negative weights, their
uniform constants, preservation of arbitrary quadratic contacts and
the positive-node requirement. Independent quadrature of a continuous
noncommuting matrix spectral measure tests both interpolation bounds
and their exact integral remainders. A separately scaled finite
spectral family tests the annular error-transfer mechanism as the
thermal index grows; it is not substituted for the exterior spectrum.
A negative-weight control detects invalid use of only lower anchor
bounds. The continuum result remains the ordinary proof above, with
its imported thermal and state/domain hypotheses. No Lean is selected.

The next step must control the actual two-row retarded metric problem
beyond its fixed-profile normal limit, with physical gravitational
boundary/initial data and variable-curvature mixing. Both off-shell
contacts must remain until an exact background or a justified limit
removes them. Use the neutral electric decoupling and Ward recovery;
the two-zero-trace boundary test still prevents a silent fast-branch
exclusion. No physical ultraviolet prescription or finite-time
stability result is supplied. S168.3-S168.4 remain open.

Remaining to-do list: actual charged metric boundary, quantitative
retarded transfer or justified physical exclusion; S167 signed
coefficient, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, full-sector stability, transport,
preparation, separate creation and full audit.
