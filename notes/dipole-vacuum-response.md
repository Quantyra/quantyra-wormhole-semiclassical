# Dipole vacuum propagation and the distinction from causal preparation

S072, 2026-10-13. The complete leading bosonic dipole source now has a
regular particular response in both polar quadrupole branches on the
global product. Its positive two-particle spectral decomposition has a
uniform separation from both master masses. Spatial overlaps cancel
apparently small global frequency denominators. This bounds a specified
time-smeared sphere deformation in the stationary free vacuum comparison.

This calculation supplies propagation of the previously derived source.
It does not select the interacting state or causal initial data. The
homogeneous quadrupole field, induced response, fermions, interactions,
physical mouth boundaries and evolving geometry remain required. In
particular the small number below is not a physical capture-error bound.

## State, normalization and spectral support

Use the fixed product and canonical triplet of
[the dipole source derivation](dipole-quadrupole-source.md), with
`C=b^2/a^2>=1`, `kappa^2=8pi G`, `kappa^2 E0^2=(C+1)/b^2` and E0>0.
On the universal global AdS2 cover choose the Friedrichs boundary domain
and H2 vacuum. In `ds2^2=a^2 sec^2(sigma)(-d tau^2+d sigma^2)`, put

```
Delta=(1+sqrt(9+16/C))/2,        2<=Delta<=3,
omega_n=n+Delta,
u_n^Delta(sigma)=N_n^Delta cos^Delta(sigma) C_n^Delta(sin(sigma)),
(N_n^Delta)^2=n!(n+Delta)Gamma(Delta)^2
                 /[pi 2^(1-2Delta) Gamma(n+2Delta)],
W_Delta(X,X')=sum_n u_n(sigma)u_n(sigma')
                 exp[-i omega_n(tau-tau'-i0)]/(2omega_n).
```

Real normalized sphere quadrupoles obey
`<Q_m(X)Q_m'(X')>=(16pi/15)delta_mm' W_Delta(X,X')^2`, with the
traceless Wick bilinear Q_m fixed in the preceding note. This ordered
vacuum kernel contains only vacuum-to-two-particle matrix elements.
The full operator also has number-conserving matrix elements, relevant
to occupied states, which are not covered by this spectral gap.

The discrete AdS product expansion is established prior art. We use
[Cacciatori, Epstein and Moschella, section 3.1 and Theorem 7.2](https://arxiv.org/html/2403.13142v2#S7),
specialized to spacetime dimension d=2, delta=(d-1)/2=1/2 and their two
mass parameters equal to our Delta-1/2. Every positivity condition in
their (7.135) holds. Their opposite metric sign gives the same positive
mass squared. The original coefficient proof and noninteger convergence
argument are in [Bros et al., sections 6 and 7](https://arxiv.org/html/1107.5161v1#S6).
Those sections were reviewed, including their stated appendix dependencies;
the appendices were not newly reviewed in full. No new identity or
exhaustive novelty review is claimed. The specialization is

```
W_Delta^2 = sum_(j>=0) rho_j W_Dj,
D_j=2Delta+2j,       lambda_j=D_j(D_j-1),
rho_j=(4Delta+4j-1)/(4pi)
       Gamma(j+1/2) Gamma(j+Delta)^2 Gamma(j+2Delta-1/2)
       /[Gamma(j+1) Gamma(j+Delta+1/2)^2 Gamma(j+2Delta)] >0.  (1)
```

The following Fock-space argument independently identifies the support
and positivity needed for the response bound. One-particle lowering and
raising operators have matrix elements sqrt[n(n+2Delta-1)] and
sqrt[(n+1)(n+2Delta)]. A lowest-weight vector at total excitation N is
`sum_(p=0)^N c_p |p,N-p>`, where

```
c_(p+1)/c_p=-sqrt[(N-p)(N-p-1+2Delta)/((p+1)(p+2Delta))],
sum_p c_p^2=1.                                             (2)
```

The lowering recurrence fixes one vector, with exchange parity (-1)^N.
The symmetric local quadratic field selects N=2j. Its lowest weight is
D_j and its Casimir is lambda_j. Raising produces mutually orthogonal
representations because their Casimirs differ. At total excitation M
there are floor(M/2)+1 such descendants, exactly the dimension of the
symmetric two-particle level. They therefore span every finite energy
level and a dense subspace. Projecting Q_m(f)|0> onto these orthogonal
spaces gives positive spectral weights and convergence in its smeared
Fock norm. Covariance and the normalizable boundary behavior identify
each projected two-point function with a nonnegative multiple of W_Dj.
Center matrix elements from (2) test the normalization in (1).

This proves support and positivity without an analytic continuation of
the field into a different state. The checker independently compares
the matrix elements to the gamma formula for three Delta values and
j=0,...,8, and checks all symmetric levels through M=18. Those finite
comparisons verify conventions; the argument above and the cited
convergent expansion justify the infinite sum.

## Both master branches and a bounded spectral transfer

Keep the full mass matrix from
[the sourced product calculation](quadrupole-product-response.md).
Writing `u=(chi,Acal)^T=utilde+v Q_m` gives

```
(Box2 I-K)utilde=R Q_m,
K=[[(6+2C)/b^2,6 kappa^2 E0/b^2],[2E0,6/b^2]],
v=(9G/(4b^2),-3G E0/4)^T,
R=(3G(C-1)/b^4,9G E0/b^2)^T,
mu_plus=(C+6+sqrt(C^2+12C+12))/C,
mu_minus=24/[C(C+6+sqrt(C^2+12C+12))].                       (3)
```

Here mu are the eigenvalues of a^2 K. Let x=1/C in [0,1]. Then
`mu_plus<=2+12x`, since
`(1+6x)^2-(1+12x+12x^2)=24x^2>=0`. Also
`lambda_0=9+16x+sqrt(9+16x)>=12+16x`. Consequently

```
lambda_j-mu_plus >= 10+4/C >0.                              (4)
```

This is a Casimir separation in the vacuum pair representation, not
a frequency gap of a physical horizon or a claim about thermal states.
No small-mass elimination of the light branch was used.

On each spectral component define the regular AdS-covariant particular
response by `[v+(lambda_j/a^2 I-K)^(-1)R] Q_m,j`. Equivalently set
`chi_part=(G/b^2)F_chi(lambda_j) Q_m,j` and
`Acal_part=(G E0)F_A(lambda_j) Q_m,j`, summed over j, with

```
D(lambda)=lambda^2-(2+12x)lambda+24x^2,
F_chi=9/4+[3lambda(1-x)+36x+72x^2]/D(lambda),
F_A=-3/4+[9lambda x-12x-60x^2]/D(lambda),
9/4 <= F_chi <= 51/20,       -3/4 <= F_A <= 0.               (5)
```

The local v Q contribution is included. To prove the inequalities,
write lambda=12+16x+e, e>=0. The denominator is positive by (4).
The numerator of the correction to F_chi is nonnegative. Its upper
bound is equivalent to nonnegativity of

```
lambda^2-(12+2x)lambda-120x-216x^2
 =e^2+(12+30x)e+48x+8x^2.
```

The numerator of F_A+3/4 is `3x(32+28x+3e)>=0`. The inequality
F_A<=0 is equivalent to

```
lambda^2-(2+24x)lambda+16x+104x^2
 =e^2+(22+8x)e+120+56x+24x(1-x)>=0.
```

For indices A,B in {chi,Acal}, the complete particular covariance is

```
<u_A,m(X)u_B,m'(X')>_part
 =(16pi/15)delta_mm' S_A S_B
   sum_j rho_j F_A(lambda_j) F_B(lambda_j) W_Dj(X,X'),
S_chi=G/b^2,        S_Acal=G E0.                             (6)
```

For every real common test smear with finite Q variance, positivity
of each W_Dj(f,f) and (5) give

```
Var chi_part[f] <= [51G/(20b^2)]^2 Var Q[f],
Var Acal_part[f] <= [3G E0/4]^2 Var Q[f].                    (7)
```

Their same-smear mixed covariance is nonpositive for E0>0. This sign
is not a pointwise sign assertion about the complex ordered kernel.
Equation (6) defines the vacuum-to-pair part; it neither specifies the
full interacting vacuum nor removes intrinsic homogeneous master noise.

## Frequency coincidences, overlaps and causal initial data

The Casimir argument has a direct spatial counterpart. Let w_k be a
normalized mode with weight D_j, frequency Omega=D_j+k and spatial
operator `H_lambda=-partial_sigma^2+lambda_j sec^2(sigma)`. Let u_n
be a quadrupole branch mode with weight h, mass mu=h(h-1) and
frequency omega=h+n. Green's identity gives

```
(Omega^2-omega^2) J_nk=(lambda_j-mu) I_nk,
J_nk=integral u_n w_k d sigma,
I_nk=integral u_n sec^2(sigma) w_k d sigma.                  (8)
```

The boundary Wronskian vanishes: h>=1 and D_j>=4 give endpoint
powers at least h+D_j-1>0. By (4), lambda_j-mu never vanishes.
Thus the forcing overlap I_nk vanishes at a frequency coincidence,
and its ratio to detuning has the finite continuation
`J_nk/(lambda_j-mu)`. This holds for every descendant and both
master branches, including parity-forbidden overlaps that are zero.
It explains why treating source projections as unrelated constants
would create a spurious small denominator.

For example, the dipole ground pair has profile `(u_0^Delta)^2`,
frequency 2Delta. Its overlap with heavy master n=2 is

```
I=A(h+2-2Delta),
A=N_2^h (N_0^Delta)^2 h sqrt(pi)
    Gamma((h+2Delta-1)/2)/[Gamma((h+2Delta)/2)(h+2Delta)],
I/[(2Delta)^2-(h+2)^2] -> -A/(2Delta+h+2).                  (9)
```

At C=1, Delta=3 and h=4, both overlap and detuning are exactly zero.
That source does not force the resonant mode. The finite continuation
of (9) selects a possible homogeneous component of a particular
solution; it is not division by zero or a unique physical response.

At a finite regulator, write q_part(t)=A_nk exp(-i Omega t).
The retarded solution with zero initial data for the same forcing is

```
q_ret(t)=A_nk [exp(-i Omega t)-cos(omega t)
                             +i(Omega/omega)sin(omega t)]. (10)
```

It solves `q_ret''+omega^2 q_ret=(omega^2-Omega^2)q_part`,
has zero position and velocity initially, and vanishes identically
when Omega=omega. In the symmetrized two-master basis the corresponding
identity is

```
u_ret(t)=u_part(t)-cos(sqrt(H)t)u_part(0)
                -H^(-1/2)sin(sqrt(H)t)dot u_part(0).        (11)
```

Equations (10)-(11) hold for regulated data or data in the required
operator domains. They do not bound the continuum initial derivative
variance: the factor Omega/omega can probe ultraviolet moments absent
from (7). A switching function multiplying the original derivative
source also introduces its commutator with Box2. It must be retained.
The physical initial state might instead contain correlated dressing;
neither zero data nor the stationary particular solution is selected
by this calculation. The stationary bound (7) cannot simply be assigned
to a switched junction.

## Explicit central time-smeared sphere deformation

At sigma=0 only even dipole modes occur. Their central Wightman weight is

```
a_n(Delta)=(u_(2n)^Delta(0))^2/[2(2n+Delta)]
 =Gamma(n+1/2)Gamma(n+Delta)
   /[2pi Gamma(n+1)Gamma(n+Delta+1/2)].                     (12)
```

Each gamma ratio is a positive beta integral. Thus a_n decreases with
n and Delta, and `a_n(Delta)<=a_0(2)=2/(3pi)`. Use a normalized real
Gaussian time window of standard deviation s=1/4 in tau; its Fourier
transform is exp(-s^2 omega^2/2). The window covers all times. Then

```
Var Q_m = (16pi/15) sum_(p,q>=0) a_p a_q
                         exp[-s^2(2p+2q+2Delta)^2]
 <=64/(135pi) sum_(m>=0)(m+1)exp[-(m+2)^2/4]
 <1/10.                                                   (13)
```

This is a complete infinite-sum bound, not a cutoff extrapolation.
For m>=5 the successive ratio is at most (7/6)exp(-15/4).
Bound the first five terms and the geometric tail with
`exp(-z)<=1/sum_(r=0)^32 z^r/r!`. The rational lower bound pi>157/50,
certified by Machin's formula and alternating arctangent series, then
gives 0.0969274932<1/10. The checker evaluates these rational
inequalities exactly; truncated double-mode and KL sums separately
check normalization in three geometries.

The fractional sphere metric coefficient and local area-element change
at linear order are `h_sphere=2 sum_m chi_m Y_2m`. Rotational invariance
and `sum_m Y_2m(n)^2=5/(4pi)` make its variance identical at every angle
and equal to `5 Var(chi_m)/pi`. Equations (7) and (13) therefore imply

```
RMS(h_sphere,part) <=(51/20)sqrt[1/(2pi)] G/b^2
                   <1.018 G/b^2.                          (14)
```

For the retained family `q=k^4`, `L4/l_p=k^5`,
`C=sqrt(1+12pi k^6)`, `b^2/l_p^2=k^10(C-1)/6`, with l_p^2=G,
the k=2000 bound is less than 1.215e-43. This measures only the
stationary vacuum-to-pair particular sphere response with the specified
Gaussian window. It is neither a total metric fluctuation nor an
all-time trajectory probability. The integrated linear quadrupole
area is zero; the local observable (14) is distinct from the monopole
JT area bound. No comparison to a physical capture margin follows
without the joint state, initial data and boundary matching.

## Verification and next physical step

[The checker](../checks/dipole_vacuum_response.py) combines independent
two-particle representation matrices, nine convergent Euclidean KL
comparisons, exact rational response certificates, five direct
near-resonant spatial quadratures and 24 further Green-identity
comparisons covering both branches. It verifies (10) and its initial
data symbolically. The initial development receipt predates the added
all-branch overlap checks; both snapshots and outcomes are retained.
The [acceptance receipt](../evidence/dipole-vacuum-response.json) records
source hashes and numerical scope; the
[source review](../evidence/dipole-vacuum-source-review.json) records
the literature, conventions, review limits and raw evidence.

The next calculation must control the difference between this reference
and the constrained interacting fermion/dipole state, including H3/H4,
induced currents, mixed correlations, local contacts and massive modes.
The causal identity identifies the additional initial-data problem;
it does not resolve it. Physical horizon matching and the evolving
geometry remain required before stress, capture or stability conclusions.

Remaining to-do list: control interacting and initial-state corrections;
match physical boundaries and the driven state; obtain complete
renormalized stress and evolving geometry; establish preparation,
sufficient stability and finite-energy transport.
