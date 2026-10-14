# Complete mixed axial pair response on the global product

S074, 2026-10-13. The
[physical mixed l=1/l=3 cubic vertex](axisymmetric-axial-action.md) now
has a complete vacuum pair response bound, including the light l=3
branch and every pair representation. The relative correction is small
in the retained k=2000 family. This completes this particular nonlocal
channel in its stated cubic representative. The full local action and
renormalization, interacting state, physical boundary drive, geometry
and sufficient stability remain open.

## Fixed state, masses and representative

Use the canonical fields and nondifferentiated cubic representative of
the axial-action note, in dimensionless global AdS2 coordinates. The
free domain is Friedrichs and the state is its stationary ground state.
Put x=1/C in [0,1], gamma=G/a^2, and define

```
d=2+4x, Delta=(1+sqrt(1+4d))/2 in [2,3],
nu_minus=1+12x-sqrt(1+24x+24x^2),
nu_plus =1+12x+sqrt(1+24x+24x^2),
h_s=(1+sqrt(1+4nu_s))/2, h_minus in [1,3], h_plus in [2,5],
mu_plus =1+6x+sqrt(1+12x+12x^2),
mu_minus=24x^2/[1+6x+sqrt(1+12x+12x^2)].               (1)
```

At finite C the light masses are strictly positive. Stable rationalized
forms are used in every candidate computation. Let g_ps be the m=0
coefficient of -Psi_p q1 q3_s from the preceding action. For a fixed
external quadrupole magnetic component, summing internal magnetic
indices multiplies g_s g_s^T by 7/3. Define ghat_s=sqrt(7/3) g_s. Then

```
Qall=sum_s ghat_s^T M^-1 ghat_s
    =gamma (53+830x+944x^2)/15,
M=diag(mu_minus,mu_plus).                              (2)
```

The two internal fields in this interaction are distinct. Their Wick
pair kernel is W_Delta W_hs, with no identical-field factor of two.
The equal-dipole bilinear in the previous calculation has its own
different normalization. Derivative/EOM transformations of the cubic
action still require local contacts and transformed observables. They
have not been set to zero by adopting this representative.

## Unequal-mass spectral decomposition and an independent normalization

Specializing
[Cacciatori, Epstein and Moschella, Theorem 7.2, equations 7.135-7.138](https://arxiv.org/html/2403.13142v2#S7)
to dimension two gives

```
W_A W_B=sum_(j>=0) rho_j(A,B) W_(A+B+2j),
rho_j=(2A+2B+4j-1)/(4pi)
       Gamma(j+1/2) Gamma(j+A) Gamma(j+B) Gamma(j+A+B-1/2)
      /[Gamma(j+1) Gamma(j+A+1/2) Gamma(j+B+1/2) Gamma(j+A+B)],
rho_j>0, lambda_sj=(Delta+h_s+2j)(Delta+h_s+2j-1).       (3)
```

Their delta is 1/2 and their mass parameters are A-1/2,B-1/2. All
conditions of that theorem hold here. Its complete expansion, rather
than a finite numerical fit, supplies convergence and the spectral
support. No novelty is claimed for (3).

There is an independent Fock normalization check. A lowest-weight
vector of two *distinct* scalar species at total excitation N has
coefficients

```
c_(p+1)/c_p=-sqrt[(N-p)(N-p-1+2B)/((p+1)(p+2A))],
sum_(p=0)^N c_p^2=1.
```

At N=2j its central local product amplitude is the sum of the products
of the two one-particle mode amplitudes, weighted by c_p. Dividing its
squared value by the central lowest-mode weight for A+B+2j gives rho_j.
There is no sqrt(2) normalization for identical particles. The checker
compares this recurrence calculation with (3) for four mass pairs and
j=0 through 8, and independently compares Euclidean products of the
two global mode sums with the spectral expansion at twelve points.
These tests verify normalization; the cited theorem supplies the
infinite expansion. They do not purport to span the entire two-species
Fock space using the local scalar product.

## Uniform Casimir gap with the light branch retained

The following deliberately conservative estimate suffices. From
sqrt(9+16x)>=3+2x and sqrt(1+12x+12x^2)<=1+6x,
Delta>=2+x and mu_plus<=2+12x. Rationalization gives

```
nu_minus=120x^2/[1+12x+sqrt(1+24x+24x^2)]
        >=60x^2/(1+12x)>=60x^2/13.
```

The derivative of nu_minus is
12[1-(1+2x)/sqrt(1+24x+24x^2)]>=0, so nu_minus<=6 and h_minus<=3.
Since nu_minus=h_minus(h_minus-1), h_minus>=1+20x^2/13.
The heavy branch has a larger h. Thus D=Delta+h_s>=3+x+20x^2/13.
For either branch,

```
lambda_s0-mu_plus
 >= (3+x+20x^2/13)(2+x+20x^2/13)-(2+12x)
  = 5/2+41/452+(113/13)(x-91/226)^2
             +(40/13)x^3+(400/169)x^4 >5/2,
lambda_sj-mu_plus >=5/2+10j+4j^2.                      (4)
```

The last line uses D>=3. This is a vacuum pair Casimir gap on the
fixed global product. It is not a physical horizon-frequency gap or
a bound on the number-conserving response in an occupied state.

## Complete response, derivative and explicit tail

For 0<=z<=mu_plus define

```
beta_s(z)=sum_(j>=0) rho_j(Delta,h_s)/(lambda_sj-z),
Sigma(z)=sum_s beta_s(z) ghat_s ghat_s^T.               (5)
```

Log convexity of Gamma gives
Gamma(t)/Gamma(t+1/2)<=1/sqrt(t-1/2), t>1/2.
Apply it to the four ratios in (3), keeping Gamma(1/2)=sqrt(pi)
exact in the j=0 term. Since 2<=Delta<=3 and 1<=h_s<=5,

```
rho_0<=15/[2sqrt(6pi)],
rho_j<=(4j+15)/(4pi j^2)<=19/(4pi j), j>=1.
```

Combining these bounds with (4) and sum_(j>=1) j^-3<=3/2 gives

```
beta_s(z)<=3/sqrt(6pi)+57/(32pi)
          <17/24+19/32=125/96<4/3,
beta_s'(z)<(2/5)(4/3)=8/15,
0<=M^-1/2 Sigma(z) M^-1/2<=epsilon_mix I,
epsilon_mix=(4/3)Qall
           =(4gamma/45)(53+830x+944x^2).               (6)
```

Here pi>3 and 1/sqrt(2)<17/24 were used. Matrix positivity and the
trace in (2) give the last two lines. They retain both internal branches
and every magnetic component needed for one external m. Different
external m blocks are equivalent, so no further factor of five appears.

A sum through j=N-1, N>=2, has the explicit uniform enclosure

```
0<=beta_s-beta_s,N<=19/[32pi(N-1)^2],
0<=beta_s'-beta_s,N'<=19/[80pi(N-1)^2].                 (7)
```

This follows by integrating the j^-3 majorant from N-1 to infinity.
Thus (5) is a complete nonlocal response in this representation. Its
convergence does not fix the local terms in a different derivative
representative or the finite counterterms in the full theory.

## Positive quadratic comparison and combination with the prior channel

Use auxiliary scalar baths indexed by s,j with masses lambda_sj and
couplings ghat_s sqrt(rho_sj). Their two-point kernels match (3).
The sequences sqrt(rho_sj) need not be square summable; use the closed
mass-form domain sum_sj lambda_sj |b_sj|^2<infinity. By weighted
Cauchy-Schwarz and (6), its off-diagonal term satisfies

```
|2 sum_sj (ghat_s dot u) sqrt(rho_sj) b_sj|
 <=sqrt(epsilon_mix)[u^T M u+sum_sj lambda_sj |b_sj|^2]. (8)
```

If epsilon_mix<1, the form is closed and strictly positive. The
infinite Schur complement is M-zI-Sigma(z), uniformly convergent on
the interval in (5). Its two decreasing eigenvalues cross zero once
each, with the poles obeying

```
(1-epsilon_mix)mu_p<=mu_tilde_p<=mu_p.                 (9)
```

There are exactly two poles below the first bath mass, counted with
multiplicity. The auxiliary norm at a pole is
u^T[I+Sigma'(z)]u; the derivative estimate controls its residue.
This is the same closed-form construction used for the prior
[dipole-pair comparison](polar-dipole-dressing.md), now with two
distinct pair masses and a generally rank-two polar self-energy.
It defines a correlated stationary Gaussian comparison by diagonalizing
the positive mass operator. Auxiliary fields are a two-point-kernel
representation, not additional physical species or extra stress energy.

The free covariance between the old dipole-dipole bilinear and this
dipole-octupole bilinear vanishes: there is no free l=1/l=3 contraction.
Their selected quadratic self-energies therefore add, without an
independence claim about the full interacting fields or higher cumulants.
The preceding bound epsilon_dip=3gamma(1+25x+28x^2)/(4pi) gives

```
epsilon_selected=epsilon_dip+epsilon_mix<177gamma,
C=sqrt(1+12pi k^6), gamma=6C/[k^10(C-1)], k=2000.
```

At k=2000, epsilon_mix<2.761e-32 and epsilon_selected<2.901e-32.
The corresponding light-pole absolute shift bounds are below
1.373e-52 and 1.443e-52, respectively, in units of mass squared times
a^2. Exact candidate substitutions are retained in the receipt. These bounds
describe the two selected nonlocal channels on the fixed product;
they cannot yet be entered as the total physical error in the opening
budget. Local contacts and all other channels remain required.

## Verification and open physical matching

[The checker](../checks/mixed_axial_response.py) verifies the polynomial
gap certificate and rational bounds, distinct-field normalization,
Euclidean product expansion, gamma-coefficient recurrence, spectral
majorants and uniform tail enclosures. It compares the full matrix
Schur poles/residues with twelve independently diagonalized finite mass
matrices. The largest pole discrepancy is below 4.67e-12 and residue
normalization discrepancy is below 1.56e-15. Those finite tests support the analytic infinite-form proof
and do not replace it. Inputs, sources, versions and terminal outcomes
are retained in the accompanying receipt and review record.

Remaining to-do list: match the local action, nonlinear even constraints,
other angular/parity and fermion sectors, transformed boundary/source
terms, quantum measure and finite counterterms in one regulator. Then
construct the physical state and drive, stress and evolving geometry,
sufficient stability, finite-energy transport and preparation/creation.
