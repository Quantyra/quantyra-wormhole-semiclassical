# Curved driven boundary coefficients and prescribed-source contacts

S074, 2026-10-13. The paired free Dirac domain now has explicit local
drive-dependent heat coefficients on matched curved collars, including
the curvature term absent on a product. Their Weyl transfer is integrable.
The resulting local relative transfer vanishes at a fixed round RN-AdS
conformal boundary after a common reference subtraction. Source-only
finite contacts do not change matter evolution when the source and
boundary geometry are prescribed. This does not compute the independent
nonlocal response, interacting state, dynamical controller or opening.

## Matched collars and the squared domain

Retain the actual paired domain of
[coupled-fermion-domain.md](coupled-fermion-domain.md), one massless
minimal complex Dirac field and matched spin/gauge frames. Both collars
have identical pulled-back metric and gauge jets. Fold them onto one
boundary, with eight spinor components. Unequal collar perturbations
are outside the following theorem. Work in Euclidean signature, with
inward unit normal n and

```
L_ab=-<nabla_a n,e_b>, L=gamma^ab L_ab,
Lhat_ab=L_ab-L gamma_ab/3.
```

This L convention agrees with the inward-normal mixed heat coefficients;
it is minus the second fundamental form defined as nabla n. Write
Gamma_n,Gamma_a for Hermitian Clifford matrices. In the folded bundle
there is a parallel anti-Hermitian matrix G with G^2=-1 and
{G,Gamma_i}=0. An explicit representation, in flavor/radial/angular
order, is

```
Gamma_n=tau_3 tensor sigma_1 tensor 1,
Gamma_0=-1 tensor sigma_3 tensor 1,
Gamma_1=-1 tensor sigma_2 tensor sigma_1,
Gamma_2=-1 tensor sigma_2 tensor sigma_2,
G=i tau_2 tensor sigma_1 tensor 1.
```

Set V=exp(theta G/2), chi=V Gamma_n V^dagger and Pi=(1+chi)/2.
The source theta=2 atan h may now vary over the boundary; extend it
normally constant to describe the boundary operator. Define
u_a=D_a theta and C_a=Gamma_n Gamma_a G. The C_a are Hermitian,
obey {C_a,C_b}=2gamma_ab and commute with Pi.

For the covariant Dirac operator D, Green's formula gives the adjoint
allowed projector Pi_adj=1-Gamma_n Pi Gamma_n. Thus the domain of
D^dagger D requires both Psi=Pi Psi and D Psi in Pi_adj. The latter
condition is Pi Gamma_n D Psi=0. Since
Pi Gamma_n Gamma_a Pi=0, differentiation of the first condition gives
mixed Dirichlet/Robin data

```
(1-Pi)Psi=0,  Pi(nabla_n+S)Psi=0,
S=Pi Gamma_n Gamma^a(D_a Pi)Pi
 =-L Pi/2 + u_a C^a Pi/2.                             (1)
```

For the adjoint squared problem the base projector and angle reverse;
the scalar -L/2 term is the same and the u.C term changes sign. Both
have strongly elliptic mixed principal boundary symbols. This derivation
includes the geometric derivative of the projector and does not simply
carry the product-collar Robin term unchanged to a curved surface.

Constant theta is a local left/right unitary intertwiner of D and its
adjoint, because G is parallel and anticommutes with D. Hence both
squared operators have identical constant-angle local heat coefficients
on matched collars. This extends the preceding static spherical argument
locally; it does not imply global spectral equivalence.

## Full local derivative contribution through a4

Use Tr(f exp(-s D^dagger D))~sum_j s^((j-4)/2)a_j(f).
Delta denotes the terms containing derivatives of theta, relative to
a constant angle on the same collar. The result is

```
Delta a_0=Delta a_1=Delta a_2=0,
Delta a_3(f)=1/[4(4pi)^(3/2)] integral_boundary sqrt(gamma) f u^2,
Delta a_4(f)=1/(4pi)^2 integral_boundary sqrt(gamma)
              [(4/5) f Lhat^ab u_a u_b+(4/15)(n f)u^2]. (2)
```

The boundary integral is over the common folded boundary, not a further
sum over two copies. The adjoint squared operator gives the same (2).
The result includes angularly varying sources and arbitrary symmetric
L_ab in three dimensions. It assumes a matched minimal Dirac connection;
additional potentials, unmatched electric/spin frames and boundary
degrees of freedom need their own calculation.

For an ordinary proof, evaluate traces in a frame with theta=0. The
geometric derivative is D_a chi=-L_ab Gamma^b+u_a G Gamma_n. In
particular its contracted square contains the matrix cross term

```
(D_a chi)(D^a chi)=(L_ab L^ab+u^2)1-2L_ab u^a C^b.
```

With rank eight, tr(S_1^2)=u^2 and
Delta tr[(D chi)^2 S]=-2L u^2-4L_ab u^a u^b, where
S_1=u.C Pi/2. Substitution into the mixed coefficient formulas leaves
(4/5)(L_ab-L gamma_ab/3)u^a u^b in a4 and (4/15)u^2 in its
normal-smearing term. Omitting that cross term would leave a spurious
nonzero unsmeared contribution on an umbilic surface. The curvature/
gauge terms linear in u vanish by the folded Clifford traces, and
tr S=-2L is independent of theta, so the traced S_:aa term adds none.

The imported formulas are the mixed heat theorem in
[Vassilevich (5.27)--(5.33)](https://arxiv.org/html/hep-th/0306138v3).
Its normal convention and localized-versus-integrated distinction
are essential here. Equation (2) is a direct specialization with the
derived paired domain, not a substitution of reflecting-wall data.

## Weyl transfer and finite curvature contacts

Let Gamma_E=-log|det D|=-1/2 log det(D^dagger D), with one common
zeta subtraction and zero modes treated separately if present. The
following is its local drive-dependent Weyl sector, not its phase or
conformally invariant state functional. Under a unitary identification
of the spinor norms, D changes infinitesimally by -{omega,D}/2.
The determinant variation averages a4(omega) of D^dagger D and
D D^dagger; their drive terms (2) agree. Thus it is +Delta a4(omega)
in this Euclidean effective-action convention. This sign is distinct
from the Lorentzian trace convention in the bulk anomaly note.

Put b_theta=4/[5(4pi)^2], d_theta=4/[15(4pi)^2]. For g=e^(2sigma)h,
both densities sqrt(gamma)Lhat^ab u_a u_b and
sqrt(gamma)(n omega)u^2 are Weyl invariant. Consequently

```
Gamma_E,theta[g]-Gamma_E,theta[h]
 = integral_boundary sqrt(gamma_h)
       [b_theta sigma Lhat_h^ab u_a u_b+d_theta(n_h sigma)u_h^2]. (3)
```

This local transgression obeys the two-step Weyl cocycle exactly.
It retains the normal term found on the product collar and supplies
the independent traceless-curvature term. An umbilic boundary has no
unsmeared derivative logarithm, but generally still has the second
term in (3).

The functional J_theta=integral sqrt(gamma)L u^2 changes by
-3 integral sqrt(gamma)(n omega)u^2 under a Weyl variation.
Thus the d_theta term is the variation of the finite local term
-d_theta J_theta/3. Its coefficient is scheme dependent. The first
term instead multiplies the Weyl-invariant boundary functional
I_theta=integral sqrt(gamma)Lhat^ab u_a u_b. In the displayed zeta
convention mu partial_mu Gamma_E=b_theta I_theta; a finite coupling
c_theta(mu) multiplying I_theta has beta function -b_theta when the
complete action is made scale independent. This statement concerns
this free boundary sector; it does not match other EFT coefficients.

The a3 power subtraction is unchanged by curvature: it renormalizes
a dimension-one coefficient of integral sqrt(gamma)u^2. The finite
coefficient of that source term, the finite J_theta choice and the
finite c_theta must remain explicit. No one-loop heat computation
determines the apparatus implementing the source.

## Fixed round RN-AdS conformal boundary

For the explicit static exterior f(R)=1-2m/R+Q2/R^2+R^2/ell^2,
write the Euclidean metric as f times its optical metric. Its inward
optical coordinate has d rho=-dR/f, and its optical sphere radius is
B=R/sqrt(f). For a spherical prescribed theta(tau),

```
B'/B=f'/2-f/R=-1/R+3m/R^2-2Q2/R^3=:w,
Lhat_opt,00=2w/3,  sigma=log(f)/2,
n_opt sigma=-f'/2.
```

Here the prime on B is d/d rho, and on f it is d/dR. Direct substitution
into (3), including the sphere integral, gives

```
W_boundary(R)=4pi integral d tau (theta')^2
 [b_theta (log f/2)(2/3)(R^2/f)w
                    -d_theta (R^2/f)f'/2].           (4)
```

At fixed ell,m,Q2, the second bracket has leading term -d_theta R;
the first tends to zero. Subtract the same local functional for pure
AdS at the same asymptotic clock and round sphere. The differences
of (2/3)(R^2/f)w and -(R^2/f)f'/2 start respectively as
2m ell^2/R^2 and -3m ell^2/R^2. Therefore the relative (4) tends to
zero, with at most O(log R/R^2). A finite-cut clock rescaling to match
induced lapse changes this limit by zero as well. The power-divergent
individual term must not be interpreted as finite injected energy.

This closes this local relative Weyl-transfer limit for the stated
fixed round RN-AdS reference. It is not a simultaneous large-flux
error estimate, a proof for unequal driven collars, a full holographic
subtraction or a computation of the finite nonlocal determinant.

## Which source ambiguities affect matter evolution

For a prescribed real-time boundary source theta(U) and fixed boundary
metric, any source-only local functional adds a multiple of the identity
to the matter Hamiltonian. It changes the evolution operator by a common
phase. Density matrices and normalized matter correlators are unchanged.
Compactly supported interior metric variations likewise see no boundary
contact term. This applies to finite source contacts, not to unmatched
vertices involving dynamical matter, metric fields or charge operators.

For example, specify a real finite coefficient zeta_s by

```
I_s=-(zeta_s/2) integral sqrt(-gamma)(D theta)^2,
T_s,ij=zeta_s[theta_i theta_j-gamma_ij(Dtheta)^2/2],
O_s=zeta_s Box_gamma theta,
D^i T_s,ij=O_s theta_j.                              (5)
```

With gamma=-dU^2+ell^2 dOmega^2, its stored energy is
E_s=2pi zeta_s ell^2 theta_dot^2 and
E_s_dot=-4pi ell^2 O_s theta_dot. A smooth switch whose derivative
vanishes at both endpoints has zero net change in this contact energy,
even if its final constant angle differs from its initial one. Thus
this particular finite ambiguity cannot change the endpoint matter
excitation for an otherwise identical prescribed drive. It still changes
intermediate stress/work assigned to the source and must be kept in
the Ward identity.

The physical protocol has not yet established that its requested
theta(t(U)) is an externally fixed trajectory. If it is feedback on a
dynamical clock, or if the boundary/controller is dynamical, these
functionals are no longer c-numbers on the enlarged system and the
argument above does not remove their dynamics. The initial state,
nonlocal driven response, supporting sources and actual trajectory
remain requirements of S075. No capture margin is inferred from (5).

## Checks and decision

The checker evaluates the Clifford identities, induced Robin matrices
and full derivative heat terms for arbitrary symmetric L_ab, all three
source gradients and both adjoint domains. It checks conformal weights,
the Weyl cocycle and RN asymptotic coefficients. Independent curved
surface checks differentiate the unit normal and finite-angle projector
of a coordinate graph, then reconstruct the Robin and heat terms.
Independent integrated
variations of a source action with time-dependent lapse and sphere
radius check its stress and source terms, together with its Ward identity.
The omission control detects the lost geometry/drive cross term.

Retain (2)--(4) in the common local boundary subtraction. Distinguish
source contact choices from dynamical vertex matching, and specify the
actual clock/controller before assigning their physical work. This
advances S074.4/S074.5 and S075 inputs without closing their full gates.

Remaining to-do list: complete unmatched/unequal boundary and interacting
response, physical source/clock matching, other S074 local terms and
observable remainder; S075 state/drive; S076 stress/geometry;
S077 stability; S078 transport; S079 preparation/creation; S080 audit.
