# Stationary exterior: geometric mass and the full boundary charge

S139/S074, 2026-10-15. Ordinary charge audit; no new numerical result.
The accepted fixed-classical-RN outer theorem makes the consistently
vacuum-subtracted geometric mass finite. It does not identify that mass
with the Hamiltonian of the complete renormalized action. This note
derives the Einstein contribution, gives an exact counterexample to that
inference in another well-posed theory, and isolates the unevaluated
Dirac charge response in the selected theory. It adds no finite coupling
to the selected zero-extra-coefficient model.

## The Einstein contribution under the retained clock normalization

Use one exterior, the physical area radius r and Lorentzian metric

```
ds^2=-exp(2alpha) F dU^2+dr^2/F+r^2 dOmega_2^2,
gamma_(0)=-dU^2+L_eff^2 dOmega_2^2, alpha(infinity)=0.
```

Keep the magnetic bundle and boundary angle theta=0 fixed. G is the
specified renormalized Newton coefficient. Assume the outer vacuum has
already been matched to L_eff and

```
F=r^2/L_eff^2+1-2G M_infinity/r+O((1+log r)/r^2),
alpha=O((1+log r)/r^4).
```

These are conditional metric asymptotics, not existence of the full
corrected exterior. For the Einstein action, its Gibbons--Hawking term
and the minimal AdS4 intrinsic counterterms, the charge of partial_U at
the spherical cut r=R is

```
E_EH(R)=exp(alpha(R)) sqrt(F(R))/G
        *[R^2/L_eff+L_eff/2-R sqrt(F(R))].             (1)
```

Its sign is fixed so a positive Schwarzschild--AdS mass has positive
energy and pure AdS has zero limiting energy. In fact

```
R sqrt(F)=R^2/L_eff+L_eff/2-G M_infinity L_eff/R
           +O((1+log R)/R^2),
E_EH(R)=M_infinity+O((1+log R)/R).
```

Thus this particular Einstein contribution tends to M_infinity, or
2M_infinity for the identical pair. The local Weyl action, its complete
boundary transgression and the Weyl-invariant determinant must still
be varied in the same boundary prescription. An action's finite value
is not its lapse derivative. The previously established outer source
cancellation is in [the radial-source note](stationary-weyl-radial-sources.md)
and [the massive boundary-layer note](stationary-massive-boundary-layer.md).
Its fixed-reference scope is unchanged.

## Exact counterexample with all boundary data specified

This subsection concerns neutral Schwarzschild--AdS in a different
classical theory. It is not an RN solution with an unchanged Maxwell
source and is not a correction assigned to the retained Dirac model.
Take

```
I_bulk=(16pi G)^(-1) integral sqrt(-g) f(R4),
f(R4)=R4-2Lambda+beta R4^2, Lambda=-3/L^2.
```

The metric field equation and its constant-curvature reduction are

```
f_R R_mu_nu-(f/2)g_mu_nu
 +(g_mu_nu Box-nabla_mu nabla_nu)f_R=0,
R0 f_R(R0)-2f(R0)=0,
-R0+4Lambda=0.                                      (2)
```

Consequently every Einstein metric with R_mu_nu=Lambda g_mu_nu,
including every member of

```
F_M=1+r^2/L^2-2G M/r, alpha=0,
R0=4Lambda=-12/L^2,
phi0=f_R(R0)=1+2beta R0=1-24beta/L^2,                (3)
```

solves (2). Choose beta nonzero and small with phi0>0. The scalar
curvature is constant, although the Schwarzschild--AdS Riemann tensor
is not constant when M is nonzero.

The higher-derivative variational problem has an extra boundary datum.
For beta nonzero introduce phi=f_R and the scalar-tensor form

```
I=(16pi G)^(-1) integral sqrt(-g) [phi R4-V(phi)]
 +(8pi G)^(-1) integral_boundary sqrt(|h|) phi K,
V(phi)=2Lambda+(phi-1)^2/(4beta).                     (4)
```

The usual orientation/sign of the timelike Gibbons--Hawking term is
understood. Its variation cancels the metric normal-derivative term,
but also contains K delta phi. Fixing the induced metric alone does
not remove that term: fix phi as well, equivalently delta R4=0.
[Guarnizo--Castaneda--Tejeiro, (3.21)--(3.25)](https://arxiv.org/html/1002.0617v2)
explicitly retain this condition. All variations of M in (3) preserve
the same phi0 and the same conformal boundary metric. At a finite cut
normalizable metric variations need not vanish; the comparison is
between states with the same asymptotic Dirichlet sources.

Here the charge normalization can be derived without interpreting an
on-shell action value. The derivative of the bulk Lagrangian with
respect to R_abcd is phi times its Einstein counterpart. In the
covariant symplectic potential and Noether two-form, all additional
terms contain a derivative of phi. On the family (3), these terms
vanish, delta phi=0 and phi0 is fixed. Therefore

```
Theta_f=phi0 Theta_EH,
Q_f[partial_U]=phi0 Q_EH[partial_U],
delta H_f=integral_(S_infinity) [delta Q_f-i_(partial_U)Theta_f]
         =phi0 delta M.                             (5)
```

This uses the Hamiltonian variation of
[Iyer--Wald, section 6](https://arxiv.org/html/gr-qc/9403028).
The Dirichlet completion and its counterterms can likewise be chosen
as phi0 times the Einstein ones on this constant-scalar sector. Indeed
V(phi0)=2Lambda phi0, so the restricted scalar-tensor action is precisely
phi0 times the Einstein action. Fixing the AdS reference to zero gives

```
H_f=phi0 M,  Delta H_f=phi0 Delta M.                  (6)
```

Independently, [Deser--Tekin, (36)](https://arxiv.org/html/hep-th/0212292v2)
give E=[1+2Lambda kappa(4c1+c2)] r0/(2G), where the
action is written R4/kappa+c1 R4^2+c2 Ric4^2. Comparing the ratio
of the R4^2 and R4 coefficients gives kappa c1=beta. Setting
c2=0 thus reproduces the factor 1+8Lambda beta in (6), independently
of the notation for the Einstein prefactor. If the added density is
c R4^2, beta=16pi G c and phi0=1-384pi G c/L^2.

The geometric mass coefficient in (3) remains M, while the energy
difference at fixed sources is phi0 Delta M. A boundary-source-only
energy offset cannot change this coefficient. One could define
G_phys=G/phi0 and then write F_M=1+r^2/L^2-2G_phys H_f/r;
that is a different normalization statement, not a proof that H_f=M
with the original fixed G. No such redefinition is silently made here.

This counterexample disproves the general implication from a finite
geometric mass to the Einstein-normalized full Hamiltonian. It neither
proves a nonzero correction in the selected c_i=0 model nor identifies
any new unknown coupling that must be added to that model.

## Absolute source offsets and state-dependent response are different

For one fixed, fully specified renormalized action and phase space,
charge prescriptions generating the same asymptotic symmetry can
differ by a function of the nondynamical boundary fields. This is the
scope of [Hollands--Ishibashi--Marolf, sections II--III](https://arxiv.org/html/hep-th/0503105v3).
It does not compare the Einstein part of one action with the full
charge of another action, nor prove that the Einstein part alone
generates the symmetry for the complete effective action.

If a finite intrinsic functional C[gamma_(0),A_(0),theta_(0)] is allowed
by the chosen boundary scheme, adding it shifts the absolute energy
by a source-dependent E_C. At the same boundary metric, magnetic
bundle, angle and laboratory clock it cancels from Delta H. Choosing
the reference energy or fixing all finite source contacts fixes that
remaining convention. Extra finite bulk coefficients being set to zero
does not by itself specify every finite boundary-source convention.
Conversely, an arbitrary finite extrinsic-curvature functional is not
automatically an admissible Dirichlet ambiguity: its normal-jet
variations can change the boundary problem unless a compatible bulk
completion or additional boundary data are supplied.

The factor in (6) changes with the state variable M even while all
sources and couplings stay fixed. It is therefore not E_C. In the
selected model, the corresponding issue is whether the already fixed
Dirac determinant and its required subtraction contribute a finite
state-dependent charge. This is a calculation, not permission to tune
another coefficient. Recentring at the same-q quantum reference removes
source-only offsets, but does not establish the coefficient of Delta M.

## What the full Dirac domain actually proves at theta=0

Use the massless vectorlike Dirac field and the identical-pair domain
of [the coupled-domain note](coupled-fermion-domain.md), retaining every
angular mode. In a normal Euclidean frame

```
chi=tau_3 gamma_n, Pi=(1+chi)/2,
chi gamma_n=gamma_n chi, chi gamma_U=-gamma_U chi.
```

The canonical radial Hamiltonian Green matrix is proportional to
gamma_U gamma_n. Its anticommutator with chi vanishes, so

```
Pi (gamma_U gamma_n) Pi=0.                            (7)
```

Thus the first-order classical radial surface bilinear vanishes for
arbitrary allowed angular traces at theta=0. This proof uses neither
the angular eigenvalue nor an LLL approximation. It removes a putative
classical matter surface Hamiltonian of that form. The explicit
calculation in [Bintanja--Freivogel--Rolph, section 3.1](https://arxiv.org/html/2308.00871v2)
instead restricts (3.2) to the LLL in footnote 8. Its on-shell classical
action argument cannot set a quantum determinant's metric derivative
to zero. A state expectation, common subtraction and the boundary
charge variation still have to be matched.

There is also an exact useful lapse identity before renormalization.
On the Euclidean metric with a=exp(alpha), first work at a finite cut
and fixed Euclidean period beta_T. The original spinor L2 measure is
a r^2 dr dU_E dOmega. Multiplication by a^(1/2) r maps it to fixed
coordinate L2 measure. In this representation

```
Dtilde=gamma_r [sqrt(F) partial_r+F'/(4sqrt(F))]
       +gamma_U/(a sqrt(F)) partial_(U_E)+D_S,E/r,
D_time=gamma_U/(a sqrt(F)) partial_(U_E),
delta_alpha Dtilde=-j(r) D_time,  delta alpha=j(r).    (8)
```

D_S,E includes the appropriate angular Clifford factor and magnetic
connection. The scalar rescaling cancels alpha'/2+1/r from the radial
connection; it does not discard the F' term. The physical projector is
unchanged in this orthonormal frame. The first-order complementary
adjoint domains remain the ones in
[the stationary phase note](stationary-paired-dirac-phase.md).
The induced square must be varied through its composition, rather than
holding a separately guessed Robin condition fixed. Equation (8) has
no j' vertex. It is a canonical lapse insertion, not a statement that
the renormalized lapse response vanishes.

## Precise remaining same-source loop-contact calculation

One can formulate the missing calculation without an unknown finite
coefficient. Keep beta_T, the inner state prescription and the magnetic
bundle fixed. Choose a smooth j_R equal to zero below R and one above
2R, with the outer regulator cut beyond 2R. Choose a smooth mass-direction
variation equal to delta F_m=-2G/r in that outer region and zero near
the bolt. This is an off-shell variation; it is not a change of the HH
period along a family of smooth classical bolts. At alpha=0, (8) gives

```
V_R=-j_R gamma_U F^(-1/2) partial_(U_E),
V_m=gamma_r [delta F_m/(2sqrt(F)) partial_r
             +delta F_m'/(4sqrt(F))
             -F' delta F_m/(8F^(3/2))]
       -gamma_U delta F_m/(2F^(3/2)) partial_(U_E),
W_mR=+j_R gamma_U delta F_m/(2F^(3/2)) partial_(U_E).
                                                               (9)
```

At a fixed finite-dimensional matrix regulator in the canonical L2
space, away from zero modes, write G_D=Dtilde^(-1) and
Gamma_D=-log|det Dtilde| for the regulated paired field; Dtilde and G_D
here denote the regulated matrix and its inverse. Differentiating this
definition twice gives the exact real modulus kernel

```
partial_m partial_(epsilon_R) Gamma_D
 =Re Tr[G_D V_m G_D V_R-G_D W_mR].                   (10)
```

The second term is an essential mixed vertex, not an optional contact.
The positive-square implementation must equivalently retain the paired
adjoint/modulus prescription. Zero modes, if present, require the same
separate projection/state prescription on both sides; a formal inverse
does not dispose of them. The accepted identical-pair phase argument
does not compute (10).

In the continuum (10) specifies a regulated calculation, not an
unsubtracted trace identity for divergent operators. A metric-dependent
covariant heat regulator has additional regulator-variation terms before
removal; these must be included when matching to the retained scheme.
The finite matrix identity alone does not justify omitting them.
Sum the Matsubara
and full angular radial-Green products in the common four-dimensional
regulator, include all local bulk and boundary counterterm derivatives,
then take the renormalized outer limit with the prescribed state and
clock. A different determinant convention must carry its local contacts
through that same operation. The lapse variation contains both the
integrated stress response in the transition/outer region and a possible
surface charge. They must be separated by the full differentiable
variational principle; (10) alone is not already a boundary energy.

The finite-cut local Weyl part of this operation is explicitly available
in [the boundary Weyl note](stationary-exterior-boundary-weyl.md).
The missing part is its combination with the Weyl-invariant functional
and the renormalized charge prescription. An endpoint a4 coefficient
fixes local ultraviolet subtraction but does not determine the finite
Green-product trace (10). Nor does convergence of the radial mass
integral justify interchanging the outer limit with the ultraviolet
angular/frequency sum. This note proves neither that the residual
state-dependent contact is nonzero nor that it vanishes.

The decisive matching statement would prove, in the selected scheme,
that the full same-source Hamiltonian variation equals delta M_infinity,
or compute its finite correction from the full action. An absolute
boundary-source offset can then be fixed separately. Until that
statement and the global differentiated stress/corrected geometry are
established, S139 and the full S074-S079 physical gates remain open.

Remaining to-do list: evaluate the common-scheme charge response with
the complete boundary prescription; quantify the global subtracted
stress and horizon-weighted metric response; solve the corrected
stationary exterior/state and transfer the driven results.
