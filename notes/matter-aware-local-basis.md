# Local operator reduction with charged matter retained

S074, 2026-10-13. The parity-even bosonic four-derivative action can be
reduced explicitly while retaining charged stress, current, cosmological
and boundary terms. In particular, eliminating the squared Maxwell
divergence generates a current-current interaction. Its complete
lowest-Landau-level angular spectrum is calculable, but a small angular
coefficient does not bound its unrenormalized continuum vacuum response.
This fixes a coefficient map needed for the common local action; it
does not close the independent fermion operators or renormalization.

## Local span and order

Use signature (-+++), kappa^2=8pi G, the canonical field F=Fbar, and

```
L0=(R-2Lambda)/(2kappa^2)-F^2/4+L_Dirac,
I=F_(mu nu) F^(mu nu),
X_(mu nu)=F_(mu rho) F_nu^rho,
S_(mu nu)=X_(mu nu)-g_(mu nu) I/4,
I4=F_(mu nu) F^(nu rho) F_(rho sigma) F^(sigma mu).
```

Let T be the remaining matter stress, tau=tr T, and define the canonical
matter current j so that D^nu:=nabla_mu F^(mu nu)=j^nu is the leading
Maxwell equation. For the charged Dirac field j contains its gauge
coupling g. Keep the following span with dimensionless coefficients:

```
L4=c1 R^2+c2 Ric^2+c3 Riem^2
   +kappa^2(c4 R I+c5 Ric.X+c6 Riem.FF+c9 D.D)
   +kappa^4(c7 I^2+c8 I4).                              (1)
```

The counting assigns weight two to curvature and weight one to F and
to a derivative. This is the bosonic four-derivative span, not a complete
operator basis with external Dirac fields. Coefficients and lower-order
parameters are understood at a common but still unspecified matching
scale. Every reduction below is first order in L4. Products of these
coefficients require the induced higher-order action.

[Cheung, Liu and Remmen, equation (4) and Appendix B](https://arxiv.org/html/1801.08546v3#A2)
give the eight nonderivative bosonic terms and their source-free metric
redefinition shifts. Our c4-c6 contain an extra kappa^2 and c7-c8 an
extra kappa^4 relative to their c notation. Their removal of D terms
assumes no charged sources. We retain c9 and perform the reduction
with the full matter equations. Their positivity theorem is not used.

For a closed two-form, the remaining derivative invariant satisfies

```
(nabla_rho F_mu_nu)^2 = 2 D.D-2 Ric.X+Riem.FF+nabla_mu B^mu,
B^mu=F^(alpha beta) nabla^mu F_(alpha beta)-2 F^(mu nu) D_nu.
                                                               (2)
```

Contracting the differentiated Bianchi identity gives
F.box F=2F^(mu nu)nabla_mu D_nu+2Ric.X-Riem.FF;
the product rule then gives (2). Other placements of two derivatives
reduce by integration by parts and that identity. Pure box R is also a
divergence. These are retained boundary terms, not zero in a driven
ensemble. The checker independently contracts three curved metric/closed
field families, including their boundary divergence, for all r>0.

## Complete coefficient map for the displayed span

Define

```
A=c1-c3, B=c2+4c3, U=c4-c6/3, V=c5+2c6,
u=c7-(B+V)/4, v=c8+B+V,
E4=Riem^2-4Ric^2+R^2,
WFF=Weyl_(mu nu rho sigma) F^(mu nu) F^(rho sigma).
```

The four-dimensional curvature identities give
Riem.FF=WFF+2Ric.X-R I/3. The leading Einstein equation specifies

```
R0=4Lambda-kappa^2 tau,
Ric0=Lambda g+kappa^2(S+T-g tau/2).
```

Using X.X=I4, S.S=I4-I^2/4 and tr S=0, the reduced density is

```
L4red = c3 E4+kappa^2 c6 WFF+kappa^4(u I^2+v I4)
      +(16A+4B)Lambda^2
      -(8A+2B)Lambda kappa^2 tau
      +(4U+V)Lambda kappa^2 I
      +kappa^4[A tau^2+B T.T+(2B+V)S.T-(U+V/4)tau I]
      +kappa^2 c9 j.j.                                  (3)
```

This retains both four-matter and mixed matter/electromagnetic contacts.
The Lambda-dependent terms also survive; absorbing them into lower-order
couplings requires matching those couplings and the physical metric.
Classical masslessness does not set T.T, S.T or j.j to zero. Further
Dirac field redefinitions could reduce some trace terms, but would need
their own contacts, sources and boundary terms; (3) makes no such step.

The two electromagnetic quartics may equivalently be written
kappa^4[(u+v/2)I^2+(v/4)J^2], where J=F.dualF. In the normalization of
[the finite electromagnetic input](finite-local-electromagnetic-input.md),
the J^2 coefficient is cJ=4v with the other reduced directions fixed.
The previous isolated J^2 term is reproduced by c7=-cJ/8, c8=cJ/4.

The map has an explicit off-shell verification. Put

```
E_(mu nu)=G_(mu nu)+Lambda g_(mu nu)-kappa^2(S+T)_(mu nu),
E=tr E_(mu nu), E_A=D-j,
H=A(-2R0+E)g+B(2Ric0-R0 g+E_tensor)
                      +kappa^2[V X-(U+V/2)I g].
```

Then the exact algebraic identity for (1) is

```
L4-L4red=E_(mu nu) H^(mu nu)+kappa^2 c9 E_A.(D+j).       (4)
```

Consequently the first-order change of variables
delta g^(mu nu)=-2kappa^2 H^(mu nu) and
delta A_nu=-kappa^2 c9(D+j)_nu removes this residual from L0+L4,
up to its boundary variation and terms quadratic in the local
coefficients. This uses the complete L0 equations, including matter.
It is a perturbative bulk redefinition, not an exact solution map for
a resummed higher-derivative theory.

The source-free Lambda=0 limit reproduces the four invariant
combinations c3, c6, u and v under the shifts in Cheung--Liu--Remmen
equation (95) at D=4. With charged matter, those four numbers alone
do not specify the action: the contacts in (3) must accompany the map.
If the original action already contains kappa^2 d_jj j.j, the current
coefficient after this step is d_jj+c9. Neither c9 in isolation nor
its removal is a physical matching prescription.

## Boundary data and quantum matching

The variational identity is

```
delta S0 = integral_bulk sqrt(-g)
 [E_(mu nu) delta g^(mu nu)/(2kappa^2)+E_A^nu delta A_nu
  +Dirac equations times delta spinors]
 +integral_boundary Theta0 +delta S0,boundary.
```

For the Maxwell part Theta0^mu=-F^(mu nu)delta A_nu. The gravitational,
spinor and chosen ensemble boundary variations must be carried as well.
Equation (2)'s flux is additional when that integration by parts is
used. E4 has no local bulk variation in four dimensions at fixed
topology, but its boundary completion and topology dependence cannot be
discarded. An electric source, physical clock or induced metric must
be transformed with the fields before observables are compared.

[Criado and Perez-Victoria, sections 2-3](https://arxiv.org/pdf/1811.09413)
explain why source transformations and the regularized Jacobian matter,
and why replacing operators with leading equations controls only first
order. We retain those requirements. This calculation does not set the
Jacobian to one or choose a continuum subtraction. Independent Dirac
operators, the boundary action, anomaly contacts and finite matching
conditions remain part of S074.

## The induced current contact in the sphere LLL

Use the retained sphere lowest Landau level with q states, spin
j_s=(q-1)/2, normalized sections eta_i and real orthonormal harmonics
Y_la. With psi4=b^(-1)sum_i eta_i chi_i, its two-dimensional current
is j4,a=g b^(-2)sum_ij eta_i^* eta_j (bar chi_i gamma_a chi_j).
The projected transverse current vanishes within this spin-polarized
level. At constant b, the contribution generated by c9 is

```
Delta S9=alpha9 integral sqrt(-h) sum_(l,a) j_Tla^b j_Tla,b,
alpha9=kappa^2 c9 g^2/b^2,
T_la=P_q Y_la P_q,
beta_l=(q-1)!q!/[(q-1-l)!(q+l)!], 0<=l<q.              (5)
```

Unlike angular photon exchange, (5) has no inverse l(l+1). It includes
l=0; this is an added local Wilson interaction and does not replace
or duplicate the retained monopole photon. Its physical coefficient
must be combined with any independently present j.j term.

The prior [angular normalization](weak-coupling-window.md) gives
Tr(T_la T_mb)=q beta_l delta_lm delta_ab/(4pi). Completeness gives

```
sum_la T_la^2=q/(4pi) I_q,
K_conj=sum_la T_la tensor conjugate(T_la),
spec K_conj: q beta_l/(4pi), multiplicity 2l+1,
||K_conj||=q/(4pi).                                    (6)
```

For a direct derivation, write R(Ω)_ij=eta_i^* eta_j. It is rank one.
Then R tensor conjugate(R)=|vec R><vec R|, and its harmonic expansion
gives K_conj=sum_la |vec T_la><vec T_la|. Orthogonality proves (6).
The unconjugated kernel is
q^2/[4pi(2q-1)] times the projector onto total spin q-1: the product
of two spin-j_s coherent states lies in that highest-spin space,
and its rotation average is the identity there divided by 2q-1.

The Legendre expansion of ((1+x)/2)^(q-1) and Parseval give the exact
second identity

```
sum_l(2l+1) beta_l=q,
sum_l(2l+1) beta_l^2=q^2/(2q-1).                        (7)
```

In the Euclidean current convention
S_int=(1/(2pi))integral sum G_AB J_A Jbar_B,
J_A=2pi :chi_R^dagger E_A chi_R:,
E_la=T_la/sqrt(q beta_l/(4pi)), the induced diagonal coupling is

```
G9_l=alpha9 q beta_l/(2pi^2),
||G9||_F^2=alpha9^2 q^4/[4pi^4(2q-1)].                 (8)
```

The factor follows from j_Lorentz.j_Lorentz=-4rho_R rho_L:
the Euclidean interaction is +4alpha9 sum rho_R,T rho_L,T.
This agrees with the preceding photon-exchange current convention
when alpha9 is replaced by g^2/[2l(l+1)].

In the retained family g=k^-4, q=k^4,
C=b^2/a^2=sqrt(1+12pi k^6), gamma=G/a^2=6C/[k^10(C-1)],
alpha9/c9=8pi gamma g^2/C. For bounded fixed c9,
the pair coefficient norm is O(k^-17) and ||G9||_F is O(k^-15).
These are finite angular/contact coefficients. They are not errors in
the physical vacuum, opening clock, energy or area.

For every nonzero resulting G9, the retained free-current cylinder
diagnostic has

```
||V_N Omega||^2=||G9||_F^2 N(N+1)(2N+1)/6,
||(H0+1)^(-1/2)V_N Omega||^2
                 =||G9||_F^2 sum_(n=1)^N n^2/(2n+1)
                 ~||G9||_F^2 N^2/4.                   (9)
```

The [existing current-vacuum proof](angular-current-renormalization.md)
applies also to the identity current, whose central term has the same
normalization. Thus arbitrarily small nonzero coefficients do not yield
a cutoff-uniform bare form estimate on the free vacuum. This excludes
that estimate, not a renormalized interacting state or the wormhole.

## Verification and next action

[The checker](../checks/matter_aware_local_basis.py) verifies the symbolic
map and its equation-of-motion residual. Eight independent algebraic
curvature/field/matter tensors each test all nine coefficient directions
and a mixed direction exactly. Three curved differential examples test
(2), and the source-free field-redefinition regression passes. Sixteen
exact flux cases verify the harmonic polynomial identities. Twelve
independent sphere quadratures, through q=8 with two resolutions, recover
both complete pair spectra. The continuum statements follow from the
displayed identities; these finite tests supplement them.

The reduced action must retain (3) and any independently matched
fermion terms. A prescription that drops the induced current contact
because the magnetic reference has j=0 is insufficient. Continue by
matching these composite operators with the other local/nonlocal terms
and physical boundary ensemble. The finite angular norm is useful
normalization data, not permission to omit their continuum contribution.

Remaining to-do list: complete the independent Dirac and boundary basis,
common regulator and finite matching, nonlinear even constraints and
remaining interactions, then S075-S080 state, geometry, stability,
transport, preparation/creation and final audit.
