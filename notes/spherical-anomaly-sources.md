# Spherical anomaly sources and consistent nonlinear constraints

S074, 2026-10-13. The accepted spherical anomaly functional now gives
its separate metric stress, radius source and electric current, with
their full bulk variational boundary current. They enter a nonlinear
Einstein--Maxwell constraint system with an explicit propagation
identity. This also corrects two inherited source-convention errors:
the sign of the angular source equation and the normalization of the
electric charge flux. The independent state functional, physical mouth
ensemble and other required interactions remain unfinished.

## A form suitable for all field variations

Use the conventions of [the anomaly action](spherical-anomaly-action.md):
X=b0^2 exp(2sigma), g2=exp(2sigma)h, y=b0^-2, r=R[h],
u=Box_h sigma, v=(nabla_h sigma)^2 and I=F2^2+q^2 y^2/2.
Let a=11/360, c=1/20, retain arbitrary b, and d=b-2a/3.
The local heat representative has b=1/30. All derivatives in this
section use h. Set N=b0^2/(4pi). The same functional can be written

```
W=N integral sqrt(-h) [L+div Bstar],
L=sigma[c(r+2y)^2/3-4ayr-2I/3]
  +b r u-3b u^2+[br+(2b-4a)y]v-6d u v-3d v^2,
Bstar^a=B^a+b(sigma nabla^a r-r nabla^a sigma).          (1)
```

B is the retained flux in the preceding note, equation (4). The
additional divergence exchanges sigma Box r for r Box sigma; it is
not discarded at a mouth. Define the partial derivatives of L by

```
P=L_r=sigma[2c(r+2y)/3-4ay]+b(u+v),
Q=L_u=br-6bu-6dv,
V=L_v=br+(2b-4a)y-6d(u+v), Z=L_I=-2sigma/3,
L_sigma=c(r+2y)^2/3-4ayr-2I/3,
E_sigma=L_sigma+Box Q-2nabla_a(V nabla^a sigma).        (2)
```

This E_sigma equals the previous Euler derivative. In particular,
E_sigma=16pi^2 exp(4sigma) A4[g4]. The transgression is established
bulk prior work: [Mazur--Mottola, equations (44), (50) and (56)](https://arxiv.org/html/hep-th/0106151v1).
The expressions below vary its spherical representative; no claim
of a new general anomaly construction or boundary ensemble is made.

## Separate stress and current

At fixed sigma and A, write delta W_bulk=N integral sqrt(-h)
[M_ab delta h^ab+E_sigma delta sigma+J^a delta A_a+div Theta].
Then

```
M_ab=-h_ab L/2+(r/2)h_ab P+(h_ab Box-nabla_a nabla_b)P
     +V sigma_a sigma_b-sigma_(a nabla_b)Q
     +h_ab[(nabla Q).(nabla sigma)+Q u]/2
     +2Z F_ac F_b^c,
tau_ab^W=-2N M_ab, s_W=N E_sigma,
j_h,W^b=-4N nabla_a(Z F^ab).                          (3)
```

Parentheses symmetrize with weight one half. The Q terms come from
varying Box sigma before integration by parts. In particular they
cannot be recovered by varying only the explicit inverse metrics.
The magnetic part of I contributes to L and E_sigma, but is independent
of h; only F2 appears in the last term of M. These sources obey

```
nabla^a tau_ab^W=s_W sigma_b+F_ba j_h,W^a,
nabla_a j_h,W^a=0.                                    (4)
```

For k^ab=delta h^ab, k=h_ab k^ab, z=delta sigma and alpha_b=delta A_b,
the boundary current for the displayed L is

```
Theta^a=P(nabla^a k-nabla_b k^ab)-k nabla^a P+k^ab nabla_b P
       +Q sigma_b k^ab-Q sigma^a k/2
       +Q nabla^a z-z nabla^a Q+2V sigma^a z
       +4Z F^ab alpha_b.                              (5)
```

The full density variation also includes N partial_a delta(sqrt(-h)
Bstar^a). Formula (5) alone therefore is not the complete variation
of W on an open boundary. A joint metric/radius/Maxwell boundary
action and its anomaly must be specified for the physical ensemble.

## Return to the physical metric and radius

Write Gamma_D=W+Gamma_hat[h,A] for the same conformally related
state/domain data as in the preceding note. Let tau_ab and j_h^a include
the respective derivatives of Gamma_hat; s=s_W because Gamma_hat has
no independent sigma argument in this representation. A variation at
fixed physical g changes h as well as sigma. The exact chain rule is

```
T_ab^(2)=tau_ab,
J_X=(2G/X) exp(-2sigma)[s-h^ab tau_ab],
j_g^a=exp(-2sigma)j_h^a.                              (6)
```

Thus the anomaly's separate contribution is explicit by inserting
(3); the full J_X also contains the trace of the unknown Gamma_hat
stress. Neither the h-frame stress nor its trace may be set to zero
when transferring back to g. Equations (4) and (6) give

```
T^(2)a_a+X J_X/(2G)=4pi X A4,
nabla_g^a T_ab^(2)=J_X partial_b X/(4G)+F_ba j_g^a.     (7)
```

Extra finite local terms have their own variations; they are not absorbed
by pretending (7)'s Dirac anomaly is the total EFT trace. The transformation
does not determine Gamma_hat or the physical initial state. Its utility
is that the entire local anomalous contribution to all three equations
can now be included consistently, not just their trace combination.

## Coupled metric, angular and Maxwell equations

Use the physical two-metric g, an unrescaled gauge field and coupling g_c.
Include a cosmological constant, fixed magnetic q and arbitrary electric
F_ab. With r_m^2=pi G q^2/g_c^2, put

```
V0(X)=2-2Lambda X-2r_m^2/X,
I0=1/(4G) integral sqrt(-g)[X R+(nabla X)^2/(2X)+V0(X)]
   -pi/g_c^2 integral sqrt(-g) X F2^2.                 (8)
```

Let Gamma contain the Dirac, finite local and other retained sources,
excluding the Maxwell term already in I0. Use T=-2 delta Gamma/delta
g^ab/sqrt(-g), J_X=4G delta Gamma/delta X/sqrt(-g), and
j=delta Gamma/delta A/sqrt(-g). Define

```
E_ab=g_ab Box X-nabla_a nabla_b X
     +(X_a X_b-g_ab(nabla X)^2/2)/(2X)-g_ab V0/2,
D0=R-Box X/X+(nabla X)^2/(2X^2)+V0'(X),
T_ab^M=(4pi X/g_c^2)(F_ac F_b^c-g_ab F2^2/4),
C_ab=E_ab-2G(T_ab+T_ab^M),
C_X=D0-(4pi G/g_c^2)F2^2+J_X,
C_A^b=(4pi/g_c^2)nabla_a(X F^ab)+j^b.                 (9)
```

The equations are C_ab=C_X=C_A=0. In particular, without the electric
term the angular equation is D0=-J_X. This sign follows immediately
from varying I0+Gamma and is independently forced by propagation:

```
nabla^a C_ab=-C_X X_b/2-2G F_ba C_A^a.                 (10)
```

If the retained source has Ward defect
R_b=nabla^a T_ab-J_X X_b/(4G)-F_ba j^a, the right side of (10)
has the additional term -2G R_b. Source conservation is therefore a
necessary input to the constraints, not a check that can be postponed.
These are the standard dilaton-gravity/Noether mechanisms; see
[Grumiller--Kummer--Vassilevich, (2.9) and section 5.2](https://arxiv.org/html/hep-th/0204253v3).
Their sign/normalization conventions are not substituted for the direct
variation of (8).

For ds2^2=-exp(2rho) du dv, the two null constraints and the mixed equation
are

```
C_uu=-X_uu+2rho_u X_u+X_u^2/(2X)-2G T_uu,
C_vv=-X_vv+2rho_v X_v+X_v^2/(2X)-2G T_vv,
C_uv=X_uv+exp(2rho)V0/4-2G(T_uv+T_uv^M).              (11)
```

The classical electric Maxwell stress has zero uu and vv components.
If C_uv=C_X=C_A=R_b=0, (10) implies partial_v C_uu=0 and
partial_u C_vv=0. Hence imposing the null constraints on the two
initial characteristics propagates them throughout a smooth patch.
This does not prove existence or uniqueness of the full semiclassical
initial-boundary problem. It identifies the equations and compatible
source data that such a solution needs.

With R=sqrt(X), the exact null focusing equation is
partial_u[exp(-2rho) R_u]=-G exp(-2rho)T_uu/R, and likewise in v.
It connects the separate null stress in (3),(6) to area evolution.
The anomaly trace alone does not determine either focusing integral.
Fourth derivatives in local EFT sources must be treated at the declared
perturbative order with matching/order reduction; the displayed system
does not add independently specifiable higher-derivative initial data.

## Charge normalization and correction of inherited formulae

Choose epsilon_uv=+sqrt(-g), define F_ab=f epsilon_ab explicitly, and
let Q_e=(4pi/g_c^2)X f. Then the Maxwell equation in (9) gives

```
partial_a Q_e=epsilon_ab j^b,
partial_u Q_e=-j_u, partial_v Q_e=+j_v.
For D=X f: partial_a D=(g_c^2/(4pi))epsilon_ab j^b.     (12)
```

Reversing orientation reverses both f and Q_e. The factors for Q_e
and D are different at general coupling, even though they coincide
in units g_c^2/(4pi)=1.

For the same charge definition, the geometric charged mass aspect is

```
m=R/(2G)[1-(nabla R)^2-Lambda R^2/3
         +(r_m^2+G g_c^2 Q_e^2/(4pi))/R^2],
partial_a m=(T_a^b-delta_a^b T)partial_b R
             +g_c^2 Q_e partial_a Q_e/(4pi R).         (13)
```

Here T excludes the classical electric Maxwell stress already removed
by the charge term in m. The metric equation in (9) proves (13) on
arbitrary radius and charge jets; the angular equation is not needed
for that algebraic identity. Source consistency and (10) are needed
for an actual solution. The last term is the electric work, which
would be misnormalized by using D's derivative for Q_e. This extends
the [retained magnetic mass balance](sourced-spherical-mass-charge-balance.md)
to Lambda and variable electric charge in one consistent convention.

At predecessor commit dbab9e409a8f7b23106869af9798e2da622fa983,
`exact-spherical-em-reduction.md` (3) prints +J_X on the right despite
defining J_X with a positive functional derivative. It is corrected
to -J_X. The simplest regression is Gamma=integral sqrt(-g) f(X):
the effective potential is V0+4G f, so D0+4G f'=0, not D0-4G f'=0.
The previous metric/JT check had no nonzero independent J_X and could
not expose this error.

The same predecessor's `sourced-spherical-mass-charge-balance.md` (2)
and `retarded-dirac-spherical-response.md` (2) attach D's factor to
the scaled Q_e. Their current text is corrected. The historical
`check_current_and_gauss` calculation in
`checks/retarded_dirac_spherical_response.py` verifies the unscaled D
law and chiral integrability; its Q label must be read as D. Its mouth
charge quadratures explicitly set g_c^2/(4pi)=1 and cannot test this
normalization distinction. The separate Schwinger-mass calculation
starts from the unscaled Maxwell equation and is unchanged.

The frozen clock manuscript states the Ward identity and qualitative
flux split, but does not print either erroneous equation. Its mode-count
argument uses the source-free JT modes, mass variation and normalized
clocks; these two corrections do not change that argument. Historical
receipts and frozen artifacts are preserved, with their scope correction
recorded here. No publication or prior workstream is reopened.

## Verification and remaining scope

The checker derives the stress and source identities, compares direct
action variations including independent metric components, and tests
the corrected gravitational and Gauss conventions at general coupling.
It also checks the physical-frame Ward transformation and charged
mass-work identity independently of the anomaly coefficient choices.
The accompanying receipt records exact and supplementary numerical
checks separately. The full boundary current is retained in (5), not
used as a substitute for a physical boundary action.

This supplies the local anomaly's separate nonlinear sources and fixes
the consistent spherical constraint system needed by S074 and S076.
It does not supply the independent nonlocal determinant, full angular
interactions, finite physical matching, physical state or mouth data.
No complete opening, stability, transport or creation claim follows.

Remaining to-do list: complete S074's other interactions, finite and
boundary matching and controlled remainder; S075 state/drive; S076
stress/geometry; S077 stability; S078 transport; S079 preparation/creation;
S080 final audit.
