# Coupled outer limit and separation of the curvature scalar

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The outer limit of the retained action separates into a curvature-scalar
form and a form with zero linearized scalar curvature. This separation
uses only the inverse of a negative auxiliary factor with a spectral gap;
it does not invert or remove the fast scalar factor. The full direct
shape block is unbounded below on outer annuli at each fixed thermal
frequency, despite the canonical logarithmic bound. That direct-block
sign is not a physical instability verdict. At zero scaled frequency,
the curvature-free block after separation is strictly positive.

These results determine which block can be used in the next coupled
reduction. They do not impose the actual charged exterior's retarded
constraints or supply a nonlinear background. S168.3-S168.4 remain open.

## 1. Coordinates, scaling and retained data

Keep the centered seed, paired HH domain, common subtraction and all
coefficients from [the preceding comparison](uniform-logarithmic-frequency-comparison.md).
Write

    g=G lambda, x=lambda B^2, C=B^2/ell_o^2,
    alpha=C/x, ell_A=1/sqrt(lambda),
    Lbar=-log(mu0 ell_A), lambda0/lambda=1+11g/(360pi).

Here g is the dimensionless coupling, not the metric tensor. Normalize
the canonical outer time coordinate by its constant seam warp. The
canonical orbit is gbar=dY^2+a(Y)^2 d tau^2 with a(0)=1 and

    h=log a=O(Y^2), Phi=log(r/B)=-log Y+.5log(alpha)+O(Y^2).

Indeed u=B/r satisfies u'=sqrt(calB(u)/C), calB(0)=x,
calB'(0)=0, so u=sqrt(x/C)Y+O(Y^3). The same expansions hold with
the differentiated orders needed below on the actual polynomial seed.
The fixed charge and sphere-curvature terms are lower order in this
limit. They have not been removed from the original action.

Take a compact profile on 1<eta<2 and a nonzero thermal mode,

    v_d(Y,tau)=d^-1/2 v(Y/d)cos(nu_d tau), eta=Y/d,
    nu_d d -> kappa>=0.

Divide every form by the actual time norm beta_tau/2. A fixed nonzero
Matsubara frequency is allowed and corresponds to kappa=0. Frequencies
are in the displayed canonical coordinates, not unconverted proper
physical units. Extend profiles evenly to the complete line whenever
using scalar Fourier calculus, and use half of full-line pairings.

## 2. Exact ultraviolet constant of the canonical multiplier

For fixed angular data, the existing flat multiplier has asymptotic

    Kflat(p)=C/(30pi) p^4[log(p/(mu0 ell_o))-31/30]+o(p^4). (1)

This is an asymptotic coefficient in the prescribed scheme. It does
not change a finite counterterm or replace an all-frequency estimate.
To derive it, use the established

    F(t)=-2[psi(1+nu_+)+psi(1+nu_-)],
    nu_+/-=(q+/-sqrt(q^2-4t))/2.

The sectorial [digamma expansion](https://dlmf.nist.gov/5.11.E2) gives
F(t)=-2log t+o(1), t F'(t)=-2+o(1), for fixed q. The roots with
their added 1 remain in the permitted right-half-plane sector; the
coincident-root point is removable in their symmetric sum. Uniform
analytic asymptotics also justify differentiating the remainder.
The existing Laplace/count bounds dominate the alpha endpoints by
an integrable multiple of a_alpha^2(1+|log a_alpha|), where
a_alpha=alpha(1-alpha). With K(z)=Kflat(sqrt z),

    K''(z)=-C/(2pi) integral_0^1 a_alpha^2
                [F(C a_alpha z)+C a_alpha z F'(C a_alpha z)
                                      +4log(mu0 B)]d alpha
      =C/(30pi)[log(z/(mu0^2 ell_o^2))-17/30]+o(1).

Here integral a_alpha^2=1/30 and integral a_alpha^2 log a_alpha=-47/900.
Integrating twice subtracts 3/2 inside the bracket. Thus
K(z)=C z^2/(60pi)[log(z/(mu0^2 ell_o^2))-31/15]+o(z^2),
which is (1). An exact useful alternative, with K(0)=K'(0)=0, is

    K(z)/z^2=-C/(2pi) integral_0^1 a_alpha^2
                       integral_0^1 u[F(C a_alpha z u)+4log(mu0 B)]du d alpha.

No numerical ultraviolet crossover for the actual large q is claimed.

## 3. Passing from the actual seed form to the dilation limit

The two-dimensional gain-two estimates used in the preceding notes
can take ANY positive fourth-order loss epsilon, with finite lower-order
constants depending on epsilon. The actual constant contact retains
its 1/Y bound; the zero-trace inequality controls it. Thus, for every
epsilon>0, the actual canonical form and auxiliary scalar multiplier
differ by at most

    C[epsilon J2+A_epsilon J1+B_epsilon J0].

On the displayed annuli J1/J2=O(d^2), J0/J2=O(d^4). Taking d->0
first and then epsilon->0 shows that this difference is o(C d^-4).
Using only the previously rounded .001 bound would not justify this
exact limit; the arbitrary-epsilon gain-two estimate is essential.

For the auxiliary scalar, after half-density conjugation and dilation,

    H_d=-partial_eta^2+d^2 V(d eta)+(nu_d d)^2 exp(-2h(d eta)),
    V=h''/2+h'^2/4,
    H_kappa=-partial_eta^2+kappa^2.

H_d converges strongly in resolvent sense to H_kappa on the even line.
This follows directly from bounded multiplication perturbations of
-partial_eta^2: the potentials are uniformly bounded and converge
pointwise, so their resolvents converge strongly by the resolvent
identity. The half-density profile sqrt(a(d eta))v tends to v in H3.
The potentials and their first weak derivatives are uniformly bounded.
For w_d=H_d[sqrt(a(d eta))v], both its H1 norm and the third spectral
moment of the original profile are therefore uniformly bounded.
This gives uniform integrability for the spectral function
z^2 log(sqrt z); at zero that function extends continuously by zero.
Strong resolvent convergence consequently passes its quadratic forms
to the limit. No fourth seam derivative is assumed.

Equation (1) can be passed through the same spectral measures. For each
epsilon its remainder is bounded by epsilon p^4 plus a finite constant
on the complementary compact momentum interval. After multiplying by
d^4 the latter constant vanishes. The term -log d multiplies the EXACT
J2 norm on both sides of the canonical/Weyl cancellation; it is cancelled
before taking the limit. This avoids an unjustified product of a small
norm error with the divergent logarithm.

The nonlocal limiting logarithm is therefore log(sqrt(H_kappa)), with
the even-line realization inherited from the actual doubled core. It
is not a freely chosen boundary operator or a replacement quantum state.

## 4. Full direct shape block in the outer limit

Put w=H_kappa v. The literal V,V block of the FULL retained Hessian
has the following limit, divided by C and with the time norm removed:

    lim_(d->0) d^4 Q_VV[v_d]/C = F_kappa[v],
    F_kappa[v]=1/(30pi){<w,log(sqrt(H_kappa))w>
                 +integral log eta |w|^2+(Lbar-31/30)||w||^2}
          +2/(15pi) integral eta^-2(v'^2+kappa^2 v^2)
          -(6/g+11/(30pi)) integral eta^-4 v^2.            (2)

All integrals on positive eta use physical-half normalization. The
Fourier logarithm is the half of its even full-line pairing.

The local terms in (2) can be checked without a kernel approximation.
In the leading flat canonical chart, Phi=-log Y+.5log(alpha). The
classical direct shape Hessian is

    -6C/g (lambda0/lambda) integral Y^-4 v^2.

Expanding the exact Weyl density at fixed Phi gives

    Q_W,VV=-C/(30pi) integral Phi(Delta v)^2
            -2C/(15pi) integral Y^-2 v Delta v
            +13C/(60pi) integral Y^-4 v^2
          =-C/(30pi) integral Phi(Delta v)^2
            +2C/(15pi) integral Y^-2 |grad v|^2
            -11C/(60pi) integral Y^-4 v^2.

The cosmological correction supplies another -11C/(60pi) weighted
potential. The canonical logarithm cancels -log d in Phi, and
ell_o sqrt(alpha)=ell_A fixes the remaining scale Lbar. Bounded orbit
curvature, sphere, magnetic and regular finite-contact terms vanish
after this dilation. Their global returns were already retained in
the comparison of section 3. This proves (2) for the actual seed,
including the literal off-shell chart convention of its Hessian.

There is a concrete obstruction to treating Q_VV as a positive block.
For kappa=0 choose the H3 profile v=(eta-1)^3(2-eta)^3 on [1,2],
zero elsewhere, with smooth compact approximations. Its derivative
norms Ij=integral |partial_eta^j v|^2 satisfy

    I0=1/12012, I1/I0=78/5, I2/I0=3432/5, I3/I0=61776.

The retained benchmark has 9<Lbar<12. For the upper bound, C/q<.00352,
q=2000^4, C=1+6x+d and |d|<1e-12 give x<1e10; combine this with
|log(mu0 B)|<1e-21. Using log|D|<=D^2/2 by scalar Fourier calculus,
log eta<=log2<.7, pi>3 and eta^-4>=1/16 on the profile support,

    F_0[v]/I0 < 97396/225-3/(8g) < -1/(4g).               (3)

The last inequality follows from g<(4e-17)^2. The logarithmic form is
continuous under H3 compact approximation, so a smooth profile retains
a strict negative margin. At every FIXED nonzero thermal frequency,
nu d->0. Its actual seed direct-shape Rayleigh quotients in the
canonical orbit L2 norm therefore tend to minus infinity as d^-4.
Any self-adjoint realization of that block containing this core, if
such a realization is chosen, cannot be semibounded below.

This is a statement about a coordinate block. It is not a physical
constraint-satisfying disturbance, a gravitational energy or a charged
wormhole instability theorem. It rules out using positivity of Qcan
to infer a positive inverse of the full direct Q_VV block.

## 5. Coupled normal action and a stable scalar separation

The leading metric is the Poincare AdS normal geometry. On radial
Fourier profiles define, in the measure dmu=eta^-4 d eta,

    D_kappa=eta^2(partial_eta^2-kappa^2)-2eta partial_eta,
    A_kappa=D_kappa-4,
    B_kappa=eta^2(partial_eta^2-kappa^2)-12,
    e=E/lambda=60pi/g+11/3.

For a four-dimensional conformal amplitude p and orbit shape v,
the linearized scalar curvature is

    delta R=-2lambda(3 A_kappa p+B_kappa v).               (4)

The complete p,p and p,v columns are local. Direct expansion of the
retained action, or its exact conformal-column identity at the AdS
normal geometry, gives

    Q_pp/C=1/(20pi)<p,(D_kappa+e)A_kappa p>_mu,
    Q_pv/C=1/(60pi)<p,(D_kappa+e)B_kappa v>_mu,            (5)
    Q_normal=Q_pp+2Q_pv+C F_kappa[v].

The leading local first variation vanishes with the prescribed
lambda0/lambda; no seed residual is dropped before taking its normal
limit. The canonical functional has no p dependence. Formula (5)
does not rely on a self-adjoint retarded kernel; it concerns the local
columns and the Euclidean normal form.

Use the negative self-adjoint realization of A_kappa defined by its
closed Dirichlet form in L2(dmu). With eta=exp(t) and
p=exp(3t/2)u(t), its unitarily transformed expression is

    partial_t^2-25/4-kappa^2 exp(2t).

Hence A_kappa<=-25/4 and ||A_kappa^-1||<=4/25. Define

    T_kappa v=(1/3) A_kappa^-1 B_kappa v,
    s=p+T_kappa v.

This auxiliary inverse is well defined on the compact core. It does
not invert D_kappa+e, select a phase for that fast factor or impose
physical gravitational boundary data. Equation (4) becomes
delta R=-6lambda A_kappa s. Since A_kappa commutes with D_kappa+e,
the normal action separates exactly:

    Q_normal/C=1/(20pi)<s,(D_kappa+e)A_kappa s>_mu
                       +R_kappa[v],
    R_kappa[v]=F_kappa[v]
       -1/(20pi)<T_kappa v,(D_kappa+e)A_kappa T_kappa v>_mu. (6)

The curvature-free component has p=-T_kappa v. It has not been
identified with a fully constraint-reduced physical degree of freedom.
This is a legitimate normal-variable separation that avoids the
unsupported positive-inverse assumption refuted by (3).

The scalar factor remains D_kappa+e. In the exact AdS comparison,
solutions of it lift through the previously established tensor
-2epsilon Q_ab(Box+E)s, satisfying the tensor constraints there.
That comparison lift is not new evidence for the actual charged
geometry. No branch has been eliminated by (6).

## 6. Positive curvature-free block at zero scaled frequency

At kappa=0 the auxiliary inverse is explicit:

    (A_0^-1 f)(eta)=-1/5[eta^-1 integral_0^eta f(r)dr
                         +eta^4 integral_eta^infinity f(r)r^-5 dr],
    T_0 v=v/3+2/(3eta) integral_0^eta v(r)dr.              (7)

The tail of T_0 v is proportional to eta^-1 beyond the profile; it
must be kept. Its stable form norm is finite. Integration by parts gives

    -<T_0v,A_0T_0v>_mu=(1/9) integral eta^-2 v'^2
                                      +2 integral eta^-4 v^2,
    ||B_0v||_mu^2=||v''||^2+24 integral eta^-2 v'^2
                                      +72 integral eta^-4 v^2.

Substitution into (6) cancels the negative eta^-4 potential exactly:

    R_0[v]=1/(30pi){<v'',log|D|v''>
                    +integral log eta |v''|^2+(Lbar-6/5)||v''||^2}
              +(1/(3g)+23/(540pi)) integral eta^-2 v'^2.  (8)

The established even-line logarithmic inequality has constant D_H>-4.
With Lbar>9, (8) implies the explicit positive lower bound

    R_0[v] > 19/(150pi)||v''||^2
               +(1/(3g)+23/(540pi)) integral eta^-2 v'^2 (9)

for nonzero compact profiles. Thus the direct-block negative directions
and this positive curvature-free normal block are consistent. The
change of variables separates their conformal content. It is not an
order reduction and does not discard the scalar fast factor.

## 7. Consequence and remaining physical problem

Use (6), rather than a presumed positive direct-shape inverse, in the
next outer coupled analysis. Formula (9) proves only the kappa=0
normal curvature-free bound. The form R_kappa for nonzero scaled
frequency is specified by (2),(6), but its global constrained retarded
realization and transfer estimates have not been proved here.

The auxiliary A inverse is not a solution of the charged initial or
incoming Ward constraints. A variable-curvature exterior has the
gauge-invariant scalar delta R-R' delta r and the additional couplings
already derived in the charged conformal note. Those couplings,
Maxwell/current response, boundary domains and the actual state must
be treated. No exact nonlinear exterior, stability tolerance, admitted
ultraviolet prescription or physical fast-branch exclusion is obtained.
S167's signed coefficient, higher response, nonlinear solve and
compatibility remain open.

The checker verifies exact local-action columns and cancellation,
the fixed ultraviolet constant, auxiliary inverse and weighted
identities, the direct-block negative margin and separated positive
margin. Numerical digamma checks are finite fixtures, not an evaluated
actual spectrum or an effective-theory validity test. Continuum limits
and domain statements remain ordinary arguments; no Lean run is selected.

Remaining to-do list: nonzero-frequency curvature-free control and
actual retarded constraints; transfer or justified exclusion at an
admitted physical scale; signed coefficient, cold/higher response,
nonlinear solve and compatibility; exterior/state, charge/return, opening,
stability, transport, preparation, separate creation and full audit.
