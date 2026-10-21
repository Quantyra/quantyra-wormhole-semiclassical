# Uniform frequency bound for the curvature-free normal form

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The separated curvature-free normal form from
[the coupled outer calculation](outer-coupled-scalar-separation.md) is
strictly positive at every real scaled Euclidean frequency kappa, with
constants independent of kappa. For the normalized form R_kappa there,

    R_kappa[v] >= a ||H_kappa v||^2+b_g G_kappa[v],
    a=19/(150pi), b_g=1/(27g)+23/(4860pi),                 (1)
    H_kappa=-partial_eta^2+kappa^2,
    G_kappa=integral eta^-2(v'^2+kappa^2 v^2).

Its positive core operator has a closed Friedrichs realization. For
kappa>0 its inverse on unweighted L2(0,infinity) obeys

    ||R_kappa,F^-1|| <= 150pi/(19 kappa^4).                (2)

The full normal-action block is C R_kappa, so its inverse carries an
additional factor 1/C. This is a normal comparison inverse on the
positive frequency axis. It is not a charged-exterior causal inverse,
a claim about complex-frequency poles, or a physical boundary condition.
The separate fast scalar factor remains in the action.

## 1. Fixed form and notation

Keep the retained action, state, paired domain, subtraction and outer
scaling from the preceding note. No coefficient is varied. Initially
let v be real C_c^infinity(0,infinity), extended evenly when using
Fourier multipliers. Half of the full-line pairing is understood.
The normal formulas extend from the earlier profile interval [1,2]
to this core by the same dilation argument for any fixed compact
positive interval. Set dmu=eta^-4 d eta and

    D_k=eta^2(partial_eta^2-kappa^2)-2eta partial_eta,
    A_k=D_k-4, S_k=-A_k,
    B_k=eta^2(partial_eta^2-kappa^2)-12,
    T_k=(1/3)A_k^-1 B_k,
    P=integral eta^-4 v^2,
    Z=integral eta^-2 v^2, G_0=integral eta^-2 v'^2,
    G_k=G_0+kappa^2 Z,
    E_k=-<T_kv,A_kT_kv>_mu
       =(1/9)<B_kv,S_k^-1 B_kv>_mu.

S_k has the previously proved gap 25/4 and its closed form includes
the positive multiplication term kappa^2 eta^2. The shorthand k in
subscripts means kappa; it is not the benchmark integer 2000.

The new estimate that resolves the nonzero-frequency question is

    G_k/81 <= E_k-2P <= G_k/9, for all kappa>=0.            (3)

At kappa=0 the upper bound is an equality. Both constants in (3) are
uniform in profile support. No frequency-dependent potential is
treated as a uniformly bounded perturbation.

## 2. Resolvent representation and the exact zero-frequency tail

Define the first-order operator C_op=2eta partial_eta-8. Then

    B_k=-S_k+C_op,
    <v,C_op v>_mu=-5P,
    <v,S_k v>_mu=G_k+4P.

Expanding the square in the definition of E_k gives

    9 E_k=G_k+14P+J_k,
    J_k=<C_op v,S_k^-1 C_op v>_mu.                        (4)

Let F(eta)=integral_0^eta v(r)dr. The explicit stable inverse at zero
frequency gives

    u_0=S_0^-1 C_op v=-2F/eta,
    J_0=4P.                                              (5)

The eta^-1 tail in (5) is part of this identity. Omitting it would
change the resolvent comparison. The pairings are finite for compact
v, since F is zero near zero and constant sufficiently far out.

Because S_k=S_0+kappa^2 M with M the nonnegative multiplication
operator eta^2, positive-resolvent ordering gives

    0<=J_k<=J_0,
    0<=J_0-J_k<=kappa^2 <u_0,M u_0>_mu.                  (6)

This remains valid for unbounded M. One proof uses its closed quadratic
form in S_0^-1/2 M S_0^-1/2 and the scalar inequality
K(1+K)^-1<=K; alternatively minimize the positive quadratic energies
and test the S_k variational problem with u_0. The latter is permitted
because the right side of (6) is finite, as shown next.

Two weighted Hardy estimates suffice:

    G_0 >= (9/4)P,
    integral eta^-4 F^2 <= (4/9)Z.                       (7)

The first follows by expanding
integral eta^-2(v'-3v/(2eta))^2. For the second, integration by parts
gives integral eta^-4 F^2=(2/3) integral eta^-3 Fv, followed by
Cauchy-Schwarz. The boundary term vanishes at both ends, including
the constant tail of F. Hence

    <u_0,M u_0>_mu=4 integral eta^-4 F^2 <= (16/9)Z.

All arguments take place in the fixed auxiliary positive problem;
there is no inverse or boundary prescription for the fast scalar factor.

## 3. Uniform two-sided bound

From (4), positivity of J_k gives the first lower bound

    9(E_k-2P)>=G_0+kappa^2 Z-4P.                         (8)

Using J_k>=4P-(16/9)kappa^2 Z from (6)-(7) gives a second one:

    9(E_k-2P)>=G_0-(7/9)kappa^2 Z.                       (9)

Neither lower bound alone is sufficient at every frequency. Average
them and use P<=4G_0/9:

    9(E_k-2P)>=G_0-2P+(1/9)kappa^2 Z
                 >=(1/9)(G_0+kappa^2 Z).

This proves the lower half of (3). The upper half follows immediately
from J_k<=4P in (4). In particular, the nonlocal contribution E_k-2P
is a positive form equivalent to the local weighted Dirichlet norm.
It is not a freely adjustable local term.

## 4. The exact separated form and logarithmic positivity

The weighted curvature norm identity extends to every frequency:

    ||B_k v||_mu^2=||H_k v||^2+24G_k+72P.                 (10)

Insert (10) into the previously derived separation. The remaining
form is exactly

    R_k[v]=1/(30pi){<w,log(sqrt(H_k))w>
                    +integral log eta |w|^2
                    +(Lbar-6/5)||w||^2}
                +(3/g+23/(60pi))(E_k-2P),
    w=H_k v, Lbar=-log(mu0 ell_A)>9.                      (11)

Thus the zero-frequency cancellation was not an isolated cancellation
of the weighted potential. It leaves the same positive nonlocal
quantity (3) at arbitrary frequency.

On the even full line, log(sqrt(-partial_eta^2+kappa^2)) is at least
log|partial_eta| in Fourier calculus. The established
[Beckner logarithmic inequality](https://arxiv.org/pdf/math/0612486v2),
Corollary 2, equation (10), therefore implies

    <w,log(sqrt(H_k))w>+integral log eta |w|^2
       >= D_H||w||^2, D_H=log2+psi(1/4)>-4.

The same half-normalization applies to every term. Combining this with
Lbar>9 and (3) proves (1). More precisely the fourth coefficient is
strictly larger than 19/(150pi); the non-strict stated bound is useful
for completion. The sharper kappa=0 gradient coefficient in the
preceding note remains valid; the factor 1/81 is the uniform estimate.

## 5. Closed normal realization and the positive-axis inverse

Here is the domain argument, separate from the pointwise inequalities.
On the compact core, (11) is the quadratic form of a symmetric operator
in unweighted L2(0,infinity). The Fourier term is represented by the
even-line multiplier H_k^2 log(sqrt(H_k)); it maps the core to L2,
also at kappa=0 since p^4 log|p| is locally square integrable.
The weighted logarithmic term is H_k[(log eta)H_k v], which is smooth
and compact away from zero for a core vector. All displayed local
weighted terms likewise map the core into L2.

For the nonlocal term the only additional point is the endpoints of
S_k^-1 B_k v. At kappa>0 the homogeneous solutions are proportional to

    eta^(3/2) I_(5/2)(kappa eta),
    eta^(3/2) K_(5/2)(kappa eta).

The stable realization selects O(eta^4) at zero and O(eta exp(-kappa eta))
at infinity. These asymptotics also follow directly from the elementary
half-integer formulas, recorded in the checker. The unweighted adjoint
of B_k is

    B_k^*=eta^2 partial_eta^2+4eta partial_eta-10-kappa^2 eta^2.

Consequently B_k^*[eta^-4 S_k^-1 B_k v] is bounded near zero and
exponentially decaying at infinity. At kappa=0 the explicit tail from
the preceding note gives the same L2 conclusion (the outside term
cancels exactly). Thus the full symmetric core operator maps its
domain into L2. Symmetry follows from the positive auxiliary resolvent
and the fixed weighted pairings; it is not an assertion about a general
retarded stress kernel.

For clarity, semiboundedness here gives a closable form, rather than
an arbitrarily chosen completion. If v_n tends to zero in L2 and is
Cauchy in the positive form norm plus L2, its limiting form pairing
with any compact core test psi is
lim <v_n,R_k psi>=0. Density of the core in that form completion
forces its form limit to be zero. This proves closability. The closed
form defines the positive Friedrichs operator R_k,F. Its domain is
the closure selected by this comparison core; it is not imposed as
the charged spacetime's gravitational boundary domain.

For kappa>0, Fourier calculus gives ||H_k v||>=kappa^2||v||, so
R_k,F>=a kappa^4. Equivalently the Riesz representation of
v-><f,v> in the closed form space gives the unique weak solution
u=R_k,F^-1 f, proving (2) and the useful bounds

    ||H_k u|| <= ||f||/(a kappa^2),
    sqrt(G_k[u]) <= ||f||/(sqrt(a b_g) kappa^2).            (12)

These inequalities extend by closed-form approximation. They are normal
comparison estimates, not bounds for the full charged retarded inverse.

There is no uniform L2 gap as kappa->0. In fact the normal scaling is

    R_k[d^-1/2 v(eta/d)]=d^-4 R_(d k)[v].                 (13)

At k=0, letting d->infinity preserves L2 norm and drives the positive
form to zero, so the spectral bottom of R_0,F is zero. The divergence
kappa^-4 in (2) is consistent with this exact scaling. No zero-frequency
inverse or real-time decay theorem is inferred from a positive-axis bound.

## 6. Verification and next stability step

The companion checker verifies the weighted differential identities,
Hardy constants, two-bound combination, exact cancellation and inverse
constants. Independent Green-kernel quadratures evaluate E_k both from
B_k and from (4), retain the resolvent tail information, and check (3)
at several nonzero frequencies on two quadrature grids. The elementary
stable Green solutions and endpoint action are checked symbolically.
These finite fixtures do not compute an actual exterior spectrum.
The continuum inverse/domain statement is an ordinary proof, not formal
certification. No Lean invocation is selected.

The remaining scalar normal factor is still D_k+e with
e=60pi/g+11/3. This work makes the curvature-free normal block available
for a controlled coupled reduction. It does not by itself transfer the
AdS fast branch to a charged exterior: variable-curvature terms,
initial/incoming Ward constraints, Maxwell/current response, the actual
state and the physical boundary realization must still be included.
Analytic continuation and possible complex-frequency poles also require
their own argument. No nonlinear background or physical stability
tolerance is established. S168.3-S168.4 and all S167 tasks remain open.

Remaining to-do list: actual charged retarded constraints and coupled
transfer or justified exclusion at an admitted physical scale; signed
coefficient, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, stability, transport, preparation,
separate creation and full audit.
