# Fast scalar boundary coefficients and dynamical solvability

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

Removing both fast scalar asymptotic coefficients is not an ordinary
reflecting boundary condition for the retained dynamical scalar branch.
In the exact AdS comparison it imposes one nontrivial compatibility
condition on the forcing at EACH Laplace frequency. This condition
excludes arbitrarily small compact initial velocities. Therefore that
boundary prescription cannot by itself give an exponentially bounded
retarded evolution on an open set of the full scalar initial-data space.

This is a precise test of a proposed boundary exclusion, not a transfer
of instability to an actual charged solution. It does not invalidate
the stationary weighted range construction, which already retains its
source compatibility. The actual gravitational boundary realization
and the physical disturbance class still have to be established.

## 1. Retained branch and the stationary boundary question

Use the [exact vacuum comparison](ads-fast-scalar-stability-preflight.md),
including its fixed local coefficients and conformal-vacuum metric lift.
Its gauge-invariant curvature perturbation psi=delta R satisfies

    (Box+E)psi=0, E=60pi/G+11lambda/3,
    psi=6(E+4lambda)s, h_ab=2s g_ab.                    (1)

The conformal lift in (1) satisfies every linearized tensor component
in THAT comparison. It is not a charged-background lift. With global
AdS boundary coordinate z, center L=pi/2, and u=cot(z)psi, its equation is

    u_tt+A u=0,
    A=-partial_z^2-a/sin(z)^2,
    a=E/lambda-2, k=sqrt(a-1/4)>0, 0<z<L.              (2)

The Hilbert measure is dz and regularity at the center is u(L)=0.
The time t in this note is the dimensionless global AdS time.
For a homogeneous equation (A+s_L^2)u=0, s_L>0 is a real Laplace
parameter, not the conformal amplitude s in (1). The checker calls it gamma.

The actual stationary construction uses the fast remainder space
X_(7/2), plus its physical mass column, and excludes both leading
z^(3/2) scalar oscillations; see the
[weighted map](stationary-global-weighted-map.md) and
[adjoint phase](stationary-adjoint-phase.md). The additional signed
range condition is already part of that construction. The question
here is whether imposing the same removal on dynamical scalar data
would automatically remove the fast branch. It does not.

## 2. Boundary traces and flux

Since sin(z)^-2=z^-2+1/3+O(z^2), the homogeneous radial equation is
regular singular. Its two independent real solutions, normalized by
their leading terms, have the form

    phi_c(z)=sqrt(z)[cos(k log z)+O(z^2)],
    phi_s(z)=sqrt(z)[sin(k log z)+O(z^2)].               (3)

The O(z^2) notation allows bounded sine/cosine combinations and their
conormal derivatives, with constants depending on the fixed a,s_L.
For example complex Frobenius exponents 1/2+/-ik have noninteger
difference and give (3) by the usual convergent power-series recurrence.
Their Wronskian phi_c phi_s'-phi_c' phi_s is exactly k: it is constant
for this equation and its leading value is k. Thus a homogeneous
near-boundary field u=A_c phi_c+A_s phi_s has two independent traces.
For two such fields the limiting Hermitian Green form is

    lim_(z->0)(conj(u)v'-conj(u')v)
       =k(conj(A_c)B_s-conj(A_s)B_c).                  (4)

Both branches are L2(dz). With the regular center fixed, a self-adjoint
separated realization selects one maximal isotropic complex line in
this two-dimensional trace space. It can be written

    cos(theta) A_c+sin(theta) A_s=0, theta real mod pi. (5)

The condition A_c=A_s=0 is the minimal trace domain, not a maximal
isotropic line; its adjoint permits both traces. It is therefore not
a self-adjoint boundary realization. This conclusion is about this
second-order curvature branch. No self-adjointness of a general
retarded quantum metric kernel is assumed.

The established below-BF result remains relevant: every self-adjoint
extension of (2) is unbounded below. It is already proved with compact
tests in the vacuum preflight and is the general result of
[Ishibashi--Wald, Proposition 3.1 and Theorem 3.2(iv)](https://arxiv.org/pdf/hep-th/0402184v4).
Choosing a different phase in (5) does not restore semiboundedness.
The new solvability argument below also treats the stronger two-zero-
trace prescription without assuming it is self-adjoint.

Fixed leading conformal boundary metric does not itself remove (3).
In the explicit conformal lift, psi=tan(z)u and s is proportional to
psi, so both modes vanish as z^(3/2). Their first z derivatives vanish,
but a nonzero mode has an unbounded oscillating second derivative.
Thus a C2 condition on this compactified conformal representative, or
the stronger stationary fast-decay condition, removes BOTH amplitudes.
This statement concerns the displayed conformal chart; it is not an
assertion that an arbitrary non-smooth coordinate change is admissible
or that the physical gravitational boundary class has been selected.

## 3. Exact compatibility for the stronger boundary prescription

Fix any real s_L>0. Define the unique real homogeneous solution w_s by

    (A+s_L^2)w_s=0, w_s(L)=0, w_s'(L)=1.              (6)

It exists on the open interval, has the traces (3), is in L2, and is
not identically zero. Let f be smooth with compact support in (0,L).
Consider the inhomogeneous problem

    (A+s_L^2)U=f, U(L)=0, A_c(U)=A_s(U)=0.             (7)

Because f vanishes near z=0, (7) is homogeneous there. Both zero traces
therefore force U to vanish identically on a boundary neighborhood.
Starting with U=U'=0 to the left of supp(f), ordinary ODE evolution
gives a unique U_0 throughout the remaining interval. The Green
identity, with no boundary contribution at zero, gives exactly

    M_s(f):=integral_0^L w_s(z) f(z)dz=U_0(L).          (8)

Indeed w_s(A+s_L^2)U_0-U_0(A+s_L^2)w_s is the derivative of
w_s' U_0-w_s U_0', and its value at L is U_0(L).
Consequently (7) has a solution IF AND ONLY IF M_s(f)=0, and that
solution is unique. This is stronger than only a necessary estimate.
The center is then regular: with f zero near L and the potential
even in L-z, the solution with U(L)=0 has the regular odd expansion.

For every s_L the moment is a nonzero bounded linear functional on
L2(0,L). To see it is nonzero even on compact smooth sources, choose
any nonnegative smooth bump chi supported where w_s is nonzero and
set f=chi w_s. Then M_s(f)=integral chi w_s^2>0.
Multiplying f by any nonzero small constant preserves the failure.
On that compact interval ker(M_s) is a codimension-one closed
hyperplane and has empty interior in L2 and in the usual finite
Sobolev initial-data norms. No bound on the absolute size of the
source forces this signed compatibility condition to vanish.

One can construct compatible sources without an approximation:
for any compact smooth q, f=(A+s_L^2)q has M_s(f)=0 and U=q.
Thus the statement is not that every source fails, but that the
stronger boundary restriction is not onto the unrestricted source space.
The functional M_s depends on frequency; the stationary signed source
condition is not a substitute for all of these dynamical conditions.

## 4. Initial-data consequence and its exact scope

Suppose a proposed realization of the FULL scalar equation admits an
exponentially bounded retarded evolution for all compact initial data
in a neighborhood of zero, with sufficient regularity that its Laplace
transform solves the spatial equation and preserves both zero traces.
The latter means the trace/falloff domain is closed in the evolution
norm; pointwise decay without such control would not justify a Laplace
boundary condition.

For u(0)=u_0 and u_t(0)=u_1, Laplace transformation gives

    (A+s_L^2)U=s_L u_0+u_1.                            (9)

Choose s_L larger than the hypothesized growth bound, set u_0=0,
and take u_1=epsilon chi w_s as above. These initial data vanish on
neighborhoods of both endpoints and can be arbitrarily small in every
fixed finite Sobolev norm. Yet M_s(u_1)>0 for epsilon>0, contradicting
(7)-(9). Therefore no such retarded evolution exists on an open set
of the full scalar initial-data space with the two-zero-trace rule.

The comparison metric lift (1) makes these scalar initial data into
compact conformal metric jets; the higher time jets are fixed by the
same KG equation. The local comparison tensor constraints are then
satisfied. This does not construct charged-exterior initial data.
Restricting the initial-data manifold to a separate perturbative branch
could exclude these jets, but that would be an additional dynamical
prescription needing physical justification. It is not proved by
static fast falloff or by projecting a stationary equation source.

This argument excludes an exponentially bounded global response with
the stated strong traces. It does not exclude local evolution before a
boundary is reached, nor prove failure of every possible finite-time
physical survival criterion. S077's finite-time nonlinear requirement
cannot be replaced by this comparison result.

## 5. Physical frequency and verification fixture

The oscillating phase is k log z. The physical orthonormal radial
derivative near this AdS boundary is sqrt(lambda)sin(z)partial_z.
Its limiting proper phase wavenumber is therefore

    k_proper=sqrt(lambda) k,
    G k_proper^2=60pi+17G lambda/12>180.                (10)

It exceeds 13 inverse Planck lengths when l_p^2=G. Sending support
toward the boundary does not lower this proper phase scale. This
observation does NOT supply a physical cutoff or justify order
reduction: no admitted ultraviolet prescription has been established.
It prevents interpreting coordinate dilation alone as a low-energy
validity argument.

For independent finite verification one can evaluate (6) explicitly.
Let p=1/2+ik and c=3/2. Then

    w_s(z)=-cos(z) sin(z)^p
       2F1((p+1+i s_L)/2,(p+1-i s_L)/2;3/2;cos(z)^2). (11)

The expression is real by uniqueness of its real differential equation
and center data, even though its representation uses complex powers.
It follows by substitution into the hypergeometric equation; the
general radial hypergeometric form appears in Ishibashi--Wald section
3.2. Here substitution directly verifies its continuation to imaginary
BF parameter. No positivity theorem from the real-parameter case is
continued into this case.

The focused checker verifies the boundary Wronskian, Green identity,
hypergeometric differential equation and normalization, moment signs,
and proper-scale coefficient. Independent Runge--Kutta integration of
the forced ODE with zero left data compares U_0(L) with quadrature of
(8) using (11), on two grids. Compatible compact sources are tested as
well. The moderate coefficient fixtures test the identity; they do not
simulate the enormous benchmark frequency or an actual exterior.
The continuum/domain result is the ordinary proof above, not formal
certification. No Lean run is selected.

## 6. Research decision

Do not identify the stationary fast remainder space with a completed
physical dynamical boundary prescription. Its exclusion of both scalar
amplitudes retains a signed range condition, and the dynamical analogue
cannot handle an open ball of the full retained scalar initial data.
The stationary S167 compatibility problem is preserved, not solved or
contradicted. The general below-BF theorem is credited to prior work;
the explicit frequency-dependent range test is the present diagnostic.

The neutral current decoupling and Ward recovery remain available for
the actual charged metric problem. Next derive its physical metric
boundary realization and control the variable-curvature retarded
coupling, or justify an admitted disturbance/ultraviolet prescription
that genuinely excludes the fast branch. The reference equations,
action coefficients and preparation data remain unchanged.
S168.3-S168.4 and all S167 tasks stay open.

Remaining to-do list: actual charged metric boundary and constrained
retarded transfer or justified physical exclusion; S167 signed
coefficient, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, full-sector stability, transport,
preparation, separate creation and full audit.
