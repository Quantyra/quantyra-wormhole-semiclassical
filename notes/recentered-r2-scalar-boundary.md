# Recentered scalar response and the paired Dirac boundary test

2026-10-22, E003/S169. Ordinary scalar-row result and finite-cut boundary
obstruction. The coupled semiclassical solution and S169.2 remain open.

Continue [the charged trace block](charged-r2-trace-block.md) at the SAME
additional c=10^22. The preceding goal turn was progress: technical
0997995 and planning 868453d supplied the charged comparison and exposed
the nondecaying source and unresolved boundary realization. This note
executes both checks rather than assigning zero to either remainder.

The exact quantum-vacuum recentering makes the actual scalar trace source
an admissible energy-dual source. Its unique linear conformal response
has a dimensionless radial energy bound <1.82e-29. Its forced leading
tail is r^-4 log r, with explicitly bounded coefficient. The paired Dirac
normal heat terms do NOT vanish on the finite-cut metric/curvature
Dirichlet conditions alone. They do vanish in the actual AdS limit for
this response class. This is a scalar-sector result, not a full boundary
domain, metric admission bound, quantum state, or nonlinear wormhole.

## 1. Recenter the geometry without changing the theory

Put g0=G lambda0=k^-10, k=2000, and let

    g+11g^2/(360pi)=g0, lambda=g/G,
    epsilon=G/(60pi), D=96piG c-epsilon,
    C_v=1+11epsilon lambda/3,
    F_lambda=1-384piG c lambda,
    Z=piG q^2/(g_R^2 F_lambda), q=k^4, g_R=k^-4.       (1)

Use the new reference f=1-2M/r+Z/r^2+lambda r^2, with the original
fixed flux, coupling and prescribed positive temperature. Its extremal
radius B obeys Z=B^2+3lambda B^4. The same monotonicity proof as before
gives one r_+>B for each prescribed positive temperature. M is then
(r_++Z/r_++lambda r_+^3)/2. These are new reference data, not the old
finite geometry or a recomputed quantum product.

The local tensor defect is the constant -3(lambda0-lambda)g_ab; it is
exactly the Dirac vacuum stress 8piG T_vac,ab at infinity. The rest of
the quantum stress is not determined by this observation. The reference
remains off shell for the full Dirac theory. No coupling, mass condition
or temperature has been adjusted to kill a signed compatibility moment.

On this geometry R=-12lambda, J has eigenvalues (-t,-t,t,t), t=Z/r^4,
and C^2=48W^2, W=M/r^3-Z/r^4. Equation (1) cancels the constant trace
residual EXACTLY. The remaining bulk trace is

    T0=epsilon H,
    H=28M^2/r^6-56MZ/r^7+(106/3)Z^2/r^8-10q^2/r^4.   (2)

This is the actual trace of the full angular Dirac anomaly on the
specified reference. The conformally invariant, state-dependent shape
stress remains uncomputed. A local trace identity does not construct
the physical state or verify its gap.

## 2. Exact scalar operator, source norm and actual linear response

The derivative of the fixed-volume row exp(4s)T[exp(2s)g] is now

    A=6(-Delta+4lambda)(-D Delta+C_v)
                                   -(22epsilon/3)div(J grad).        (3)

The corrected zeroth-order coefficient is 24lambda C_v. This differs
from the preceding bare-lambda0 comparison. The density contact still
cancels all Weyl and gauge zeroth-order potentials; T0 remains a source.
The gradient coefficient K is
(6+24D lambda+22epsilon lambda)g+(22epsilon/3)J.

Rational bounds using both root endpoints give

    B^-2/lambda<1.223e-10,
    K>=(6+7.06e-8)g, C_v>=1, D>0.                    (4)

Therefore the complete-end energy form is coercive in the explicit norm

    ||s||_E^2=int[6D(Delta s)^2+6|grad s|^2+24lambda s^2]dvol.

As before the complete-end space is the E-completion of smooth compact
functions regular at the bolt. The bounded positive form gives a unique
linear solution s_tr to A s_tr=-T0. Rotational and Euclidean time
symmetry of the form, source and domain imply that this solution is
radial. This is the actual response of this scalar comparison row, not
an assumed value for an unknown inverse constant or for the full metric.

Remove the common factor 4pi times the Euclidean time period from both
radial integrals. Write E_rad for this radial norm. Since T0 now decays,

    ||s_tr||_E,rad <= ||T0||_L2,rad/sqrt(24lambda),
    ||T0||_L2,rad^2=epsilon^2 sum_(p,q in {4,6,7,8})
                         a_p a_q r_+^(3-p-q)/(p+q-3),                (5)
    (a4,a6,a7,a8)=(-10q^2,28M^2,-56MZ,(106/3)Z^2).

The sum in (5) retains signed cross terms and is an EXACT integral,
not a quadrature certificate. For a bound uniform in prescribed positive
temperature put beta=1.223e-10 and j=G lambda. The geometric inequalities
M/r_+^3<=(2+beta)lambda and Z/r_+^4<=(3+beta)lambda imply

    U=q^2/(lambda^2 r_+^4)
                    <=g_R^2(3+beta)/(pi j),
    N_tr:=||s_tr||_E,rad/sqrt(lambda r_+^3)
       <= (epsilon lambda)/sqrt(24) [10U/sqrt(5)
           +28(2+beta)^2/3
           +56(2+beta)(3+beta)/sqrt(11)
           +(106/3)(3+beta)^2/sqrt(13)]
       <1.82e-29.                                       (6)

The factor sqrt(lambda r_+^3) is displayed to connect the energy norm
with a dimensionless metric amplitude scale. It does not make the
underlying physical bound smaller: (5) is the unrescaled bound, and no
admission conclusion is drawn from choosing units. N_tr is not a cold
supremum, boundary-jet norm or bound on the actual coupled correction.
Such embeddings, especially near a very cold bolt, still need proof.

## 3. The forced tail and its allowed outer boundary data

Set y=1/r. The radial fourth-order equation from (3) has analytic
regular-singular coefficients at y=0. Its indicial polynomial is

    P(p)=6lambda[4-p(p-3)][C_v-D lambda p(p-3)].          (7)

Its four roots are -1, 4, and
Delta_+/-=3/2 +/- sqrt(9/4+C_v/(D lambda)); Delta_+>18000.
No zero or repeated root is hidden: C_v-4D lambda>0.99999998.
The coefficient perturbations from 1, M/r, Z/r^2 in f, and from J,
are analytic in y and start at positive powers. Frobenius reduction
therefore gives powers with finite logarithmic multiplicity; possible
integer resonances do not change the leading admissible power. The
energy condition excludes the -1 and Delta_- branches. The source's
first term is +10epsilon q^2 y^4 on the right side of A s_tr=-T0.
Because -P'(4)=30lambda(C_v-4D lambda), the asymptotics are

    s_tr(r)=a r^-4 log(r/r_+)+b r^-4
                      +O(r^-6 log^p r)+O(r^-Delta_+ log^p r),
    a=epsilon q^2/[3lambda(C_v-4D lambda)],
    0<a/r_+^4<6.60e-30.                                (8)

Only a is fixed locally. b and the fast amplitude are determined by the
global bolt-regular solution, not set to zero. The differentiated bounds
through any fixed finite physical radial order follow from the same
regular-singular system. In particular the first three derivatives in
a smooth canonical coordinate Y proportional to 1/r tend to zero.

The r^-4 log r term generally is not C4 in Y. Thus old smooth-collar
C5 estimates for the full shape determinant CANNOT simply be imported.
A conormal/polyhomogeneous version or an explicitly sufficient weaker
regularity theorem is required for the coupled problem. The boundary
calculation below is done at smooth finite cuts and uses only the finite
jet orders occurring in a4 before taking the limit.

## 4. Execute the actual paired Dirac normal-jet test

Use the retained [paired domain](curved-fermion-boundary.md), with fixed
constant boundary angle and identical static collars. Let n point inward,
L_ab=-<nabla_a n,e_b>, L=tr L_ab, L2=L_ab L^ab. The squared Dirac
operators have equal-rank mixed projectors Pi_+,Pi_-, chi=Pi_+-Pi_-,
S=-L Pi_+/2. Let N be the bundle rank: N=8 for the folded pair, or
N=4 per identical exterior. Both partner-square coefficients agree;
their determinant average adds no extra factor of two.

Specialize the actual mixed heat coefficient
[Vassilevich (5.33)](https://arxiv.org/html/hep-th/0306138v3#S5.E33).
For a scalar smearing omega vanishing at the boundary, its COMPLETE
normal-derivative contribution is

    a4_boundary(omega)=(4pi)^-2 int_boundary sqrt(h) N
            [(2L^2-3L2)omega_;n/140-L omega_;nn/60].    (9)

The three-derivative smearing term has coefficient tr chi=0. The
curvature and gauge normal-smearing terms contain tr chi or tr(chi E),
both zero in this matched minimal domain. The remaining trace identities
are tr S=-NL/4, tr S^2=NL^2/8 and
tr(chi_:a chi_:a)=N L2. Substitution of ALL normal rows, including 72SL,
240S^2 and -18chi_:a chi_:a, gives (9). This is not the previously
calculated angle-derivative contribution alone. Constant-angle unitary
conjugation proves its angle independence locally; no global spectral
equivalence is asserted. Unequal collar perturbations are not covered.

For the proposed finite-cut local metric/curvature Dirichlet tangents,
omega=0 and delta R=-6Delta omega=0. Since
Delta omega=omega_;nn-L omega_;n there,

    a4_boundary(omega)=-(4pi)^-2 int_boundary sqrt(h)
                   N(L^2+9L2)omega_;n/420.             (10)

At a finite RN cut L=f'/(2sqrt(f))+2sqrt(f)/r and
L2=f'^2/(4f)+2f/r^2. Equation (10) is nonzero for arbitrary allowed
normal first derivative. With only the local EH+cR^2 Dirichlet completion,
the paired determinant therefore leaves a boundary term. Intrinsic
counterterms fixed by the induced metric cannot cancel it on these
tangents. An additional prescribed boundary completion would need its
own derivation; none is introduced here. This is a failure of that
finite-cut physical-domain inference, not of the bulk inverse theorem
or of all AdS boundary realizations.

## 5. Actual AdS limit for the scalar response class

For (8), omega and all fixed-order physical normal derivatives are
O(r^-4 log^p r). At large r the physical RN curvature, its required
normal derivatives, L and L_ab are bounded. The boundary volume is
O(r^3). Every term in the full a4 boundary coefficient, including the
omega row, therefore tends to zero as O(r^-1 log^p r). For two such
conformal directions the differentiated boundary variation is
O(r^-5 log^p r). The local EH+cR^2 boundary variation has the same
vanishing linear power, with its fixed finite c retained. Usual local
intrinsic AdS counterterms in this retained sector also vanish on these
directions. No cut-dependent finite term is selected to force the limit.

These estimates are a direct smooth-cut limit. They do not apply the
smooth-boundary heat theorem to a C3 compactification or exchange an
unproved full heat-kernel remainder limit. The massless conformal
variation is the anomaly, with the paired state/domain transported as
specified; its independent shape dependence is still present. The
physical boundary scalar value is unchanged because delta R decays.
The r^-4 log r conformal correction also preserves the outer leading
metric, clock and mass coefficient. Its bolt value can still change
physical radius and other observables; (6) alone does not bound those.

Thus the unique linear scalar response has no leftover local/Dirac
conformal boundary variation at the actual AdS end. This supplies a
necessary boundary check for this response sector. It is not a theorem
that every E vector, arbitrary shape perturbation, source drive, or
finite-amplitude rough determinant family has the full physical domain.
The remaining shape boundary row and its mixed form must be derived.

## 6. Research decision and verification

Continue the selected theory using the exact-vacuum-centered reference
and the actual AdS limit. Do not promote finite-cut Navier conditions
to the physical paired boundary. Use the explicit source and scalar
response bound (5)-(8) when constructing the coupled solution. Next
derive the shape/mixed boundary form on a class allowing the forced
r^-4 log r tail, and bound the complete shape response. The previous
targets ||E||<=1/4 and ||L^-1 C|| ||B||<=3/8 still concern the unproved
transfer/coupled problem; this calculation assigns them no values.

The ordinary argument precedes verification. The checker derives the
normal heat traces with the actual rank-eight Clifford projectors for
both partner domains, checks the exact finite-cut obstruction, and
independently contracts the recentered RN curvature and trace. Exact
integration gives (5); rational bounds certify (4),(6),(8). At 60/90
digits direct radial source integration agrees with the signed formula,
and unexpanded charged boundary geometry verifies the finite-cut and
tail scaling. Indicial substitution checks the resonant coefficient.
No quantum state is sampled and no full coupled solution is certified.
Freeze all scientific sources and preserve every run, source hash and log.

Remaining to-do list: shape/mixed physical boundary form with the forced
tail and complete coupled inverse; finite geometry/state, full residual,
admission and nonlinear map; controlled solution and exact compatibility;
matching, validity, all-sector stability, charge/return, opening,
transport, preparation, separate creation and full-result audit.
