# Joint charge and spin correlation from the fermion reservoir

S072, 2026-10-13. **The physical thermal boundary correlator of a
specified global electric/angular Casimir model has an exact finite-Fock
formula and a controlled optical-volume limit.** Its compact factor,
including the charge-spin selection rule and dimension-squared weights,
emerges from the actual fermion trace. It is not an added independent
rotor. This replaces a missing state assumption for this comparison
model and reaches the point field, rather than the normalized canonical
mode average of [fermionic-compact-transfer.md](fermionic-compact-transfer.md).

The Hamiltonian below is not the full radial angular Hamiltonian or a
matched black-hole state. In particular its electric term is a *global*
charge square, whereas [radial-electric-correlator.md](radial-electric-correlator.md)
retains the enclosed-charge profile. Their two answers must not simply
be multiplied. Physical horizon gluing, angular radial dynamics, massive
interactions and the driven stress remain unresolved.

## Finite Hamiltonian and the operator to be matched

Take radial orbitals r with dimensionless energies u_r>=0. Each contains
q=2j+1 particle modes of electric charge +1 and spin j, and the conjugate
antiparticle multiplet. The vacuum is a neutral spin singlet. On their
actual finite fermion Fock space set

```
Hcal = beta H = H0 + x Q^2 + a J^2,
H0 = sum_r u_r (N_r + Nbar_r),       x,a>0, q>=2.       (1)
```

H0 commutes with Q and all J_a. No extra charge or spin Hilbert space is
introduced. Individual multiplet occupation numbers commute with (1),
but the spin components within them need not be separately conserved.
For 0<t<1 define the spin-averaged particle insertion

```
b_r(t) = Tr[exp(-(1-t)Hcal) c_r^dagger exp(-t Hcal)c_r]/(q Z),
b_r,0(t) = exp(-(1-t)u_r)/(1+exp(-u_r)),
r_r(t) = b_r(t)/b_r,0(t),                              (2)
```

where the first trace includes the sum over the q spin entries. At
t=1/2 this is the positive TFD midpoint matrix element in the previous
normalization. Away from the midpoint r_r(t) need not equal r_r(1-t)
at finite volume. The limit below does have that symmetry.

The one-particle representation is R(g)=exp(i alpha)D^j(g_s). Its faithful
group is G=(U(1)xSU(2))/Z2 for even q and G=U(1)xSO(3) for odd q.
Write d_l=2l+1 and E_Ql=x Q^2+a l(l+1). With normalized Haar measure,
the positive heat kernel of this group is

```
p_s(g) = sum_(Q,l allowed) d_l exp(-s E_Ql) exp(i Q alpha) chi_l(g_s),
allowed: 2l=(q-1)Q mod 2,
P_s = p_s(e) = sum_allowed d_l^2 exp(-s E_Ql).          (3)
```

The kernel is an integral representation of exp[-s(xQ^2+aJ^2)] on
*any* representation of G, including the Fock representation in (1).
Using (3) in a trace therefore does not assign rotor multiplicities to
the fermions. [Iliesiu, section 2.2](https://arxiv.org/html/1909.05253v2)
gives the separate group-particle heat-kernel construction. Here the
Fock determinant, retained below, determines the finite-state weights.

## Exact determinant and insertion formula

For a tagged orbital set z=exp(-u_r), and for every orbital s put

```
A_z(g) = (1+z R(g))/(1+z),
f_z(g) = |det A_z(g)|^2,
Phi_rest(g) = product_(s != r) f_(exp(-u_s))(g),
B_z(g) = det(A_z(g))^* adj(A_z(g)),
F_z(g,h) = tr[B_z(g) R(h)]/q.                           (4)
```

B_z=f_z A_z^{-1} where the inverse exists. The cofactor definition is
regular when z=1 and an eigenvalue of R(g) is -1, and is the definition
used in the checker. The exact result is

```
C_z,t(g) = integral_G p_(1-t)(h) p_t(h^{-1}g) F_z(g,h) dh,
D_z(g) = p_1(g) f_z(g),
r_r(t) = [integral Phi_rest(g) C_z,t(g) dg]
         /[integral Phi_rest(g) D_z(g) dg].             (5)
```

Derivation. The normalized free Fock character is product_s f_zs(g):
the particle and conjugate antiparticle determinants pair into a
nonnegative modulus square. This is the finite trace identity of
[Klich, section 2, equations (8)--(12)](https://arxiv.org/pdf/cond-mat/0209642).
Inserting the two heat kernels gives two group arguments h and h^{-1}g.
For exterior-algebra operators Gamma(A), Gamma(B), the CAR trace with
one insertion is

```
sum_m Tr[Gamma(A)c_m^dagger Gamma(B)c_m]
   = det(1+AB) tr[(1+AB)^{-1} A].                       (6)
```

Move Gamma(A) through c^dagger and use
Tr[Gamma(AB)c_n c_m^dagger]=det(1+AB)[(1+AB)^{-1}]_nm.
For the tagged particle block A=z^(1-t)R(h), AB=zR(g).
Division by (2), including its free thermal denominator, gives (4)--(5).
All other blocks give their determinants. Equation (6) also extends by
continuity through a singular inverse. This derivation keeps the full
fermion representation, not merely its necessary center selection rule.

The trace is real and positive even though an individual integrand in
(5) can be complex. Treating that numerator as a positive classical
orientation average at finite volume would discard its resolvent.

## Uniform bounds on the insertion, including zero radial energy

Let delta(g)=||R(g)-1||_F^2 and eta(g)=||R(g)-1||_op. If exp(i theta_m)
are the eigenvalues of R, put p=z/(1+z), v=p(1-p). Then

```
f_z(g) = product_m [1-v |exp(i theta_m)-1|^2],
0<=f_z<=1,       1-f_z<=delta/4,
||B_z||_op<=1,
||B_z-1||_op <= ((2q-1)/2) eta.                         (7)
```

To see the last two statements, diagonalize R. Each eigenvalue of B_z
is a product of q conjugate eigenvalues of A_z and q-1 unconjugated
ones. Each factor has modulus at most one and differs from one by at
most p eta<=eta/2. Telescoping the product proves (7). Thus |F_z|<=1,
uniformly in every energy including u_r=0. In particular

```
|C_z,t(g)|<=p_1(g)<=P_1,        0<=D_z(g)<=P_1,
C_z,t(e)=P_1 C_G(t),           D_z(e)=P_1,
C_G(t) = integral p_(1-t)(h)p_t(h^{-1}) tr R(h) dh/(q P_1). (8)
```

These uniform estimates permit the ultraviolet point-field limit; no
divergent equal-point CAR norm appears in the error bound.

## Explicit finite-reservoir certificate

For the *remaining multiplets*, not their individual spin entries, set

```
V = sum_(s != r) exp(-u_s)/(1+exp(-u_s))^2,
S_J = q(q^2-1)/12,
c_q = 1/[48 e pi^4 sqrt(q) S_J^(3/2)].                  (9)
```

For V>=1, the probability measure dmu=Phi_rest dg/integral Phi_rest
satisfies, for any epsilon>0,

```
mu{delta>epsilon^2} <= T(V,epsilon)
  := min(1, c_q^{-1} V^2 exp(-V epsilon^2)).             (10)
```

Proof. The product in (7) gives Phi_rest<=exp(-V delta). On the cover
U(1)xSU(2), choose |alpha|<=1/(2 sqrt(q V)) and rotation angle
0<=theta<=1/(2 sqrt(S_J V)). Since sum m=0 and sum m^2=S_J,
delta<=q alpha^2+S_J theta^2<=1/(2V). Each factor then obeys
log(1-v_s |exp(i theta_m)-1|^2)>=-2v_s |exp(i theta_m)-1|^2,
so Phi_rest>=exp(-1). The normalized class measure on SU(2) is
sin^2(theta/2)dtheta/pi, bounded below by theta^2 dtheta/pi^3 on this
box. Its product with d alpha/(2pi) integrates to at least
1/[48 pi^4 sqrt(q) S_J^(3/2) V^2]. Pushing the cover measure onto G
does not change integrals of these invariant functions. This proves
integral Phi_rest>=c_q/V^2, and hence (10).

All constants needed to control the heat kernels are explicit convergent
representation sums. Define

```
L_s = 2 sum_allowed d_l^2 (|Q|+l/j) exp(-s E_Ql).        (11)
```

If delta(g)<=epsilon^2 with 0<epsilon<1/4, then on the faithful quotient
there is a lift near the identity with |alpha|<=2epsilon and
|theta|<=2epsilon/j. Indeed all eigenphases of R are within
2 asin(epsilon/2) of zero. Adjacent spin weights differ by theta;
these phase differences have no 2pi wrap in this neighborhood. Their
average bounds alpha, and their total spread bounds j theta. A lift
near the other central element represents the same point of G.
Termwise character estimates in (3) therefore give, uniformly in h,

```
|p_s(h^{-1}g)-p_s(h^{-1})|<=epsilon L_s.                (12)
```

Use (7) and (12) on the near-identity region and (8) on its complement.
Writing T=T(V,epsilon), obtain

```
U = epsilon [L_t/P_1+(2q-1)/2]+2T,
W = epsilon L_1/P_1+epsilon^2/4+T,
|r_r(t)-C_G(t)| <= (U+W)/(1-W),       if W<1.            (13)
```

For the numerator, shift its second heat kernel first; the other kernel
integrates to one. Then bound B_z-1 at fixed g using (7), with the
unshifted kernel product integrating to P_1. The resulting normalized
error is at most U. The corresponding denominator differs from P_1
by at most W. Finally 0<=C_G(t)<=1, either by its thermal spectral
formula below and Holder's inequality or by the rotor trace identity,
so division gives (13). This is uniform in the tagged energy and in
the number and energies of all other orbitals once V is specified.

For example epsilon=sqrt[(log(c_q^{-1})+4 log V)/V] gives T<=V^{-2}
and an error tending to zero at fixed q,x,a,t. Equation (13), not an
unspecified exchange of q and volume limits, is the finite-volume
statement. It is conservative and need not be small at modest V.

## Emergent compact factor and physical point field

Character multiplication evaluates (8) as

```
C_G(t) = [1/(q Z_G)] sum_(Q,l,l' allowed) d_l d_l' N_(l,j)^l'
              exp[-(1-t)E_Ql-t E_(Q+1,l')],
Z_G = sum_allowed d_l^2 exp(-E_Ql)=P_1.                 (14)
```

N is the SU(2) triangle/parity multiplicity. This is precisely the
group-particle thermal factor, but here it is *derived as a reservoir
limit of (5)*. It is not the finite Fock partition function. In
particular forbidden finite-occupation branches were never inserted.
Equations (8) or (14) also give C_G(t)=C_G(1-t).

For the antiperiodic optical circle, E_n=2pi(n+1/2)/ell and each boundary
spin entry has spatial normalization 1/sqrt(ell). There is one paired
spin multiplet per positive radial energy. Thus V_total~ell/(4pi beta),
and removing a tagged multiplet changes V by at most 1/4. Formula (13)
therefore applies uniformly to *all* radial energies as ell tends to
infinity at fixed q,x,a. Different radial contributions have no mixed
two-point functions because their separate particle/antiparticle
occupation numbers are conserved. Rotation invariance makes the spin
two-point matrix a scalar. Particle annihilation and antiparticle
creation contribute the t and 1-t thermal weights respectively.

After spatial smearing and removal of the ultraviolet regulator,

```
G_point(beta t) = C_G(t)/(2 beta sin(pi t)),   0<t<1.   (15)
```

The free integral is
integral_0^infinity du [exp(-t u)+exp(-(1-t)u)]/(1+exp(-u))
=pi/sin(pi t). At fixed ell the ultraviolet determinant product
converges because sum_n exp(-beta E_n)<infinity. The bound
|C_z,t|<=P_1, Z>=1 and finite Z_free bound every r_r(t) by
P_1 Z_free at that ell, so the decaying free thermal weights justify
the regulated sum limit. For increasing ell, (13) supplies a bound
uniform in r once V_total-1/4 is large. This permits the point limit
and proves (15), including either order of a fixed physical energy
cutoff and the optical-volume limit. Nothing here proves the analogous
limit for a fluctuating horizon or nonlocal radial spin Hamiltonian.

At the midpoint the already proved compact affinity inequality from
[quantum-junction-state.md](quantum-junction-state.md) now yields for
this physical point-field comparison

```
0 <= 1-2 beta G_point(beta/2) = 1-C_G(1/2)
  <= [x+a j(j+1)]/4.                                  (16)
```

If the earlier static candidate coefficients are the coefficients of
(1), they obey x+a j(j+1)<1/(3k), k>=2000. The limiting midpoint loss
in this model is then below 1/(12k), or 4.167e-5 at k=2000. At a
finite optical cut one must add the explicit error (13). We have not
shown that a physical horizon state supplies those coefficients or a
sufficient reservoir, and do not assign this number to the actual
wormhole source. The existing finite-cut *radial electric* bound has
a different scope and does not establish the angular premise here.

## Verification and next physical calculation

The checker compares (5) against direct diagonalization of (1), including
q=2,3,4, paired radial energies, a zero-energy orbital, and non-midpoint
insertions. The integral uses heat kernels on the full cover; its Fock
determinant enforces the center quotient. The independent limit uses
only the allowed-representation and fusion sum (14). Increasing radial
reservoirs test that limit, with a separate quadrature refinement.
The cofactor estimates and explicit certificate constants are also
checked. These are supplementary numerical tests of the ordinary proof,
not a numerical extrapolation to q=k^4 or a certification of geometry.

All six finite Fock/group comparisons agree to less than 2.66e-14.
Nine reservoir cases approach the independently summed compact limit;
the three quadrature refinements differ by less than 3e-10. The initial
trial V=1e8 did not achieve a requested two-percent analytic error budget
(the worst bound was about 0.505). That failed test and its source are
retained. Increasing V to 1e12 makes all three trial bounds smaller than
0.0061, and V=1e16 reduces them further. The original V=1e8 rows remain
in the receipt explicitly marked as missing the two-percent target.
These are numerical evaluations of the convergent sums in (11), not
interval-arithmetic certificates or a physical finite-cut assertion.
The ordinary bound (13) and the Hamiltonian did not change.

The next physical step is to derive the angular radial Gauss Hamiltonian
and its dressed insertion in the *same* gravity/matter gluing used for
the electric sector. Test whether its nonzero-mode contribution can be
controlled relative to (1) with an error below the capture margin; the
current/coset algebra identifies its nonlinear content but does not
justify replacing it by J^2. Keep the bare RN Hessian as a matching
condition, not an independent partition-function multiplier. Only then
use the resulting initial correlators in a driven covariance/stress
calculation with common subtraction and evolving geometry.

Remaining to-do list: physical gravity/matter and radial angular matching;
driven state, complete stress, massive effects and geometry; preparation,
quantum fluctuations, sufficient stability and finite-energy transport.
