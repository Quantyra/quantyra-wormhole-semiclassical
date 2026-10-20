# The TT mass functional and an adjoint residual certificate

2026-10-20 HST, E003/S166. Ordinary proof and contract before execution.
**Exact mass extraction and a certified approximation method. The actual
adjoint defect norm and global fast fourth/log inverse remain unbounded
numerically. S166.2 is not complete.**

Keep the same seed, action, physical pairing, clock, state, magnetic
flux and paired domain. This is a consequence of the actual
[weighted TT inverse](stationary-weighted-tt-inverse.md) and
[global weighted map](stationary-global-weighted-map.md), not an inverse
of an order-reduced equation.

Write u=B/r, u_h=B/r_h<1, W=u^-2, g=G lambda and

    dmu=gamma0 u^-4 du, gamma0=4pi beta B^3,
    Dbar=D/lambda^2, G=Dbar_crit^-1.

The symbol G here denotes the dimensionless critical TT inverse;
Newton's constant occurs only through g. In this note every source f
is the scalar TT source in the physical pairing, not a disk source.
The final conversion to the coupled disk problem is stated separately.

## 1. One actual adjoint solve

The constant sphere component 1 is smooth at the bolt and belongs to
the auxiliary TT space W_delta for every -1<delta<0. Its TT lift and
fixed outer return are the same global maps as for decaying inputs.
This auxiliary test field changes no physical boundary datum.

Define, using the actual operator,

    h=Dbar 1, v=G h, psi=1-v.                       (1)

These definitions are legitimate in the following domains. The normal
operator annihilates 1, since N_D(0)=0. The actual-minus-normal estimate
gains two powers, with finite logarithms; the separated compact-to-end
tail has degree four. Thus h belongs to the critical physical L2 space.
This assertion can also be obtained on an outer cutoff of 1: its normal
commutator has degree four, while the coefficient error has degree two.
The smooth compact part is regular at the bolt.

The constant is NOT put into Dom(Dbar_crit). Apply that inverse only
to h, which is in its target. Then v belongs to the actual critical
fourth/log domain, and distributionally

    Dbar psi=0, psi-1 in W_(3/2).                   (2)

If two such regular fields existed, their difference would lie in the
critical kernel. The accepted positive gap makes that kernel zero.
Hence (2) uniquely normalizes the actual zero-root adjoint field.
It includes its forced outer corrections and all global state terms.
Replacing psi by the leading constant 1 would generally be incorrect.

The critical gap in the scalar sphere pairing is

    Dbar >= .99/(8pi g), ||G||_L2->L2 <=8pi g/.99<26g. (3)

Since ||W^-1||_L2->L2=u_h^2<1 and
||W^-1 1||^2=gamma0 u_h, psi has the required reciprocal weighted norm.
No numerical full fourth/log inverse estimate is used here.

## 2. The boundary constant and exact mass

For any fast source f, the already established mass-extended solution
has the unique expansion

    c=Gf=m_c chi(u)u^3+c_fast, c_fast in W_(7/2).   (4)

The cutoff is one near infinity. Changing it does not change m_c.
Retain the normal scalar TT symbol

    N_D(s)=9lambda^2 N_T(s),
    N_T(s)=q[d0(q+2)H(s)-A_T], q=s(s-3),
    d0=1/(480pi^2),
    A_T=1/(48pi g)+23/(8640pi^2),
    H(3)=1-gamma_E+kappa_M, 39/5<kappa_M<54/5.

Define the dimensionless boundary constant

    K=-9gamma0 N_T'(3)
     =27gamma0[A_T-2d0(1-gamma_E+kappa_M)].        (5)

Its sign and normalization are those of the source/mass pairing in
section 4 of the weighted TT inverse. In particular

    gamma0/(6g)<K<gamma0/(5g).                    (6)

For the lower bound use 0<gamma_E<1, H(3)<12, 3<pi<22/7
and g<2000^-10; for the upper bound drop the negative H(3) term.
The leading coefficient is 9/(16pi g). This is the scalar-c mass
normalization. The disk mass coefficient is (3/4)m_c, not m_c.

The Green pairing gives the exact identity

    <1,f>-<h,c>=K m_c.                            (7)

Here is why no uncomputed interior term is missing. Apply formal
symmetry with an outer cutoff of the constant and use (4). The regular
bolt contributes zero. The cutoff commutator limit for the mass cubic
is the normal root pairing -9gamma0 N_T'(3). Terms from the fast
remainder tend to zero. The two-power actual-minus-normal gain makes
its boundary pairing vanish; separated degree-four tails also vanish.
Every interior contribution is already in h=Dbar 1. This uses the
finite-strip logarithmic Green/contour argument proved in the global
weighted map, not an assumption that the actual operator is local.
The underlying Mellin residue construction is also described in
[Melrose, sections 5.17 and 6.2](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/melrose.pdf).

Both integrals in (7) converge: Wf is L2 and W^-1 is L2 as a
function, while h and c are critical L2. Critical symmetry gives
<Gh,f>=<h,Gf>=<h,c>. Therefore

    m_c(f)=<psi,f>/K.                             (8)

As a sign check, the leading local normal operator is
-9A_T[(u partial_u)^2-3u partial_u].
Its lower-end outward flux on m_c u^3 is +27gamma0 A_T m_c,
in agreement with (5). The finite quantum part of (5) is retained.

## 3. Exact operator norm and a residual error bound

Equip the full fast TT source space with ||Wf||_L2. Its mass map
has the exact norm

    ||m_c||=||W^-1 psi||/K.                       (9)

Cauchy--Schwarz proves the upper bound. Equality is attained by a
multiple of f_*=W^-2 psi, since Wf_*=W^-1 psi is L2. This source
belongs to the fast TT target. The statement concerns the whole TT
source space; it does not assert that f_* is the preparation source
or independently satisfies the coupled scalar compatibility condition.

Now let v_a be any real approximation in Dom(Dbar_crit) whose actual
strong residual is known:

    psi_a=1-v_a, r_a=Dbar v_a-h.

The domain includes regular bolt behavior and the original outer
realization. Then

    psi-psi_a=G r_a,
    ||W^-1(psi-psi_a)|| <=26g ||r_a||.             (10)

Thus a certified L2 residual for ONE critical adjoint solve gives

    |m_c(f)-<psi_a,f>/K|
       <=(26g/K)||r_a|| ||Wf||,                   (11)

and, putting A_a=||W^-1 psi_a|| and E_a=26g||r_a||,

    max(0,A_a-E_a)/K <=||m_c||<=(A_a+E_a)/K.       (12)

The exact identity also gives the sharper source-specific error
<r_a,c>/K if an independent bound for c is available.
Equations (10)-(12) require the full actual operator residual,
including its nonlocal tails. A sampled differential residual, a
truncated angular determinant or a normal-model residual is insufficient.
Any quadrature or finite-interval error in A_a or the pairings must
also be enclosed before using a numerical certificate.

For the uncomputed zero trial v_a=0, define the finite number

    Xi=||Dbar 1||/sqrt(gamma0).

It follows that

    ||m_c|| <(6g/sqrt(gamma0))[sqrt(u_h)+26g Xi]
             <(6g/sqrt(gamma0))[1+26g Xi].         (13)

This is a reduction to a single explicit-probe norm, NOT a numerical
upper bound: Xi has not been enclosed in this note. Positivity and
normal data alone do not supply its value. A nonzero trial with a
certified residual can instead use (12).

## 4. Consequence for the coupled problem

For the actual compatible physical solution, use the
[strong effective source](strong-effective-tt-source.md):

    f=fbar_c-Bbar* tau, ||Wf||<6.6e7 Fsrc
                                  <=8.184e9 lambda^-2 Fdisk.

Hence, conditional only on an eventual enclosure for Xi, (13) gives

    |m_c| <=(3.96e8 g/sqrt(gamma0))(1+26g Xi) Fsrc
          <=(4.9104e10 g lambda^-2/sqrt(gamma0))
                                      (1+26g Xi) Fdisk. (14)

For an approximate adjoint, use (11) with the same source conversion.
The physical disk cubic has coefficient (3/4)m_c. No second mass,
clock adjustment or preparation datum has been introduced.

This separates the mass estimate from the global fast fourth/log
inverse estimate. The next mass calculation can bound Dbar 1 directly,
or construct v_a and certify its actual L2 residual. The fast remainder
still needs quantitative compact, bolt, coefficient, commutator and
nonlocal-tail estimates. Higher reconstruction, full residual/inverse
composition, the nonlinear family and signed preparation condition
remain open.

## Frozen verification contract

Freeze this ordinary proof and checker before execution. Revalidate
the exact receipts and current input hashes of the global weighted
map, weighted TT inverse, quantitative normal inverse and strong
effective source.

Check symbol differentiation at the mass root, the boundary sign,
physical measure conversion, reciprocal norm, gap reserve, K interval,
residual identity and coupled-source constants by exact arithmetic.

For an independent manufactured diagnostic, use on 0<u<1

    L_b=-u^2 partial_u^2+2u partial_u+V_b,
    psi_b=1-b u^2+(b/2)u^4,
    V_b=2b(u^2+u^4)/psi_b, b=0,1/10,1/2,9/10.

Use measure u^-4 du, outer critical L2 and Neumann data at u=1.
This is a toy Sturm--Liouville endpoint, not the physical bolt.
The local gap is at least 9/4: conjugation to x=-log u gives the
nonnegative derivative energy, 9/4 times the norm and the positive
3/2 endpoint term. V_b>=0 and psi_b>=1-b/2>0.
It has h_b=V_b, v_b=1-psi_b, L_b v_b=h_b and boundary constant 3.

For n=3,4,5,7,9 test c_n=u^n-n u^(n+1)/(n+1), whose Neumann
data vanish. Verify exactly that both mass pairings are 3 for n=3
and zero otherwise. Test v_a=t v_b at t=0,1/2,1,3/2.
At 80 decimal digits integrate the nonsingular squared norm
(V_b/u^2)^2, the reciprocal adjoint norms and all source pairings.
Use identity tolerance 1e-65; inequalities need their correct sign
without a tolerance. Verify (10)-(12) with the model gap constant 4/9.
These fixtures test pairing, normalization and error accounting;
they do not evaluate the actual Dbar, Xi, psi or preparation source.

Continuum mapping, Green limits, inverse existence and density remain
ordinary proofs. Preserve every run and failed scientific source.
No actual numerical mass upper bound, complete S166 inverse,
nonlinear existence or physical stability is certified by this check.

Remaining to-do list: enclose the actual adjoint defect or solve it
with a certified residual; quantify the global fast TT remainder and
full inverse/residual composition; nonlinear family and signed
compatibility; exterior/state, charge and return; opening, stability,
transport, preparation, creation and full-result audit.
