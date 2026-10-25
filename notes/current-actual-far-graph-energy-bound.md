# Actual distant-input graph return from the coupled solution's energy

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress (complete weak residual and linear mass bounded).
The full stable semiclassical wormhole objective remains incomplete.
The closed feasibility audit and original physical targets are unchanged.

OUTCOME TO TEST: the complete distant-input contribution to the ORIGINAL
localized graph is <1.24e-20 for the actual full linear solution. Its
original shape-row finite part is retained and <1.27e-34. These close
this particular unknown input return, not the complete graph, kappa_m,
independent phase derivatives, nonlinear admission or a wormhole theorem.

## 1. Attempt the complete graph identity before assigning estimates

Keep the actual two-row H, A=z^4 H=A0+R, source F, paired mass/log
domain, original common regulator and all fixed finite contacts. The
known full linear solution has E(u)<1.1e-10. Let

    u=m Phi_m+w,
    V=chi[-F-A(m Phi_m)-R w]+[A0,chi]w,
    G=(D-4)V_b-(D-6)V_nu/3, D=z partial_z.

The mass lift is the exact reference lift with its original .1/.2
transition; it vanishes before z=.2 and retains the regular-bolt
matching. Its local cancellation is not a fixed-temperature zero mode.
The ORIGINAL graph cutoff is the beta(11,11) logarithmic cutoff, one
through a=.05 and zero from b=.1, exactly as in
current-global-mass-conversion.md. No basis or norm changes here.

Split the unknown input at c=.2. On its distant part w_far=1_(z>=c)u,
chi w_far=0. Consequently the COMPLETE unknown distant contribution is

    -chi R w_far+[A0,chi]w_far=-chi A w_far.          (1)

The normal part cancels exactly. Bounding the two terms independently
would lose this cancellation. Every local differential/contact block,
including both local gravitational rows, Weyl and fixed-domain LLL,
has zero cross-support action. Thus the only surviving distant block
is the actual massive canonical shape Hessian Hm; the scalar row is zero.
Set f(z)=Hm[1_(z'>=c)T](z), where T=u_b is the ACTUAL coupled solution's
shape, not a freely selected source or independent phase input. Then

    V_far=(0,-chi z^4 f),
    g_far=z^-sigma G_far
          =-chi z^(5-sigma)f'-(Dchi)z^(4-sigma)f, sigma=4.005. (2)

The second term is essential. The constant coefficient f(0) disappears
from the differentiated row only where chi=1; it remains in the original
shape equation and its transition. We bound it separately below.

Equation (1) holds with matched regulated A,A0,R first. The input split
is used only in a separated integral, not as a discontinuous physical
metric or a new boundary condition. Smooth approximations to the split
have the same separated bounds and limit. The remaining near-diagonal
contacts and local rows stay in the UNRESOLVED complementary block.

## 2. Replace the previously unknown supremum by a proved weighted return

The fixed-seed output comparison gives ||T/z^2||L2(dz)<17.6 E.
The retained TT gradient estimate, with the same physical normalization,
gives ||sqrt(Q) T'||L2(dz)<130 E: its squared constant is
1e5 beta/(.99 F_l)<130^2. This applies to the shape component of any
physical test under the fixed York return; the conformal scalar changes
nu, not T. No energy-to-H3 or neighborhood embedding is asserted.

The actual bolt polynomial, already enclosed in the complete mass
residual calculation, gives, with delta=z_h-z and gamma>0,

    2gamma delta+delta^2/8 <= Q <=2gamma delta+1.001 delta^2,
    z_h>.99, gamma<1e-30, Q<.2.                     (3)

There is an anchor z_* in [.2,.3] with |T(z_*)|<=17.6 E/sqrt(.1).
The fundamental theorem and weighted Cauchy--Schwarz yield

    |T(.2)|/E <=17.6/sqrt(.1)+130sqrt(.8/.69^2)<300.

For delta>=gamma, integration of 1/Q from c gives

    |T(z)|/E<=300+130sqrt(8/delta)<=700/sqrt(delta). (4)

Write s=Y(z)-Y(c), dY=dz/sqrt(Q), delta_c=z_h-c>.79.
While delta>=gamma, sqrt(Q)/delta<=sqrt(3.001)<2, so
delta>=delta_c exp(-2s). Therefore, for every k>=18,

    integral_(delta>=gamma) e^-ks |T| dY
                         <=800 E/(k-1).           (5)

No bolt point-evaluation norm is substituted here. On delta<=gamma use

    integral_c^z Q^-1 dt <=gamma^-1[8+.5log(gamma/delta)],
    dY<=d(delta)/sqrt(2gamma delta).

Using sqrt(A+B)<=sqrt(A)+sqrt(B) and the exact elementary gamma integral
int_0^1 t^-1/2 sqrt(log(1/t))dt=sqrt(2pi) gives

    integral_(delta<=gamma)|T|dY
       <=[300sqrt(2)+130(2sqrt(8)+sqrt(pi))/sqrt(2)]E/sqrt(gamma)
       <1110 E/sqrt(gamma).                         (6)

From Q<=1.001(delta+gamma)^2, the distance to this last region is
s_gamma>=log[(delta_c+gamma)/(2gamma)]/1.001. Its weighted cost is
at most 1110 E gamma^-1/2(2gamma/.79)^(k/1.001). The product with
k-1 decreases for k>=18; at gamma=1e-30, k=18 it is <1e-500.
Monotonicity in gamma gives the same upper bound for the actual smaller
positive gamma. Combining (5)-(6) proves the WHOLE physical-input bound

    integral_(Y(c))^(Y(z_h)) e^(-k[Y'-Y(c)]) |T(Y')|dY'
                               <=900 E/k, k>=18.    (7)

This includes the actual regular-bolt endpoint and all intermediate
geometry. Its purpose is to estimate a return from the coupled solution,
not to optimize the already controlled cold source. The singular value
majorant in (6) is integrable; no unproved bounded trace is assumed.

## 3. Apply the actual complete paired kernel, not a normal surrogate

Use the original first-order transfer and spin traces in
current-separated-quantum-graph-bound.md:

    M_n^2=n(n+q)/C, d_n=2(q+2n), p_j=gamma(j+1/2),
    Udot=W(1-U^2)-2M U, W=p/sqrt(Q), 0<=U,V<=1,
    Bdir=-4U_A V_B exp(-2I_AB)/(1+U_B V_B)^2,
    Bimg=-4U_A U_B exp(-2I_img)/(1+U_B V_B)^2,
    Hm(z,Y')=-2G gamma/(pi B^2 sqrt(Q(z)))
                         sum_(n>=1,j>=0)d_n M_n^2(Bdir+Bimg).  (8)

The invariant intervals, orientations, both thermal signs and physical
paired image normalization are precisely the retained actual realization.
The regular endpoint branches and boundary projectors are not varied.
For k=sqrt(M^2+5p^2), the Riccati comparison gives the original bounds
exp(-2I_AB)<=4 exp[-2k(Y'-Y)] and
exp(-2I_img)<=16 exp[-2k(Y'-Y)]. Image travel is no shorter than direct
travel. These statements retain the actual thermal/bolt/domain return.

The observation collar is now z<=.1, so check its analytic germ rather
than importing the smaller-collar guard: the actual polynomial obeys
|Q/beta-1|<.01 on the entire complex disk |z|<=.11. For each real
observation continue only its one-sided U,V germs on radius
1/[1000(M+p)]. Then |Q^-1/2|<4, |Q^-1|<16. The Riccati right side and
Lipschitz constant are <36(M+p) on |U-U_real|,|V-V_real|<=.1.
Picard preserves the .036 disk, |1+UV|>=.79, and
|partial_z I|<=18(M+p). The same Cauchy bound is consequently valid:

    |partial_z^j[(Bdir+Bimg)/sqrt(Q)]|
       <=640 j! [1000(M+p)]^j exp[-2k(Y'-Y)], 0<=j<=4. (9)

This is an observation derivative of the ACTUAL kernel. It does not
differentiate an unsigned bound, continue an arbitrary mass/log input,
or replace the cubic paired double by a smooth geometry.

The proper separation is at least d=sqrt(5)(c-b). Factor the exponential
in (9) at Y(c) and use (7) with 2k. It replaces the former hypothetical
supremum by 900 E with the same 1/(2k)<=1/(2M) integration factor.
Put c_M=sqrt(2)d and c_p=sqrt(10)d. The COMPLETE sums obey

    sum_n d_n M_n^r e^(-c_M M_n)<=4C(r+2)!/c_M^(r+2), r>=1,
    gamma sum_(j>=0)p_j^l e^(-c_p p_j)
                    <=l! c_p^(-l-1)(1+gamma c_p), l>=0.       (10)

The first follows from the complete counting function <=4CM^2.
For the second the midpoint sum differs from the integral of a positive
unimodal function by at most gamma times its maximum. The inequality
(l/e)^l<=l! supplies the displayed bound. The nonzero gamma factor is
retained; no discrete thermal sum is replaced by an integral equality.

Expanding (M+p)^j before summation gives

    sup_(z<=b)|f^(j)(z)| <=900 E A_j,
    A_j=2560 GC/(pi B^2) j!1000^j(1+gamma c_p)
        sum_(l=0)^j binom(j,l)(j-l+3)!l!
                            /[c_M^(j-l+3)c_p^(l+1)], j<=4.   (11)

The positive majorants (7),(9),(10) justify all four derivatives and
the common-regulator removal on this separated block. Local seagulls
and restorations have zero cross-support action, not zero full action.
Both original equations, the mass/log completion and regular bolt remain.

## 4. Original H3 output including every cutoff derivative

Use exactly the original Fourier H3(dx), x=log z. For chi's logarithmic
beta(11,11) cutoff the Bernstein finite-difference bound gives

    h_0=1, h_l=[21!/(21-l)!]2^(l-1)/log(2)^l, 1<=l<=4. (12)

Define C_00(alpha)=1 and
C_(m+1),j(alpha)=(alpha+j)C_mj(alpha)+C_m,j-1(alpha).
The product rule and whole-interval integrals in (2) give, for 0<=k<=3,

    L1_k=sum_(l=0)^k sum_(j=0)^(k-l) binom(k,l) h_l
           |C_(k-l),j(.995)| b^(.995+j) A_(j+1)/sqrt(2(.995+j)),
    L0_k=sqrt(log 2) sum_(l=0)^k sum_(j=0)^(k-l) binom(k,l) h_(l+1)
           |C_(k-l),j(-.005)| max_(a<=z<=b)z^(j-.005) A_j,
    C_far=900 sqrt[sum_(k=0)^3 binom(3,k)(L1_k+L0_k)^2].       (13)

All maxima in (13) are endpoint power maxima; there is no sampled
derivative or quadrature tail. The original-row transition is L0 and
has not been killed by applying (D-4) before the cutoff.

The supporting evaluation will test

    C_far<1.124e-10,
    ||g_far||H3 <=C_far*1.1e-10 <1.24e-20,
    |f(0)|<=900*1.1e-10 A0<1.27e-34.               (14)

The latter is the forcing of the original z^4 shape row, whose scalar
quartic coefficient is 20beta Fstar>3.3333332. Its isolated algebraic
conversion is <3.82e-35, but this is NOT the actual quartic coefficient:
all other coupled row entries still determine that coefficient together.

The independently retained auxiliary Mellin tail is different from
(14). Substituting the newly certified |m_linear|<1.627532e-15 into
its original bound 1.02e-46+6.50e-39|m| gives <1.020000106e-46.
Neither kind of distant return is the unknown near-diagonal graph.

## 5. Decision for the complete coupled problem

The attempted full localization now has a certified distant-input term
on the actual linear response, with no assumed global supremum, independent
phase amplitude or cold norm. The original mass conversion remains <2174
per unit mass. After reserving the conservative mass bound and (14),

    1e-11-2174*(1.627532e-15)-1.24e-20 >6.4617454e-12.       (15)

This is the remaining SUFFICIENT SHAPE-graph allowance, not the actual
near graph or a bound for the original scalar graph. The conversion
d_chi=B_sigma[z^-sigma chi(b_m-z^3)] used in (15) is the shape output;
no new scalar norm comparison is inferred from its value 2174.
The first still unsupported quantitative step is the common-restored
near-diagonal operator and its transition response on 0<z<.2, with both
original local equations and their finite matching rows. Neither (7)
nor the weak inverse controls the derivatives needed for this part.
Its fixed-seed qualitative graph theorem supplies no numerical constant.
No inference from .01 local fast estimates to the full .05/.1 transition
is made. The constants for independent phase columns and nonlinear
neighborhoods remain unknown.

The conditional mass-column target kappa_m<=5e-5 is not evaluated here.
Its costs 2.14765e-12 affine and .0024731 feedback, and residual budgets
Y_other<=5.34234999994958e-12 and
eta_other+(.0195239+4.583e-13)gamma_E<.0012069, remain conditional.
The screen .028913+1e-11*1.0706e10=.135973<.1599164 keeps only its
proved two-fast-input scope. No actual independent phase coefficient is
identified with its allowance. The exact slow pair and scalar quartic
remain in the unknown coupled near/transition solve.

RESEARCH DECISION: retain the construction. Use (1) and (7)-(14) for
the complete distant response and stop assigning it hypothetical input
amplitudes. Its certified cost is negligible in the original allowance;
further improvement there will not decide feasibility. The complete
coupled estimate remains UNRESOLVED, with no excessive upper bound or
actual obstruction yet established for its near/transition block.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED HERE): enclose the actual
two-row near/transition graph response of the same full linear solution,
including common-restored coincidence contacts and original scalar/slow
matching, against the remaining 6.4617454e-12 allowance. A coupled defect
certificate must bound its strong error, rather than substitute energy
or the two-fast-input screen. This is the unexecuted part of the requested
global calculation, not a new smaller objective.

Remaining to-do list: complete coupled graph and independent phase columns;
physical-norm/neighborhood comparisons; nonlinear admission/contraction;
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
