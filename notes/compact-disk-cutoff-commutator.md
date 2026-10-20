# A quantitative cutoff estimate for the complete disk Hessian

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **One-patch localization;
the global fast inverse and physical solution remain open.**

The preceding [complete compact graph](complete-compact-disk-graph.md)
requires supported inputs. This note makes its cutoff error explicit
for a global input. Derivatives of the input are needed only on a
larger interior patch; its far part enters through a global value norm.
The full massive kernel, returns and local contacts are retained.

## 1. Statement, cutoffs and norms

Use the preceding normalized analytic conformal chart of radius 1/32,
centered at Y>=1/4, including a possible regular bolt center. Write
O={|z|<1/128}, gbar=exp(2sigma)|dz|^2, theta=S+T/3 and w=(theta,T).
The operator L is the COMPLETE fixed-gbar-area disk Hessian, with the
same action, parameters, state and domain as the preceding graph.

Let P(s)=462s^6-1980s^7+3465s^8-3080s^9+1386s^10-252s^11.
For a>0 let c_a(z) be one through radius a, zero by radius 2a,
and 1-P(|z|/a-1) in between. Choose

    chi=c_(1/1024), eta=c_(1/128).

Extend both by zero outside the chart. In particular chi is supported
strictly inside V={|z|<1/256}, as required by the preceding theorem.
Use the Euclidean Fourier norms there, and put

    X4(v)=||v_theta||_H4+||v_T||_H4log,
    U3(v)=||v_theta||_H3+||v_T||_H3,
    N0(w)=||theta||_L2(gbar)+||T||_L2(gbar),
    Fchi=||chi L_theta w||_L2(O)+||chi L_T w||_L2(O).

The N0 norm is on the WHOLE physical half and includes time.
For smooth global inputs in the retained domain, initially compact
away from the seam, the following a priori bounds hold:

    ||1_O[L,chi]w||_(L2+L2)
        <=1e-10 C X4(eta w)+1e155 C U3(eta w)+1e150 C N0(w), (1)

    X4(chi w)
        <=(2000/C)Fchi+1e-6 X4(eta w)+1e670 N0(w).           (2)

The constants are deliberately loose. Equation (2) has a small
larger-patch derivative term; it is NOT yet an estimate independent
of that term. A quantitative global cover, outer weighted estimate
and their absorption still have to be supplied. In particular the
large value constant is not a nonlinear contraction.

The cutoffs are C5 with bounded piecewise sixth derivatives.
The statements follow first for smooth approximations preserving
the strict bounds below and then by convergence in H4log and in
the displayed commutator norms. The same estimates extend to the
corresponding graph closure when the source converges in L2(O).
No new seam or bolt boundary condition is imposed.

## 2. Bounds on the separated kernel and its full tail

Write H for the actual massive entry, the only nonlocal entry of L.
The [curved response](stationary-curved-quantum-tt.md), sections 3-4,
and [strong interior remainder](strong-interior-quantum-remainder.md)
give the local geodesic bubble envelope 72C/d^6, its comparison
error 1e8C/d^4, and a separate actual gluing/far kernel with both
Schur rows at most 3e24 C. These are kernel estimates, not a
strong estimate inferred from a quadratic form.

For x in O, convert the local part to the fixed coordinate chart.
On separation r=|y-x|<=1/1024, metric/distance ratios lie between
.9 and 1.2 and the density conversions are less than four. Thus

    |K_local(x,y)| <=1e12 C/r^6.                           (3)

The bounded gluing kernel stays separate; (3) is not asserted for
that kernel pointwise. The cutoff satisfies 0<=chi<=1 and
Lip(chi)<1e6. Multiplication gives the exact separated commutator
kernel (chi(y)-chi(x))K(x,y). Consequently its local part satisfies

    |K_comm,local(x,y)| <=1e18 C/r^5.                     (4)

Fix R=1e-30. Keep (4) for r<R. Its entire complement, including
the gluing kernel even at r<R, is an ordinary bounded operator B.
Both endpoint rows, with input gbar and output Euclidean measure,
are bounded by

    1e30 C R^-4 =1e150 C.                               (5)

For clarity, an ample bound before this enlargement is
100[3e24+4pi(72/(.9R)^4+1e8/(.9R)^2)]C. It follows by
integrating the proper-distance envelopes up to the local chart
radius and retaining the actual far/gluing rows. Reverse rows
are the same adjoint rows in the retained paired realization.
The factor 100 covers coordinate, half/double and commutator
conversions. The value of chi is at most one globally.
This also bounds B on probes of supremum at most one.

In particular (5) uses the actual far response. It neither replaces
it by a compact cap nor asks for derivatives of an input near the seam.

## 3. Near remainder and the order of its local contacts

For a fixed x, subtract the Euclidean degree-three Taylor polynomial
J3_x v from the input in the near kernel and define

    A_R v(x)=integral_(|y-x|<R) K_comm,local(x,y)
                                      [v(y)-J3_x v(y)] dy.

Taylor's integral remainder, translations, Minkowski and (4) give

    ||A_R v||_L2(O)
       <=(2pi/24)1e18 C R |||D|^4 v||_2
       <=1e20 C R ||v||_H4.                            (6)

For the point version a fourth Cartesian component bound costs
at most four more, still within the reserve 1e20. Only x within
R of supp chi contributes. There eta=1, including the whole
Taylor segment, so A_R v=A_R(eta v).

The remaining local distribution has input differential order at
most THREE. This requires the retained contact allocation, not
just (4). Before multiplication the full near massive finite part
uses Taylor subtraction through degree four, together with its
fixed ordinary fourth restoration. Multiplying that distribution
by chi(y)-chi(x), whose value at y=x is zero, removes the fourth
input derivative in the product Taylor polynomial. The ordinary
fourth differential restoration commutes to order three as well.
Every other retained contact already has order at most two.
Changing the subtraction radius or coordinates changes coefficients,
but not this input order.

Equivalently, if an input has all jets through order three zero at
x, the product with chi(y)-chi(x) vanishes through degree four,
so the original subtracted distribution is the absolutely convergent
near integral. Thus the difference from A_R+B vanishes on all such
inputs. It is a combination of the ten jets of degree at most three:

    [H,chi]v = A_R v+Bv+sum_(|alpha|<=3) a_alpha(x) D^alpha v(x).
                                                               (7)

This is the point-supported distribution principle in
[Dyatlov, Theorem 4.19](https://math.mit.edu/~dyatlov/18.155/155-notes.pdf#page=53),
inspected 2026-10-20. The order-three conclusion here follows from
the preceding vanishing-jet argument and the specified physical
restoration; the general structure theorem alone does not give it.
Local coefficients in (7) vanish where chi is identically zero
in a neighborhood. Hence their derivative terms also use eta v.

## 4. Bound every coefficient by an actual commutator probe

At every x in O use

    q_alpha(y)=c_(1/2048)(y-x)(y-x)^alpha/alpha!, |alpha|<=3.

Their first four jet levels at x form the identity matrix,
their suprema are at most one, and their supports remain inside
the chart of radius 1/64. The cutoff derivative ledger in the
strong remainder proof, rescaled to these radii, gives all Cartesian
component derivatives through six of chi and every q_alpha below
1e32. For chi the first bound is much sharper:
each component derivative is <=22*512, implying Lip(chi)<1e6.

Now exploit the EXACT identity

    [H,chi]q_alpha(x)=H v_x(x),
    v_x(y)=(chi(y)-chi(x))q_alpha(y),  v_x(x)=0.       (8)

The identity includes the density contact. In fact the fixed-area
Hessian is H v=D T_mass[v]+2v T_mass. The second term vanishes
for v_x at x; it has not been discarded for arbitrary inputs.

Here is an explicit uniform reserve for the first-Cauchy norm
used in the strong remainder proof. Leibniz bounds the old-chart
jets of v_x through five by 1e70. To pass to the normalized chart
centered at x, use its analytic conformal transition. On a disk
of radius 1/512 around every point of the probe support, both
charts exist and the transition satisfies |F'|<2: its derivative
is the ratio of the two conformal scale factors. The margin to
either radius-1/32 boundary is larger than 1/512 by first-exit
length comparison. Cauchy therefore gives

    |D^j F| <=2(j-1)!512^(j-1) <=j!10000^j, j<=5.

For two real coordinates the chain-rule bound at order j is

    1e70 sum_(k=1)^j 2^k Bell_(j,k)
                        (1!10000,...,(j-k+1)!10000^(j-k+1))
       <=1e70 (2*3^(j-1)) j!10000^j <1e100.          (9)

Indeed Bell_(j,k)(1!b,...)=j! binom(j-1,k-1)b^j/k!;
discarding 1/k! gives the stated upper bound, not an equality.

The value and actual/capped gradient norms fit the same reserve.
The probe support lies where the recentered cap equals one.
Consequently U(v_x)<1e100. The retained first-Cauchy source
bound, valid for these nonstationary test directions as well,
gives

    |[H,chi]q_alpha(x)| <6e122 C.                    (10)

The near action on a probe is at most 1e20 C R times its fourth
derivative bound, hence <1e22 C. The complete tail is <1e150 C.
Subtract both from (10). The identity jets then give

    |a_alpha(x)| <1e152 C, |alpha|<=3.              (11)

This evaluates all ten coefficients, including the zero coefficient.
There is no assertion that a zero-extended constant is annihilated.
No derivative of an unsigned coefficient bound is taken.
Equations (5)-(7) and (11) prove

    ||1_O[H,chi]v||_2
       <=1e-10 C ||eta v||_H4
          +1e154 C ||eta v||_H3+1e150 C ||v||_L2(gbar). (12)

## 5. Complete local terms and application of the graph

For a coordinate differential term a_alpha D^alpha, Leibniz gives

    [a_alpha D^alpha,chi]v
       =a_alpha sum_(0<beta<=alpha) binom(alpha,beta)
                                      D^beta chi D^(alpha-beta)v.

Thus the complete lower local coefficient sum 1e60 C in the
preceding disk proof contributes at most 7e92 C U3(eta w).
Its two fourth diagonal terms contribute less than 1e35 C U3.
There is no derivative of a coefficient in this calculation:
the preceding proof already expanded the full divergence operators
in nondivergence form, including their completion terms.
Enlarge the local sum to 1e94 C U3 and add (12). This proves (1).

The exact Fourier expansion of (1+xi1^2+xi2^2)^3 has positive
coefficient sum 27. Leibniz and the cutoff jets therefore imply

    U3(chi w)<=1e40 U3(eta w).                        (13)

Apply the intermediate preceding graph to chi w, using
L(chi w)=chi Lw+[L,chi]w. Equations (1) and (13) give

    X4(chi w)<=(2000/C)Fchi+2e-7 X4(eta w)
                        +1e160 U3(eta w)+2e153 N0(w). (14)

Use U3(v)<=epsilon X4(v)+(2/epsilon)^3 U0(v), with
epsilon=1e-167, and U0(eta w)<=2N0(w). The derivative coefficient
is at most 3e-7<1e-6, and the resulting value coefficient is
at most 16e661+2e153<1e670. This proves (2).

For the physical variables and fixed physical disk source f,
the preceding measure and transpose formulas give

    Fchi <=1e4 B^2 C^-1/2 Fdisk,
    N0(w)<=C^1/2 B^-2 ||(S,T)||_L2(dmu_4),

where Fdisk=||u^-2 f||_L2(dmu_4), u=B/r. The second inequality
uses r>=B globally and the two-row transpose factor 5/3.
Thus (2)'s source term is at most 2e7 B^2 C^-3/2 Fdisk.
For the already constructed compatible linear solution the accepted
global value bound can be substituted as well. Neither substitution
removes X4(eta w) or closes the weighted outer estimate.

## 6. Frozen verification and limitations

Freeze this proof and checker before execution. Pin both preceding
receipts and their complete input hashes. Check the exact cutoff
jets and probe jets, two-dimensional identity jet matrix, transition
and Bell reserves, local kernel conversion, both tail rows, Taylor
integral formula, finite-product fourth-jet cancellation, density
cancellation on (8), all differential commutator coefficients and
orders, H3 multiplication and the final interpolation arithmetic.
Preserve every scientific run and any failed source snapshot.
Finite checks support the ordinary operator/domain argument;
they are not a formal verification of the continuum theorem.

The first retained run rejected an incorrectly stated equality
in (9): its upper envelope had been written as the exact Bell sum.
The corrected statement is the displayed inequality; the exact
sum retains 1/k!. That failed run and both failed source snapshots
remain in the receipt. The numerical reserve and theorem constants
are unchanged.

This resolves the explicit cutoff error on a single interior chart.
Summing a cover with a controlled overlap, higher York transfer,
the outer fast norm and complete inverse/residual composition remain
necessary. No nonlinear branch, signed compatibility solution or
physical stability follows from this localization estimate.

Remaining to-do list: global compact patching and higher York maps;
outer fast fourth/log inverse and residual composition; nonlinear
signed compatibility; exterior/state, opening, stability, transport,
preparation, separate creation and full-result audit.
