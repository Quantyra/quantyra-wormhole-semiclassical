# The cold quantum column along the current reference correction

2026-10-22, E003/S169. Ordinary argument before frozen verification.
This executes a directional source test within the already selected
finite-domain defect calculation. It does not reopen the closed audit or
certify the actual coupled correction. The full stable-wormhole objective
remains incomplete.

## 1. Combine the current local sources before taking absolute values

Retain the current action, c_R=1e22, q=2000^4, g_R=2000^-4,
mu0=1/b0, positive prepared temperature, physical clock and paired HH
domain. The symbols and source signs are those of
[the full trace calculation](current-seed-full-trace-defect.md).
Put z=B/r, x=lambda B^2, u=1/x, chi_R=16piG c_R lambda,
e_Q=Delta Q/(lambda B^4), and beta=delta/x. Then

    f=F/x=z^-2+u-(4+2u+beta)z+(3+u+du)z^2-2du z^3+du z^4,
    R_n=R4/lambda=-12+du(4z^5-6z^6),
    T_n=tau/lambda=(3+u+du)z^4-du(5z^5-9z^6/2),
    d_t=-z^2 d_z,
    H_r=H^r_r/lambda^2=-2R_n T_n+2[(d_t f)/2+2zf]d_t R_n,
    H_t=H_r+2f d_t^2 R_n.                              (1)

The exact effective charge relation is Q_phys=(1-24chi_R)Z and
Q_seed=Z+Delta Q. Consequently

    (Q_seed-Q_phys)/(lambda B^4)
       =24chi_R(3+u+du)+(1-24chi_R)e_Q.

Combine R2 with the NONVACUUM reference subtraction. The vacuum is still
combined with the quantum vacuum, as in the exterior source proof. Define

    J=[24chi_R(3+u+du)+(1-24chi_R)e_Q]z^4-4du z^5+3du z^6,
    P_rho=G(rho_R2-rho_ref,nv)/(lambda r^2)=(chi_R H_t-J)/2,
    P_p=G(p_R2-p_ref,nv)/(lambda r^2)=(-chi_R H_r+J)/2. (2)

The entire leading 12chi_R(3+u+du)z^4 cancels. For d=0,
P_rho=-(1-24chi_R)e_Q z^4/2 and P_p=-P_rho, exactly.
Every remaining curvature term contains chi_R du. This is a physical
cancellation in the same source and normalization, not a rescaling.

For a polynomial in z let A(p) be its absolute coefficient sum, bounding
each parameter by

    u<=1/(8e9), |d|<=8.6e-19, chi_R<=4.91e-10,
    |e_Q|<=3.4e-22, 0<=beta<=7h_max^2,
    h_max=2000^-10, h_min=h_max/8.

Exact polynomial arithmetic in (1)-(2) gives, for both profiles,

    A(P)<1.701e-22,
    A[-z^3 d_z(P/z^2)]<3.401e-22.                     (3)

The second expression is precisely (G/lambda)rho'/r or its pressure
analogue. Restore the quantum, Weyl and LLL derivative contributions:

    Qrho < (31e15+1e6)1.001/2000^8+3.401e-22 <1.3e-10,
    Qp   < (29e15+1e6)1.001/2000^8+3.401e-22 <1.3e-10. (4)

The earlier 1e-7 was a valid, insufficient upper bound obtained by
separating these local terms. Equations (2)-(4) replace it for this use.
They do not say the old actual correction was 1e-7.

## 2. Check the spatial and boundary hypotheses for the actual path

An energy first-derivative improvement alone would not establish C5.
Here the added local terms can be checked through all required source
jets. Since Z/(lambda B^4)=3+u+du-e_Q>2.99 and 1-24chi_R>.99999998,

    H_q=lambda B^4/(Gq^2)
       =pi/[g_R^2(1-24chi_R)(3+u+du-e_Q)]
       <1.052(2000^8).

The normalized local densities H_q P/z^2 are polynomials of degree at
most ten, with minimum degree two and coefficient sum <4.59e4. Hence
their absolute exterior mass integrals are <4.59e4 q^2/B. This is a
replacement for the local allocation, not another separate integral of
the uncancelled quantum vacuum or reflected pole.

For completeness the radial tensor-jet calculation does not divide by
the vanishing bolt warp. Set b=calB/C, p=sqrt(b), z_rho=-p. For any
radial polynomial L(z), the nonzero orthonormal Hessian entries are

    L_;rr=b L_zz+b_z L_z/2, L_;theta theta=b_z L_z/2.

The ordered third entries are

    L_;rrr=-p[b L_zzz+3b_z L_zz/2+b_zz L_z/2],
    L_;r theta theta=-p[b_zz L_z/2+b_z L_zz/2],
    L_;theta r theta=L_;theta theta r=-p b_z L_zz/2.   (5)

Use b<=1, |b_z|<=5, |b_zz|<=3 globally. On s=r/B-1<=.001,
b<=s^2, |b_z|<=3s, |b_zz|<=3. These follow from the retained
factorization calB_ext=(1-z)^2[x(1+2z+3z^2)+z^2+dz^4]
and the thermal term -delta z^3. In particular the first derivative
has coefficient at most 2(1-z)+3delta/C, below 3s; no bound on
N'/N at the bolt is needed. Equation (5), the degree bound, and
division by j!100^j give a reserve 1e8 q^2 a_s through j=3 for
the local difference from its product value, a_s=s[1+log(.01/s)].
The thermal offsets of P_rho,P_p at z=1 are respectively
108chi_R du beta and 24chi_R du beta. After multiplication by H_q
and division by h_min their absolute bounds are <1e-39, negligible
within this reserve. The total product value is cancelled by the CURRENT exact
product equations; the local term alone is not assigned zero.

In the outer collar these polynomials start at z^2. The already retained
analytic proper map has z/Y<.412 and jet multiplier <1.021. Subtract
the exact pole before integrating the canonical source; the new local
polynomial adds <1e8 to the existing 1e20 normalized source-jet reserve.
It creates no new mass boundary condition or nondecaying term. The
full trace polynomial in current-seed-full-trace-defect.md likewise
has exact product zero and bounded rational jets. Its coefficient
sum, after division by z^2, is <1e-24 in lambda units; (5) supplies
its three divided proper jets. The previously retained local curvature
forcing reserves are therefore ample, including the nonconstant R4
Lie term. These statements use finite spatial derivatives, not a
parameter derivative of a spatial source bound.

The current massive canonical source inputs are 5e17 q^2 a_s in the
cold C3 norm and 6e15 q^2 in the middle C3 norm, from
[the cold audit](cold-source-feasibility-audit.md) and
[the fixed-cap proof](fixed-cap-seed-source-three.md). They are stronger
than the reconstruction's 1e20 and 1e26 allowances. Its outer 1e20
allowance retains the pole, common finite restoration and correlated
mass cubic. The current geometry satisfies the same rational family
guards and actual paired-domain hypotheses. The scale exception
|log(mu0 B)|<1e-8 has already been substituted in
[the current seed source test](fixed-clock-seed-bolt-feasibility.md);
the restoration remains inside its old strict margins. Do not reset mu0.

The reconstruction also needs an ABSOLUTE lapse integral. It is not
inferred from the lapse value. With K=1e15 the actual Ward primitive
obeys |J(t)|<=Kq^2 F(t)(t-v)/(vt). Thus

    integral_v^infinity |a_mass,t|dt
      <=GKq^2/(.99 B^2 x v) integral_v^infinity dt/[t^2(t-1)]
      =GKq^2/(.99 B^2 x v)[log(v/(v-1))-1/v]
      <[80(1.001)/.99]K/2000^8 <3.2e-10.             (6)

The same positive integral for |Theta_LLL'|<=qC/4 gives <3e-28.
Differentiating the explicit Weyl lapse, rather than bounding its value,
gives total variation at most

    (G lambda)[(2.01+(19/8)3.01)/(15pi)
                        +3|d|/(2pi x)] <2G lambda.

The R2 lapse adds <9.182e-37. Their sum is <3.3e-10. Together
with the current horizon bound <2e-12, (4), the vacuum-recombined mass
bound <1e16 q^2/B, and the spatial inputs above, this supplies every
input of [the first reconstruction](stationary-centered-first-order-geometry.md)
and [its radial C5 transfer](stationary-completed-metric-response-norm.md).

That transfer uses the invariant curvature/radius Volterra equations:
three radial source derivatives give five metric derivatives, including
the regular integral at the bolt. It does not invoke a generic endpoint
C3-to-C5 elliptic sup estimate. Consequently the CURRENT specified
first reference correction h1 has the same bounds

    |a|<3.3e-10, |fbar/F|<7.5e-10, |H|<1.08e-9,
    |u|<1.2e-8, |du|proper<7.7e-6, |u_chi|<3e-6,
    |V|<1.2e-8, |s V_t|<3e-7,
    cold normal C5: orbit 1e-6 a_s, radius 1e-5 a_s,
    middle normal C5: orbit and radius 1e-4,
    outer weighted C5 after correlated cubic: each 1e-7. (7)

The derivative convention in the last three lines is j!100^j.
The cold barrier also gives |u|<=1e-6 a_s and |u_chi|<=8e-6 a_s.
One-center mixed jets cost the retained additional factor 100. The
fixed-chart inverse normalizer gives P5_40(u)<1.3e-4. The outer map
is weighted C6 on the OPEN collar, only C4 at its endpoint; allowed
Y^4 log Y tails are not promoted to closed-collar C5 regularity.

The comparison path is sigma_ell=sigma0+ell u,
r_ell=B+(r0-B)exp(ell V). Its endpoint proper identification is made
once for each ell and is identity on the cold region. The physical
paired boundary data, clock, period and cubic relation are preserved.
No parameter differentiation of a moving C5 coordinate map is used.
This proves path admissibility for this explicit h1, NOT admissibility
of the unknown full correction or a ball in the solver's X norm.

## 3. Take the endpoint derivative at the seed before spending log losses

Let S_ell=2G Theta_ell/(B^2 C) be the massive canonical scalar PER
DEFORMED AREA. The repaired actual-domain endpoint theorem in
[finite cold gluing](finite-cold-local-gluing-and-admission.md), with
d=0 and epsilon_end=1.000001(1e-6 ell), gives on this path

    |S_ell-S_0| <=1.71e-15 ell a_s,
    |(S_ell-S_0)(chi)-(S_ell-S_0)(0)|
                   <=1.71e-14 ell h L_h sqrt(chi), 0<=chi<=1. (8)

It includes the full subtracted angular tower, fixed common finite
restoration, regular bolt and actual finite paired domain. A cap inverse
or an unrenormalized operator trace cannot replace that proof. The
current geometry, fixed scale and path guards checked above are its
hypotheses. This use remains dependent on those ordinary continuum
kernel arguments; a numerical checker does not prove them.

The fixed-domain subtracted resolvent family is pointwise holomorphic
near ell=0, as in [the interior path proof](comparison-path-interior-response.md).
Divide (8) by ell and take ell down to zero. The result is

    |D S_0[u]| <=1.71e-15 a_s,
    |D S_0[u](chi)-D S_0[u](0)|
                   <=1.71e-14 h L_h sqrt(chi).          (9)

This uses more information than a uniform small bound on an interval:
the endpoint difference is proportional to ell. It avoids the
logarithmic losses in the old 3.5e-12 and 3.5e-11 uniform derivative
bounds AT THE SEED. It is not an improvement for all ell or for D2S.
Pointwise convergence plus (8) proves the bounds on the derivative;
it does not by itself prove convergence in the square-root norm.
Uniform holomorphy away from chi=0 and the uniformly vanishing tail
after division by chi^(1/4) give the weaker quarter-power norm limit,
as in the earlier path proof. The cold lapse integral is dominated
by chi^(-1/2), so it also has a valid first-derivative limit.

For the action's fixed SEED density the row is e^(2ell u) S_ell.
Its derivative includes the contact 2u S_0. The anchor and gradient
give |Theta|<=1.0001e15 q^2 on the cold patch. Since
2Gq^2/(B^2C)<2.72e-27, |S_0|<2.721e-12. Hence

    |2u S_0|<5.442e-18 a_s,
    |D[e^(2ell u)S_ell]_0|<1.72e-15 a_s.             (10)

For the bolt oscillation, write the contact difference exactly as
2[(u(chi)-u(0))S_0(chi)+u(0)(S_0(chi)-S_0(0))].
On chi<=1, a_s<=2hL_h and |u_chi|<=8e-6 a_s. The first
term is <8.708e-17 hL_h chi. The current canonical C3 gradient
and the cold proper-coordinate factor <1.2 bound the second term
by 1.4e-12(hL_h)^2 chi. With hL_h<81h_max this fits the
remaining 1e-16 allowance. Therefore the density-correct version is

    |Fdot_mass(chi)-Fdot_mass(0)|<1.72e-14 hL_h sqrt(chi). (11)

This executes a massive-canonical fixed-density cold VALUE and bolt
oscillation bound along h1. It does not supply three spatial derivatives
of this response, the independent radius/Weyl rows, or the other
off-shell row/pullback contacts in the full physical action defect.

## 4. Closed audit: what this closes and what it cannot certify

Closed by this calculation: the previously missing current local
R2/reference cancellation; the complete Qrho,Qp input; the absolute
lapse integral; the current explicit h1 path's C5 source-response
hypotheses; and the massive directional cold value/oscillation column
including its canonical density contact. The retained scalar-only
trace error <1.467e-13 and local reference-bolt feedback <3.927e-18
are unchanged. None is the full physical coupled error.

The numerical feasibility requirements remain, in the SAME solver norm,

    rho=1e-8,
    Y0+Z0 rho+Z2 rho^2/2<=rho, Z0+Z2 rho<1,
    Y0<=5e-9, Z0<=.25, Z2<=5e7;                      (12)
    1.632e-11 K_m+1.36e-9 K_c+Y_rest<=5e-9,
    K_m<=61.2745098039, K_c<=1.4705882352,
    Y_rest<=2e-9                                      (13)

as one sufficient allocation of (13). The two displayed source
coefficients in (13) are canonical SEED C3 norms. Equations (10)-(11)
are different quantities. Replacing 1.36e-9 by 1.72e-15 in (13)
would be invalid, even though the latter number is smaller.

Preservation also requires

    b45+rho E45<=1e-3 (E45<=87000 if b45<=1.3e-4),
    b_c+rho E_c<=5e-7 (b_c<=4.9e-7 if E_c<=1),
    b_o+rho E_o<=1.5e-7 (b_o<=1.4e-7 if E_o<=1).

The stronger finite-centered cold D_ell ball has radius 2.5e-10;
embedding the entire rho=1e-8 solver ball requires a distinct constant
E_D<=.025 when centered at the admitted comparison path. These
conditions concern the ACTUAL correction. The path contracts in (7)
are not their substitute. No finite N2-to-raw-C5 embedding exists;
the weighted derivative/domain requirement cannot be closed by energy
coercivity alone.

K_m,K_c,Y_rest, the required solver embeddings, the full directional
defect return and usable nonlinear Z constants are unknown in these
norms. In contrast the parked generic source-times-inverse estimates
are KNOWN upper bounds that fail: C_F=1.06963e201 versus the
2.38393e-25 sufficient zero-added-c threshold, and unsigned product
gain >301.92. Neither is an observed inverse value or nonexistence
theorem. The current action's unit form inverse is in V*->V; treating
it as a weighted C5 inverse would conceal the principal unresolved loss.

The massive first-order graph bounds, the full renormalized source map,
and the gravitational boundary realization are separate claims. The
matter projectors do not fix gravitational constraint boundary data.
The current stationary mass-extended graph realization has an ordinary
argument but no numerical weighted admission constant. Its physical
mass phase and log tails must be retained in the direct solve. Retarded
constraint boundary realization is still open; no growing or stable
mode of the actual finite geometry has been proved. Matching, exact
compatibility, semiclassical validity and all-sector stability remain
obligations; no boundary-free normal model can discharge them.

Decision: continue the SPECIFIED direct signed defect/adjoint method,
and stop cold source-value optimization. The cancellation and (10)-(11)
make the local/directional source route quantitatively plausible. They
do not make the unknown weighted return constants plausibly order one
by fiat; there is currently no numerical evidence that all of (12)-(13)
and the admission margins are attainable. Do not extend the generic
derivative chain or change construction solely because an unsigned
upper bound is large.

The single next calculation is the ACTUAL finite-domain physical defect
solve H e=-(F_eq+Hh1), equivalently its adjoint bolt evaluation WITH a
weighted remainder. Include the full quantum rows, canonical and physical
density/pullback contacts, common subtractions and physical boundary
phase in that calculation. Require coupled cold difference <1.667e-9,
shape-to-trace forcing <=1.25014e-27 lambda in the retained weight,
total Y0<=5e-9 and the displayed coefficient admissions. The full solve
is STILL UNEXECUTED. Equations (1)-(11) are an executed directional test
and reduce its unknown input; they are not a claimed solution of it.
If the actual return violates a target, test that signed result against
an alternative center/construction; if only an enclosure is too large,
sharpen the enclosure rather than assert nonexistence. S169.2-S169.5
remain unchecked and S166 remains the last completed component.

## 5. Frozen supporting verification

Freeze this note, checker and source review before running. Check exact
charge/source cancellation, polynomial coefficient and derivative
budgets, local tensor-jet and mass allocations, the positive lapse
integral, the transferred numerical contracts, endpoint limit and
density arithmetic. Use independent symbolic curvature/coordinate
identities and a manufactured fixed-density endpoint family; keep that
fixture distinct from an evaluated quantum determinant. Preserve all
runs, frozen source snapshots, input hashes and raw rows, including
failures. No Lean is selected. These checks do not replace the ordinary
actual-domain source and radial reconstruction arguments above.

The first supporting run failed an auxiliary thermal-offset assertion
<1e-50: its rational coefficient bound was 8.396e-41. The exact
offsets just displayed justify the corrected <1e-39 reserve, still
irrelevant beside the 1e8 local C3 allocation. This changes no claimed
path or response constant. Preserve that run and its original frozen
sources; replay the complete checker after this correction.

Remaining to-do list: complete the full physical defect solve with
weighted error and admission; preserved nonlinear solution, compatibility,
matching, retarded boundary realization, semiclassical validity,
all-sector stability and every remaining stable-wormhole requirement.
