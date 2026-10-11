# The vacuum reference needs the interaction across the mouths

2026-10-11, S053. A positive photon exterior and its relatively bounded
clock derivatives do not suffice to renormalize the throat-only vacuum
reference with a bulk subtraction. An exactly soluble member of that
class has a linear ultraviolet boundary term, including divergent clock
derivatives. In a uniform full-line model the omitted cross-mouth
interaction cancels the term exactly. The next task is therefore the
complete mouth interaction and its renormalization matching, before
extending the throat-only calculation to higher current diagrams.

This is an obstruction to a proposed inference from the existing
projected model. It is not a physical wormhole instability or a new
Casimir principle. The finite return-dependent difference in the
[previous note](winding-interaction-stress.md) is unaffected.

## 1. The reference and a positive exterior realization

At Euclidean frequency k>0, take one electric multipole on I=(0,T),
with constant positive core mass m and angular coefficient kappa.
Let L u=(u',-ik u) and let R be the restricted photon dual-scalar
resolvent. The potential propagator and the unwrapped current are

    D = kappa^-1 (1-L R L*),
    Pi_infinity = (tau/pi) P,   P = L S L*,
    S(s,t)=exp(-k|s-t|)/(2k).                              (1)

P is the restriction of the full-line massless transverse projection,
including the same local vector contact as before. The L* in (1) includes
the endpoint distributions. Define the normalized reference integrand

    F(k)=Tr_I[(1-L R L*) P],
    E_ref,Lambda^(2)=g^2 tau/(2 pi^2 kappa) int_0^Lambda F(k) dk. (2)

For each k>0 this trace is finite. Its integral need not be finite.

Attach two massless half-lines to the core. In the common continuous
field variable, their Euclidean form is

    A_ext,k = h int_ext (|u'|^2+k^2|u|^2) ds,  h>0.        (3)

The kinetic Hilbert norm on each exterior has the same weight h.
Eliminating its decaying solution gives the endpoint form z=k h at
both ends. This is one positive realization at every frequency, not
independently assigned Robin values. The core form is
int_I(|u'|^2+(k^2+m^2)|u|^2). Thus it obeys core energy domination
and has no global gap. It can also dominate the earlier AdS core on
a fixed truncated interval by choosing m sufficiently large. A single
angular mode already suffices for the counterexample.

For the clock family T=T0 exp(-phi), m=m0 exp(phi), h=h0 exp(phi),
pull I back to (0,1) preserving endpoint values. The coefficients of
the core gradient, core mass, core k^2 term and exterior form are
positive multiples of exp(phi) or exp(-phi). Their first and second
derivatives are bounded in absolute form by the whole form itself.
The domain and endpoint map are fixed. On a compact phi interval
this family satisfies the positivity and two-derivative assumptions
used for the finite winding diagram. It is a counterexample within
that declared mathematical class, not a derived moving-mouth geometry.

## 2. Exact trace and ultraviolet coefficient

Put r=sqrt(k^2+m^2), rho=(r-z)/(r+z), e=exp(-kT), and

    g_plus  = 1/[z+r tanh(rT/2)],
    g_minus = 1/[z+r coth(rT/2)].                          (4)

These are the symmetric/antisymmetric eigenvalues of the boundary
Green matrix G=(r[[coth(rT),-csch(rT)],[-csch(rT),coth(rT)]]+z I)^-1.
The diagonal Green integral is

    Tr R = T/(2r) (1+rho^2 exp(-2rT))/(1-rho^2 exp(-2rT))
           + rho(1-exp(-2rT))/(2r^2(1-rho^2 exp(-2rT))).   (5)

Direct contraction of the two potential components gives

    F = m^2 Tr R + z(g_plus+g_minus)
          + (k/2)[(1-e)g_plus+(1+e)g_minus] - 1.           (6)

The last three terms cannot be dropped by an interior Ward identity.
One way to verify (6) without cyclically rearranging divergent traces
is to use the exact core Green kernel

    R(s,t) = [exp(-r|s-t|)+rho exp(-r(s+t))
              +rho exp(-r(2T-s-t))+rho^2 exp(-r(2T-|s-t|))]
                /[2r(1-rho^2 exp(-2rT))].                (7)

In distributional derivatives R_st=delta(s-t)+R_st,reg. The contact
terms and ordinary smooth integrals in (2) are precisely

    int_I [-R_st,reg(s,s)+k/2] ds
      + (k/2) int_(I x I) exp(-k|s-t|)
          [R_st,reg-k^2 R+k sign(s-t)(R_s-R_t)] ds dt.     (8)

Substitution of (7) and integration of exponentials yields (5)-(6).
Equation (8) makes the delta cancellation and the remaining finite
trace explicit; no product of two unresolved delta distributions is
assigned a value. It also supplies an independent numerical check.

With z=hk, as k tends to infinity,

    F(k,phi)=c(h)+m^2 T/(2k)+O(k^-2),
    c(h)=h/(1+h).                                        (9)

The expansion holds with its first two phi derivatives uniformly on
compact clock intervals: r/k=1+O(k^-2), the denominators r+hk are
bounded below by a positive multiple of k, and the remaining exponentials
and their first two derivatives decay for T>=T_min>0. This also proves
integrability of the twice differentiated remainder. At k=0 the massive
core Neumann resolvent is finite and z=hk is smooth; (6) has a finite
limit, so there is no compensating infrared divergence.

Subtract the usual local bulk term with any fixed mu>0:

    B_Lambda(phi)=int_0^Lambda
            [F(k,phi)-m^2 T/(2 sqrt(k^2+mu^2))] dk.        (10)

Then, for j=0,1,2,

    D^j B_Lambda = D^j c(h) Lambda + b_j(phi)+O(Lambda^-1). (11)

For h=2 exp(phi) at phi=0,

    c=2/3,   Dc=2/9,   D^2c=-2/27.                      (12)

Thus the energy, first clock derivative and static curvature in this family
all have nonzero linear divergences after bulk subtraction. A fixed
finite bulk coefficient or changing mu cannot remove them. The finite
winding correction is exponentially convergent at high k and cannot
cancel (11). Taking g small before removing a regulator does not prove
a regulator-independent result at a fixed nonzero g.

The linear coefficient is stated for the explicit frequency cutoff (10).
Power divergences depend on the regulator; a prescription that drops them
does not by itself establish the physical matching or fix a finite
interface energy. The counterexample disproves sufficiency of the stated
positive-form hypotheses with this bulk prescription, not the existence
of a consistent renormalization of the complete theory.

If every angular channel has this same high-frequency h, its total
linear coefficient, using the exact LLL weights, is

    g^2 q(H_q-1)/(8 pi^3) c(h) Lambda.                    (13)

Formula (13) concerns the projected model, with q finite before taking
Lambda to infinity. It is not a four-dimensional divergence calculation
or a prediction at an unspecified physical cutoff. Fixing h under clock
variation removes this particular divergent clock coefficient; that
additional restriction does not follow from positivity, and does not
fix the full energy or physical matching.

## 3. Exact cancellation when the exterior interaction is included

There is a second useful check with a uniform massive photon on the
whole line. Now R=exp(-r|s-t|)/(2r), so the exact exterior load is z=r
and rho=0. The fermion current is still the full-line massless one.
For both vertices restricted to I, (6) becomes

    F_II = m^2 T/(2r) + b(k,T),
    b(k,T)=k/(2r)[1-exp(-(r+k)T)].                        (14)

At high frequency b tends to 1/2, independent of T. The full-line
Fourier symbols instead give

    Tr_components[D(p,k) P(p,k)] = m^2/(p^2+k^2+m^2)

after removing the factor 1/kappa. Its spatial integral assigns the
bulk value m^2 T/(2r) to an interval of length T.

The missing term can be exhibited without an infinite-volume trace.
For the ordered contribution with s in I and t outside I, the smooth
contracted kernel is

    - k(r+k)^2/(4r) exp(-(r+k)|s-t|).

Since int_I ds int_outside dt exp(-a|s-t|)=2(1-exp(-aT))/a^2,

    F_I,outside = -b(k,T),
    F_II + F_I,outside = m^2 T/(2r).                     (15)

Both integrals are finite at each k. Equation (15) is the local trace
identity with its first coordinate in I and the second on the full
line; it is not an assertion that the full infinite-volume energy is
finite. In a full partition, both ordered cross terms are included.
They cancel the two interface contributions from the diagonal blocks.

This exact example explains why adding a passive photon load while
leaving the current interaction supported only in the throat is not
enough to identify the reference vacuum stress. The missing interaction
across that artificial split can be essential. It does not imply a
universal cancellation for the actual curved two-mouth system.

In current language, restricting a conserved full current with the
indicator 1_I creates endpoint flux terms in its divergence. The earlier
projected susceptibility remains a legitimate response to a supported
source; using it alone as the entire vacuum interaction introduces an
extra assumption. Restoring the full current and any other required
four-dimensional vertices is a physical matching problem.

## 4. Source comparison, verification and decision

[Graham et al., hep-th/0207120v2](https://arxiv.org/html/hep-th/0207120v2),
sections 3.2-3.4 and 5, distinguish renormalization with fixed physical
inputs from singular boundary limits. Their example has a divergent
boundary energy while the rigid-separation force is finite. This is
established precedent for testing the actual varied observable and
retaining source details. Their scalar background is different from
the current/photon projection here; (6)-(15) are derived above, not
imported from that example. An arbitrary finite boundary subtraction
would add missing physical data, not complete the four-dimensional
renormalization. The source supplement records the limited review.

The [check script](../checks/reference_mouth_completion.py) and
[receipt](../evidence/reference-mouth-completion.json) compare direct
contact-resolved kernel integration with (6), verify the asymptotic
coefficient and clock derivatives, and check the ordered cross-mouth
cancellation against both real-space integration and the full-line
Fourier integral. Declared refinements test the numerical implementation;
the exact identities establish the obstruction. No Lean/Lake was run.

**Decision:** keep physical long-return damping as the primary and matched
finite-duration energy transfer as the fallback. Change the next action:
derive the complete interaction through the mouths and specify its fixed
renormalization conditions before further throat-only higher-diagram
work. The earlier quadratic response and finite winding estimates remain
usable components. Neither the present counterexample nor an arbitrary
surface subtraction closes physical stress, stability or novelty.

Remaining to-do list: derive the full mouth interaction and renormalized
reference, control higher-order and dynamical stress, match the coupled
geometry and state, complete novelty and the consequential result, then
independently verify, GCP-certify suitable endpoints and prepare the manuscript.
