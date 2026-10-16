# The actual finite-cut boundary Weyl functional

2026-10-15, E003/S139.1. Ordinary derivation, independently reviewed
before focused verification.
This note computes the parity-even, local boundary Weyl transgression
for the retained static paired Dirac domain. It includes the conversion
fluxes needed when using the accepted bulk representative. It does not
infer a physical charge from an on-shell boundary action value, construct
the full stationary stress, or close S139.

## Domain, normalization and what is being integrated

Use the identical classical magnetic RN exteriors and fixed bundle map
of [the conformal-domain calculation](stationary-exterior-conformal-domain.md).
More generally, put the smooth static spherical metric in the form

```
g = exp(2sigma) ghat,
ghat = dx^2 + N(x)^2 dU_E^2 + b^2 dOmega_2^2,
g_t = exp(2t sigma) ghat, 0<=t<=1.
```

Here x increases toward the outer boundary and is proper distance in
ghat; b is fixed. All calculations initially have a finite spherical
outer cut. A regular Euclidean bolt is not another boundary. Write

```
a=N'/N, s=sigma', Y=b^-2, R2=-2(a'+a^2), R0=R2+2Y.
```

Primes in this note denote d/dx, unless explicitly labeled otherwise.
The physical projector and its induced square are

```
chi=tau_3 gamma_out, Pi_+=(1+chi)/2,
Dom D: Pi_- psi=0, Dom Ddagger: Pi_+ phi=0,
L_+=Ddagger D: Pi_- psi=0,
                    Pi_+(nabla_out+K/2)psi=0.
```

The first-order operator is not declared self-adjoint on Dom D. In the
inward heat convention, L_ab=K_out,ab and S=-K Pi_+/2. In particular
chi_:a is generally nonzero. Neither a component Dirichlet condition
nor an independently chosen MIT condition replaces this domain.

The determinant modulus requires both positive partner squares. In a
fixed spinor L2 space, conformal covariance gives
`D_t=exp(-t sigma/2) D_0 exp(-t sigma/2)`. Consequently

```
delta L_+ = -(sigma L_+ + L_+ sigma)/2-Ddagger sigma D,
delta zeta_(L_+)(z)
 = z[Tr(sigma L_+^-z)+Tr(sigma L_-^-z)], L_-=D Ddagger.
```

This follows first where the resolvents are trace class, using
intertwining, and then by meromorphic continuation. Thus the Euclidean
fermion action `Gamma_E,mod=-log(det L_+)/2` has Weyl derivative equal
to the average of the two smeared a4 coefficients. Possible kernel
terms are separate finite spectral data; this local ultraviolet
calculation does not silently absorb them.

The partner exchanges Pi_+ and Pi_-, with the same K/2 Robin term on
its allowed trace. For scalar smearing, all traced a4 terms below are
unchanged: odd chi traces vanish, each projector has half rank, and
chi-derivative products are even. Thus their average equals either
computed coefficient in this static identical-pair problem. This is
an actual-domain calculation; coincidence with another boundary
condition's traced coefficient would not identify the domains.

Let r_f=4 for one exterior. For the identical pair, either sum the two
physical boundary integrals at rank four or use one representative
boundary integral at rank eight. Define

```
P_f = r_f/[360(4pi)^2].
```

The Euclidean modulus cocycle has sign `+integral_0^1 a4(sigma,L_t)dt`.
The Lorentzian action convention of the accepted spherical anomaly
notes has the opposite sign. Its bulk trace is minus the Dirac heat
coefficient. No phase/eta assertion is needed here.

## Full mixed boundary heat coefficient

Use the general mixed coefficient, not the bag specialization, in
[Vassilevich, (2.6), (3.27)-(3.28), (5.27)-(5.33)](https://arxiv.org/html/hep-th/0306138v3).
Its normal is inward; its R_anan is minus the usual sectional curvature.
With K=tr K_ab, K2=tr K_ab^2 and K3=tr K_ab^3, the traced boundary
coefficient for arbitrary scalar f reduces to

```
a4,partial(f) = P_f integral_partial sqrt(h) {
 f[-3 R_;n -6 Delta_partial K +5 R K
   +4 R_anan K +18 R_anbn K_ab +4 R_abcb K_ac
   +(34 K^3+78 K K2-232 K3)/21]
 +f_;n (36 K^2-54 K2)/7 -6 K f_;nn } .             (1)
```

Repeated tangent indices are summed. The magnetic spin term drops out
of these boundary traces; it remains in the bulk anomaly. Useful
intermediate traces fixing signs and coefficients are

```
tr E=-r_f R/4, tr chi=tr(chi E)=0,
chi_:a=tau_3 K_ab gamma_b,
tr(chi_:a chi_:a)=r_f K2,
tr(chi_:a chi_:b K_ab)=r_f K3,
tr(chi chi_:a Omega_an)=r_f K_ab R_anbn/2.
```

The last identity uses the minus sign in the spin-curvature part of
Vassilevich (3.28). For example, the bare cubic projector terms and
S terms give `(34K^3-300KK2+272K3)/21`; chi derivatives add
`18KK2-24K3`, yielding (1). On a flat round four-ball its cubic
bracket is 44k^3 and the one-Dirac integrated coefficient is 11/180.
At K=0, (1) reduces to the previously retained normal-curvature flux.

## Literal finite-cut path integral

On g_t the principal curvatures are

```
k_time=exp(-t sigma)(a+t s),
k_sphere=exp(-t sigma)t s  (twice).
```

The scalar curvature is

```
R_t=exp(-2t sigma) Rcal_t,
Rcal_t=R0-6t(s'+a s)-6t^2 s^2.                     (2)
```

The four usual sectional curvatures, after removing exp(-2t sigma),
are `-(a'+a^2)-t(s'+a s)`, `-t s'`,
`-t a s-t^2 s^2`, and `Y-t^2 s^2` for normal/time,
normal/sphere, time/sphere, and sphere/sphere respectively. Also

```
f=sigma,
f_;n=-exp(-t sigma)s,
f_;nn=exp(-2t sigma)(s'-t s^2).
```

Each term in (1) has weight minus three. Multiplication by sqrt(h_t)
cancels that factor exactly. The f bracket then has polynomial
coefficients

```
P0=-6a''+44a^3/7+10aY,
P1=-18s''-18a s'-6a's+12a^2 s+10Ys,
P2=14a s^2, P3=14s^3.
```

The remaining smearing bracket is
`18a^2s/7-6as'+t(-66as^2/7-18ss')-36t^2s^3/7`.
Integrating this polynomial proves the full literal boundary part

```
H_partial = 4pi b^2 P_f integral dU_E N Fcal,

Fcal = sigma Acal + Bcal,
Acal = -6a''+44a^3/7+10aY
       -9s''-9a s'-3a's+6a^2s+5Ys
       +(14/3)a s^2+(7/2)s^3,
Bcal = 18a^2s/7-6as'-(33/7)as^2-9ss'-(12/7)s^3.   (3)
```

This is an explicit local functional at every finite cut for arbitrary
smooth static N and sigma. No expansion near infinity or near the
horizon was used. Its Lorentzian literal-path contribution is
`-H_partial`.

## Conversion to the retained bulk representative

It would be incorrect to add `-H_partial` directly to the already
integrated-by-parts bulk Lpoly. There are three conversion fluxes.
First write the bulk anomaly as `c C^2-A Q+d Box R`, with

```
A=11/360, d=1/30-2A/3=7/540,
Q=Euler-(2/3)Box R.
```

At a finite cut the exact identity

```
integral sqrt(g) sigma Box R
 = -(1/12) partial_t integral sqrt(g)R^2
   -integral_partial sqrt(h)(R n_out.sigma-sigma n_out.R)
```

follows from two integrations by parts and the infinitesimal Weyl
variation of sqrt(g)R^2. Thus the local R^2 primitive in the accepted
Wraw differs from the literal bulk path by a boundary term. In the
proper collar above define

```
Jcal = s[R0-3s'-3as-2s^2]
 -sigma[R0'-3s''-3a's-3as'-sR0+4as^2+3s^3].         (4)
```

The required Lorentzian correction relative to Wraw is
`-4pi b^2 P_f integral N [Fcal+(7/6)Jcal]dU`.
Here 7/6=90d for one Dirac; both sides scale with the number of
identical exterior copies.

Second, the retained identity Wraw=Wpoly+flux(B) is
[the spherical anomaly action, (4)](spherical-anomaly-action.md).
Its explicit normal flux is

```
Bflux = d[(R2+2Y)s-sigma R2']
 -2A{sigma(s''+a's+as')-(s'+as)s
                         +(R2-4Y)sigma s/3}.       (5)
```

Third, to use the convenient physical-radius bulk density which has
no sigma Box R2 term, two more integrations by parts supply
`(sigma R2'-s R2)/30` in the same flux units. For one Dirac the
complete Lorentzian local boundary functional accompanying that bulk
representative is consequently

```
W_partial,total = 4pi b^2 P_f integral dU N Ccal,
Ccal=-Fcal-(7/6)Jcal+90Bflux+3(sigma R2'-s R2).       (6)
```

All quantities on the right are explicitly defined in (3)-(5).
Equivalently Ccal=-Gcal, where elementary expansion gives

```
Gcal = sigma[12aa'+44a^3/7+10aY]
       -6a's-(24/7)a^2s-6as'-(96/7)as^2
       -18ss'-(85/21)s^3.                           (7)
```

The associated bulk density, per one exterior and one Dirac, is

```
W_bulk = (1/(4pi)) integral dU dr exp(alpha)/r^2 {
 sigma[Rbar^2/60-Rbar/18+1/15-q^2/3]
 +Rbar T/30-T^2/10+11TF/90-11F^2/180-F/18 },
T=r(F'+F alpha'), Rbar=b^2 Rhat_2.                  (8)
```

Primes in (8) denote physical r derivatives. The last conversion flux
in those coordinates is
`[sigma exp(alpha)F Rbar'-exp(alpha)F Rbar/r]/(120pi)`.
Equations (6) and (8), rather than (3) and (8), form the same complete
local Weyl functional as the literal bulk-plus-boundary path. They
refer to a finite cut, with all conversion fluxes retained.

## Off-shell boundary variation and the charge distinction

The unexpanded form (6) or the polynomial (7) determines the complete
static off-shell variation. In proper collar gauge put j=delta N/N,
k=delta sigma. For any polynomial Ccal(a,a',a'',s,s',s'',sigma),

```
delta integral N Ccal
 = integral N[ Ccal j+Ccal_a j'+Ccal_a' j''+Ccal_a'' j'''
              +Ccal_sigma k+Ccal_s k'
                         +Ccal_s' k''+Ccal_s'' k''' ].           (9)
```

The derivatives on Ccal subscripts are partial derivatives with
respect to the indicated independent jet, not extra x derivatives.
For (7), the absent higher jets simply have zero coefficients.
To include the radial metric variation, write ghat_xx=M^2 before
varying and let h=delta M/M, with primes now proper derivatives.
The chain rule used in (9) is replaced by

```
delta a=j'-a h, delta s=k'-s h,
delta a'=j''-2a'h-a h',
delta s'=k''-2s'h-s h',
delta a''=j'''-3a''h-3a'h'-a h'',
delta s''=k'''-3s''h-3s'h'-s h''.
```

This gives all boundary coefficients without choosing which normal
jets an ensemble fixes. Those coefficients must be combined with the
bulk variation of (8). A Dirichlet charge cannot be read by discarding
the j', k' or h' terms. The complete action, boundary prescription and
remaining Weyl-invariant functional determine their cancellation or
the additional data to be held fixed.

For a constant time rescaling j, keeping all other jets fixed, the
partial action variation is proportional to its value. That restricted
variation is not a derivation of the physical boundary energy. It does
not vary the full state, the regular HH period, or the other action
pieces. In particular an RN action evaluated first and differentiated
only in its mass parameter is not a substitute for (9).

## Exact RN outer action tail

This paragraph concerns the value of the local action on the fixed
classical RN reference, not its off-shell charge. Put

```
u=b/r, B=X+u^2-2a_m u^3+h_m u^4,
N=sqrt(B), sigma=-log(u),
partial_x=-sqrt(B)partial_u/b.
```

The exact jets needed in (3)-(7) are

```
a=-B_u/(2b sqrt(B)), s=sqrt(B)/(bu),
a'=[B_uu/2-B_u^2/(4B)]/b^2,
a''=-sqrt(B)[B_uuu-B_u B_uu/B+B_u^3/(2B^2)]/(2b^3),
s'=[B/u^2-B_u/(2u)]/b^2,
s''=sqrt(B)[2B/u^3-3B_u/(2u^2)+B_uu/(2u)]/b^3.     (10)
```

Expansion of the literal boundary functional gives

```
N b^2 Fcal = (1/b){
 sigma[-29X^2/(2u^3)-14X/(3u)]
 -75X^2/(7u^3)-12X/(7u)-(114/7)X a_m+o(1)}.        (11)
```

Its mass-dependent sigma term cancels: there is no a_m log(u)
divergence in this boundary-only value. After all conversion fluxes,
the total boundary value accompanying (8) is instead

```
N b^2 Gcal = (1/b){
 -463X^2/(21u^3)-260X/(21u)+(232/21)X a_m+o(1)}.    (12)
```

For a direct algebraic check, the nonsigma part of N b^2 Gcal is
exactly

```
(1/b)[-3B B_uu/u-(6/7)B_u^2/u
                   +(132/7)B B_u/u^2-(463/21)B^2/u^3].
```

The remaining sigma part is O(u log(u)). This proves (12) without
exchanging an infinite expansion with a heat integral.

Thus the finite value of (6), after deleting just its displayed power
divergences, is `-4pi P_f (232/21)m/L4^2` per unit U. This specifies
an algebraic finite-part operation on this component; it is not a
complete holographic renormalization or a boundary energy theorem.
The representative, sign and conversion fluxes matter even for this
limited finite-part statement.

There is no inference here that the full remaining stress has an
integrable mass tail on the unchanged bare-L4 exterior. Bulk W4 can
have a mass-dependent logarithmic tail after its vacuum term is
subtracted, and the boundary-layer/state contribution of the remaining
functional must be included in the same scheme. The companion
[massive boundary-layer calculation](stationary-massive-boundary-layer.md)
derives the leading cancellation on the fixed RN reference with a
qualitative integrable outer remainder; it does not give its global or
temperature-uniform differentiated bound. An asymptotically AdS
metric or a local anomaly alone does not prove cancellation. The
conditional outer-vacuum centering distinction in the conformal-domain
note remains in force.

As a relevant comparison, [Taylor--Breen, (65)](https://arxiv.org/html/2007.09405v2)
retain an explicit log(r/r_+) term in the radial mass correction from
the scalar Page approximation on Schwarzschild--AdS. That is an
approximate scalar calculation, not evidence for the exact Dirac tail.
It supplies no general no-log theorem for the present Dirac problem.

Remaining to-do list: combine this local functional's off-shell
variation with the bulk and Weyl-invariant stationary functional under
the physical boundary prescription; control the subtracted boundary
layer and full stress remainder; solve full exterior/state matching
and complete the S074-S079 gates.
