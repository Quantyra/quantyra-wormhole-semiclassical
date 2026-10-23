# Coincident Dirac action assembly for the physical correction calculation

2026-10-23, E003/S169. Ordinary argument before execution.
The full stable semiclassical wormhole objective remains incomplete.

The closed audit selects the full physical H e=-r1 solve. Inspection of
the current executables finds radial Green functions and a separated
kernel evaluator, but no assembled renormalized operator. This implements
the missing coincident finite-mode action, including its tadpole and the
variation of the actual paired boundary data. It is reusable operator
assembly, not a second estimate of a separated kernel value. The finite
mode matrices below are DEVELOPMENT INPUTS for the full solve. They are
not the complete quantum Hessian, a truncation certificate or a correction.

## 1. The same physical domain and a complete finite-mode variation

Use the actual parameters and polynomial Q of
[the finite-seed Green evaluator](current-seed-radial-green-enclosure.md).
Let rho increase from the bolt toward the seam, dz/d rho=-sqrt(Q),
M_n^2=n(n+q)/C, d_n=2(q+2n), p=gamma(j+1/2) and W=p/sqrt(Q).
The canonical transformed Dirac expression and variable mass are

    D=[[0,partial_rho+W],[partial_rho-W,0]],
    m=M exp(b), Gamma_M[b]=-Tr log(D+m).

Here Gamma_M is its regulated mass-vertex part. The separate two-dimensional
Weyl allocation and the common four-dimensional subtraction/restoration
must still be added. Calling this regulated part the renormalized action
would be wrong. The angular weights and both signs of thermal frequency
are retained explicitly below; there is no extra spin or mouth factor.

For a left regular solution write its spinor as (-x-y,x-y). The homogeneous
equation is exactly

    x'=m x+W y, y'=W x-m y,
    U=y/x, U'=W(1-U^2)-2mU.                            (1)

The right solution, propagated inward, obeys the same equation for V.
Even paired transmission requires V_seam=U_seam. In particular its
variation is also transmitted; fixing the right datum during a variation
would remove a physical reflected response. A compact direction supported
between z=1/2 and z=7/8 does not change the inner regular solution at
z=15/16. Its unknown unperturbed value is enclosed by 0<=U_cut<=1 and
its first variation there is zero. This is an enclosure of the physical
bolt solution, not a boundary wall.

For a nonnegative direction u, define R=-D_u U and S=-D_u V. Differentiating
(1), on the SAME background and paired domain, gives

    R'=-(2M+2WU)R+2M u U,
    Sdot=-(2M+2WV)S+2M u V,
    R_cut=0, S_seam=R_seam.                            (2)

All four functions are nonnegative, U,V<=1 and R,S<=sup(u). The last
bound follows because the forcing divided by the damping is <=u U<=u.
The source trace at coincident points, with the first-order matrix jump
and spin trace taken before estimating, is

    T=(1-UV)/(1+UV), tr S_M(rho,rho)=T.

Thus its source and complete directional derivative, per radial thermal
mode, are

    F_M=-M T,
    D_u F_M=-M u T-2M(VR+US)/(1+UV)^2.                 (3)

The FIRST term is the exponential-mass tadpole. The second is the full
bubble, including the response propagated from the paired seam. Both
are essential. Equations (1)-(3) are valid for arbitrary smooth background
Q>0 in this band, not just at a constant mass potential W.

For two nonnegative compact test functions u_i,u_j, the contribution of
the pair of frequencies +omega,-omega to the unrenormalized Hessian is

    H_ij=-2M integral [u_i u_j T
                  +2u_i(VR_j+US_j)/(1+UV)^2] d rho.   (4)

Multiplication by d_n accounts for one angular level. The normalized
fixed-action version multiplies by -2G/B^2 as elsewhere. Neither this
integrated action matrix nor a single-mode value is a pointwise source
norm. A regulated matrix may have negative entries; its sign does not
determine that of the full, commonly subtracted quantum form.

## 2. A second calculation of the complete contact

For this even pair, the right adjoint solution at the seam is proportional
to the left solution (x,y). The half-double characteristic determinant
therefore gives, up to a direction-independent boundary factor,

    Gamma_(one frequency)=-log x_seam
                          -log(1+U_seam^2)/2.         (5)

This also follows without a determinant theorem: differentiate the right
side using the adjoint of the two-by-two transfer equation. Its derivative
is -integral m u (1-UV)/(1+UV) d rho, the exact trace variation (3).
The integration constant is independent of all compact mass variations,
so its second derivative vanishes. No new self-adjoint boundary condition
or finite counterterm is introduced by (5).

For constant W,m on a cell of length L let K=sqrt(m^2+W^2), T0=tanh(KL).
The exact normalized transfer and its logarithmic amplitude are

    Unew=[W T0+(K-mT0)U]/[K+mT0+W T0 U],
    delta log x=KL+log([(1+E)+(m+WU)(1-E)/K]/2),
    E=exp(-2KL).                                      (6)

Differentiate the product of (6) twice at m=M exp(t u) and then (5).
The executable uses independent order-two arithmetic, including the
second derivative m''=M u^2. This supplies an independent full tadpole
and boundary check for the piecewise constant coefficient approximation.
It does not turn that approximation into the continuum operator. The
continuum enclosure is obtained separately from (2)-(4).

## 3. Validated continuum integration

Use the inherited logarithmic mesh from z=15/16 to the actual seam z=0.
Each cell encloses W and its true proper length. The constant-flow
comparison from (1) gives bounds on U at both endpoints and throughout
the cell. For R in (2), let a_min,a_max enclose 2M+2WU and let
f_min,f_max enclose 2M u U. Scalar comparison gives endpoint bounds from

    Rlo'= -a_max Rlo+f_min,
    Rhi'= -a_min Rhi+f_max.                            (7)

These constant equations have exact exponential solutions. Enclose their
whole cell by their initial and final values, since each is monotone
toward its equilibrium. This also covers the enclosed true cell length.
Nonnegativity makes the coefficient ordering in (7) valid. The same
construction applies inward to V,S, starting with the paired endpoint
values and sensitivities. It is not interval evaluation of a constant
coefficient flow pretending the coefficients are constant in the geometry.

Insert the resulting whole-cell bounds in (4), integrate using the
enclosed proper length and sum. Every point of the continuous support
is included. Mass, temperature and geometry intervals are retained.
The test functions are

    a=1/2, b=7/8, w=(z-a)/(b-a),
    u1=[4w(1-w)]^6, u2=w u1, a<=z<=b,
    u1=u2=0 elsewhere.

They are C5 across their support endpoints and 0<=u_i<=1. They are
assembly probes, not asserted values of h1 or of the unknown correction.
The finite thermal labels are selected by j=floor(p_target/gamma_cl)
and the actual p=gamma(j+1/2) is retained. The selected cases are
(n,p_target)=(1,.25),(1,1),(2,1), on N=256,512,1024 logarithmic
subdivisions per factor two. None stands for an omitted angular or
thermal sum. This execution is a matrix assembly diagnostic, not another
separated-kernel observable or a claimed full-solution feasibility test.

Before execution require: negative finite-regulator diagonal entries;
overlap of independently assembled H12 and H21; successive continuum
enclosures overlap and narrow; and all entries of the independent
piecewise-transfer matrix lie within the final continuum intervals.
Require final relative widths below .15. No acceptance threshold on
the physical solution is relaxed. Preserve any failed run and its source.

## 4. The actual lowest thermal mode through the bolt

After the first successful middle-band execution, extend the SAME
primitive to the physical bolt before it is used in a global operator.
Preserve the first snapshot/run and replay the whole expanded package.
This is an actual finite-temperature endpoint calculation, not a frozen
zero-temperature or flat comparison. Use j=0, omega=1/2, n=1,2.

Direct evaluation of Q near its tiny positive-temperature root loses
the correlation of the root and thermal coefficient. Avoid that by an
exact factorization. Put v=1+h, t=r/B, A=3beta+u+epsilon and

    K(v)=(v+1)[u+3beta(v^2+1)]/(2v^3)
                         -epsilon(v-3)/(2v^5),
    H1(t,v)=beta(t+2v)+A/(t v^2)
        -2epsilon(2t+v)/(t^2 v^3)
        +epsilon(3t^2+2tv+v^2)/(t^3 v^4).

The same prescribed-temperature root has gamma=h K(v). Subtracting
the horizon value of t F_ext(t)/C, then factoring t-v, gives EXACTLY

    F/C=(t-v)[2h K(v)v+(t-v)H1(t,v)]/t.               (8)

Use y=cosh(chi)-1, t=1+h cosh(chi) and Bchi=2K(v)v+yH1. Then

    Q=h^2 y Bchi/t^3,
    d rho/d chi=sqrt((y+2)/(t Bchi)),
    W=omega K(v) t^(3/2)/sqrt(y Bchi).                (9)

The h factor cancels analytically in the last two formulas. This is
preservation of an exact geometric identity, not a change of physical
normalization or a relaxation of the bolt target. Substitute the current
h enclosure; every actual fixed-temperature root is included. The
global interval on 0<=chi<=3 verifies 1/2<d rho/d chi<2. No point sample
is used to decide the positivity of the factored geometry.

On this collar choose smooth radial probes

    v1=(1-chi^2/4)^6, v2=(chi^2/4)v1, 0<=chi<=2,
    v1=v2=0 for chi>=2.

They are smooth at the bolt and C5 at the outer edge, with 0<=v_i<=1.
Start at chi_min=1e-6 with U,R_i in [0,1]; at chi=3 start the inward
V,S_i in [0,1]. These are proved universal enclosures of the regular
HH solution and its compact-direction sensitivities. No guessed wall
value, missing initial derivative or decoupled seam is substituted.
The broad outer data include the actual paired response. Damping through
the retained buffer narrows them before the probe's support is reached.

For the omitted tip 0<=chi<=chi_min, (3), U,V,R,S,v_i<=1 and
d rho/d chi<2 give each frequency-pair matrix entry an additional
interval [-20 M chi_min,0]. Thus the true bolt is included, rather than
discarded as a measure-zero endpoint. Use 20N logarithmic cells from
chi_min to 3 at N=64,128,256. Whole-cell comparison (7) remains valid
with the lengths and W from (9). Require reciprocity, overlapping and
narrowing enclosures, and final relative widths below .15, as before.

This additionally tests the most delicate j=0 regular-bolt channel on
the current finite geometry. It does not replace its full mode sum,
fix the unknown metric at the bolt, or quantify the full Hessian inverse.

## 5. What remains before this can enter H e=-r1

The uncut frequency sum in (3)-(4) is divergent. The angular sum is also
not a physical Hessian until the common a0,a1,a2 subtraction, the Weyl
allocation and finite restorations are included. The existing a0/a2
budgets are not permission to set them to zero. The current evaluator
must expose its regulator and refuse to label a finite matrix a complete
renormalized one. The original scale is fixed. Source and Hessian must
use the same subtraction and actual-domain propagation.

Next incorporate this coincident functional into that common subtraction
and full angular/thermal summation, alongside the local scalar/mixed
rows and the physical mass/log boundary realization. The desired full
calculation remains the actual coupled correction and its weighted
residual/admission, not success of these probes. The compact support used
to validate the primitive does not settle bolt/outer source expansions,
an entire nonlinear ball or the retarded gravitational domain.

The previous audit remains closed. Its <1.667e-9 bolt difference,
<6.67e-10 exterior return, Y0/admission and nonlinear targets are unchanged.
No stability, physical state, validity or wormhole existence claim follows
from the finite matrices. The technical advance is an executable, complete
coincident finite-mode variation that can be used in the full operator,
including contacts and the differentiated paired boundary map.

Remaining to-do list: finish the commonly subtracted full operator and
source assembly and execute the physical coupled correction with weighted
error/admission; nonlinear compatibility; physical matching, causal
realization, semiclassical validity, all-sector stability, charge/return,
opening, transport, preparation, separate creation and the full goal.
