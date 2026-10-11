# Propagating photon polarization with unequal returns

2026-10-11; S053. Ordinary argument and reproducible finite checks; no Lean
certification. This closes a **quadratic fermion-influence comparator**,
not the interacting clock stress response or physical wormhole matching.
It follows the [contact obstruction](angular-current-renormalization.md).

## Model and result

Take q free massless Dirac channels in their neutral Neveu-Schwarz vacua.
The common throat coordinate is s in [-T/2,T/2]; channel j completes a
circle of length L_j=T+D_j, with any finite D_j>0. Right and left movers
are independent. Fix 0<T_min<=T<=T_max. Photon vertices are the real
Hermitian lowest-Landau-level matrices T_A, A=(l,a), 1<=l<q, with

    Tr(T_A T_B)=tau_l delta_AB,
    tau_l=q beta_l/(4 pi),
    beta_l=(q-1)! q!/[(q-1-l)! (q+l)!],  kappa_l=l(l+1).

Write P=diag(tau_l), M=diag(kappa_l), repeated over spatial and two
potential components. Neither unequal returns nor the noncommuting
vertices are diagonalized into independent scalar fields. All operator
norms below use the direct sum of throat L2 spaces.

Let chi(w) be the exact free retarded current commutator, with convention
-i theta(t)[J(t),J(0)], and let

    Sigma(w)=chi(w)+S,  S=(P/pi) tensor diag(0,1).

S fixes the local vector-current contact in the potential components
(A_0,A_1), coupled to (rho,j). It follows from local covariant point
splitting; its normalization is also checked by the Ward identity below.
Returns may mix angular modes in chi. This is an explicitly specified
one-loop projected model, not a claim of a complete four-dimensional Ward
identity at the mouths or at higher orders.

For eta=Im(w)>0 define

    m_T=max(2/pi^2, (2/pi)(1/(2 pi)+1/T_min)),
    A_T=4 pi m_T (1+T_max^2/4),
    B_J(w)=4 A_T (1+|w|/eta)+1/pi.

Then, with no fermion-frequency cutoff,

    ||P^(-1/2) Sigma(w) P^(-1/2)|| <= B_J(w).                 (1)

For fixed-radius, fixed-metric Maxwell electric multipoles on a global
AdS2 throat cut at |x s|=a<pi/2, choose Dirichlet boundary conditions for
the dual electric scalar. These are comparator boundaries, not physical
mouth conditions. Let D_0(w) be their propagator in potential variables.
If x>=x_min>0 and |w|<=W<2 x_min, then

    ||M^(1/2) D_0(w) M^(1/2)|| <= C_ph,
    C_ph=2/(1-W^2/(4 x_min^2)).                              (2)

The quadratic resummation is D_g=(D_0^(-1)-g^2 Sigma)^(-1), understood
equivalently through its bounded resolvent identity. Put

    epsilon=g^2 C_ph B_J max_l(tau_l/kappa_l),
    max_l(tau_l/kappa_l)=q beta_1/(8 pi)<=q/(8 pi).

If epsilon<1, D_g exists and

    ||M^(1/2)(D_g-D_0)M^(1/2)|| <= C_ph epsilon/(1-epsilon). (3)

For Sigma_eff=Sigma+g^2 Sigma D_g Sigma, the corresponding P-normalized
correction is at most B_J epsilon/(1-epsilon). This permits arbitrary
angular mixing at the free two-point level. In the working family
g^2=q^-2, eta comparable to delta=q^-1/16, epsilon=O(q^-15/16).
The constants are independent of the longest return, number of angular
modes and a approaching pi/2, provided T remains in the stated interval.
This is a uniform estimate at a shrinking positive imaginary resolution;
it is not a real-axis or global stability theorem.

## Proof of the current bound

A right-moving particle in channel i and hole in channel j have energy

    E_nm=2 pi [(n+1/2)/L_i+(m+1/2)/L_j],  n,m>=0.

Their current matrix element on the throat is T_A,ij exp(i E_nm s)
divided by sqrt(L_i L_j), up to the harmless convention for Fourier sign.
This follows directly by expanding the fermions and applying the CAR to
the filled sea. Define the positive measure

    mu_ij=(1/(L_i L_j)) sum_nm delta_(E_nm),  nu_ij=dmu_ij/E.

The number N(R)=mu_ij((0,R]) satisfies N(R)<=R^2/pi^2: each occupied
half-odd grid below R has at most L_i R/pi points (if any point exists).
Integration by parts therefore gives nu_ij((0,1])<=2/pi^2. For integer
k>=1, counting n up to k+1 and then m in an interval of width one gives

    mu_ij([k,k+1)) <= (k+1)/pi [1/(2 pi)+1/L_j],
    nu_ij([k,k+1)) <= (2/pi)[1/(2 pi)+1/T_min].

Endpoint conventions can be chosen consistently without double counting;
the same upper bounds hold for closed unit bins. Thus every nonnegative
integer unit bin of nu has mass at most m_T. For f supported in the
throat, use the unnormalized Fourier transform F(E)=integral e^(-iEs)f(s)ds.
On each unit interval the elementary one-dimensional Sobolev inequality
gives sup|F|^2<=2 integral(|F|^2+|F'|^2). Plancherel, F'=FT(-isf), gives

    integral |F(E)|^2 dnu_ij(E) <= A_T ||f||_2^2.             (4)

This requires no spatial derivative of f; even a sharp throat-supported
L2 source is allowed. Reflection E to -E has the same estimate. Since
E/|w-E| and E/|w+E| are at most 1+|w|/eta, Cauchy-Schwarz in (4) bounds
each signed Lehmann branch. The two branches cost a factor two.

For vector sources f_A,0 and f_A,1 use the chiral combinations f_A,0+f_A,1
and f_A,0-f_A,1. For each chirality sum the pair functions
F_ij(s)=sum_A T_A,ij f_A(s). The matrix trace identity gives exactly

    sum_ij ||F_ij||_2^2=sum_A tau_A ||f_A||_2^2.

Transposition i,j in the negative branch only relabels the sum, since
mu_ij=mu_ji. The two chiral norms add to twice the vector norm. Hence
the commutator contributes 4 A_T(1+|w|/eta); S contributes 1/pi.
This proves (1). The same domination proves convergence of all smeared
quadratic forms and strong operator convergence of frequency truncations.
Local uniform bounds and the Lehmann representation give holomorphy in
the upper half-plane. No uniform operator-norm ultraviolet tail is claimed.

The bound concerns susceptibility, whose energy denominator matters. It
does not make the contact interaction applied to the free vacuum have a
finite norm, and does not contradict the preceding contact obstruction.

## Propagator and resummation

The Euclidean Maxwell action has alpha(s)=cos^2(xs)/x^2 and

    C_E(k)=kappa I+B(k)^* alpha B(k),
    B(k)=(-partial_s, i k),
    H_l=-partial_s^2+kappa/alpha,
    C_E(k)^(-1)=kappa^(-1)[I-B(k)^*(H_l+k^2)^(-1)B(k)].      (5)

Derivatives and Dirichlet dual-scalar boundary conditions are interpreted
as closed forms; (5) also follows by solving for the electric scalar.
Continue k=-iw: B=(-partial_s,w), B_sharp=(partial_s,-w).
For physical potential components, C_M=U C_E(-iw) U with U=diag(1,-i),
and D_0=U^(-1) C_E(-iw)^(-1) U^(-1). These unitary factors preserve norms;
the continuation uses a transpose-type analytic companion, not an adjoint
at complex w. This fixes the relative sign of the spatial component.

The full global-throat Friedrichs operator obeys H_l>=(l+1)^2 x^2,
by the factorization of -partial_sigma^2+l(l+1)sec^2(sigma).
Dirichlet restriction raises this lower bound. Set lambda_*=4 x_min^2.
Both H_l^(-1/2)B and B_sharp H_l^(-1/2) have norm at most
sqrt(1+W^2/lambda_*); H_l^(1/2)(H_l-w^2)^(-1)H_l^(1/2) has norm at most
(1-W^2/lambda_*)^-1. The triangle inequality in (5) proves (2).
It applies uniformly to all l>=1; the singlet and metric perturbations
are separate sectors.

Conjugate the resolvent identity by M^(1/2). From (1),
||M^(-1/2) Sigma M^(-1/2)||<=B_J max(tau/kappa).
The Neumann series therefore converges when epsilon<1 and gives (3).
For the current correction conjugate instead by P^(-1/2) and use
||P^(1/2)D_g P^(1/2)||<=max(tau/kappa) C_ph/(1-epsilon).
No simultaneous diagonalization of Sigma and M was used.

## Known commuting baseline and independent checks

For a common translation-invariant circle, spatial Fourier momentum p,
one vertex with Tr(T^2)=tau, z=w^2-p^2 and s_g=g^2 tau/pi, the exact free
commutator and the contact give

    chi=tau/[pi z] [[p^2,pw],[pw,p^2]],
    Sigma=tau/[pi z] (p,w)(p,w)^t,
    Sigma (w,-p)^t=0,
    C_M=[[kappa+alpha p^2,alpha wp],[alpha wp,alpha w^2-kappa]],
    det(C_M-g^2 Sigma)=kappa[alpha z-kappa-s_g].              (6)

The dressed current has its gapless residue multiplied by
kappa/(kappa+s_g), and a massive pole at z=(kappa+s_g)/alpha.
Equivalently the boson propagator is

    (Q^2+kappa/alpha)/[Q^2(Q^2+(kappa+s_g)/alpha)]
      = kappa/(kappa+s_g)/Q^2
        + s_g/(kappa+s_g)/(Q^2+(kappa+s_g)/alpha).

A separate positive two-scalar Gaussian with canonical kinetic terms and
mass potential (kappa/alpha)(y-h phi)^2, h^2=s_g/kappa, gives this phi
propagator by orthogonal rotation. This independently checks the mass
shift, positive residues and sign in (6). These are known Sommerfield
results: [Georgi and Warner, v2, sections 1, 2 and 4](https://arxiv.org/html/1907.12705v2).
Their exact multi-vector solution assumes commuting fermion couplings.
Noncommuting LLL vertices and unequal returns prevent importing it as a
solution of the present full theory.

The [check script](../checks/propagating_photon_polarization.py) reproduces
(6) symbolically, checks the spectral pair formula against explicit
finite-sea CAR signs and energies, checks counting and the Fourier bound, compares
(5) with direct compatible-matrix inversion, and retains exact Fourier
integrals of piecewise-constant sources in the fermion pair sum. It tests
equal and unequal returns with every angular vertex for q=2,3, verifies
the weighted bounds, and refines both frequency cutoff and spatial mesh
for a fixed smooth photon probe. Finite checks support the stated ordinary
proof; they do not certify continuum physical matching.
Unequal returns produce angular off-diagonal mixing in the q=3 example.
For q=2 the chosen real harmonic basis stays diagonal: there is only one
diagonal traceless generator and the x/y cross term cancels between pairs.
The retained first numerical failure incorrectly required mixing at q=2;
the corrected assertion checks this structural exception explicitly.

## Consequence and next decision

The entire projected free polarization, including angular mixing, is
small after photon propagation at damping resolution in the working
family. The next obstruction is the stress observable: two metric
derivatives of the renormalized photon/fermion influence functional are
not bounded by (1)-(3). Those derivatives, the photon vacuum trace,
higher current cumulants, and a controlled interacting state remain open.

There is also a concrete matching caution. In flat constant-field QED,
[Gorbar, Hashimoto and Miransky, v1, equations 16-26](https://arxiv.org/html/hep-th/0501135v1)
find transverse LLL two-point polarization but a three-point Ward defect
repaired by the collective higher-Landau-level contribution. Their setup
has a fermion mass, many flavors and a noncompact transverse plane. It
does not prove such a defect or its size on this spherical throat; it
does rule out assuming that two-point control automatically closes all
four-dimensional gauge matching.

Continue the physical long-return damping primary and finite-duration
energy-transfer fallback. The commuting baseline is established prior
work; novelty of the physical result is unaccepted. This comparator does
not replace the target with a current-only graph paper. No construction,
physical damping, finite-disturbance survival or shortcut is claimed.

Remaining to-do list: control renormalized clock stress and higher current
terms, close metric/photon mouth and higher-mode matching plus the state,
finish novelty and the consequential result, then independently verify,
GCP-certify suitable endpoints and prepare the reproducible manuscript.
