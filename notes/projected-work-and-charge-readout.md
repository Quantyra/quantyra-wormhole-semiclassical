# Separating finite state contributions in cumulative work and endpoint charge

This provisional methods note isolates part of a finite quantum readout when two state contributions were computed in different bases. Both increments and the observables must be expressed in the same physical Hilbert space at each time.

## A common readout

Consider a covariance increment
$$
\Delta C=(I+B)C(I+B)^\dagger-C .
$$
Suppose the range of $B$ lies in a known subspace $S$. For an orthogonal projector $N$ onto a subspace of $S^\perp$, we have $NB=0$. Expanding gives
$$
N\Delta C N=N(BC+CB^\dagger+BCB^\dagger)N=0 .
$$
The calculation permits off-diagonal covariance and a nonunitary factor $I+B$. A change of output basis transports the support and projector together.

Let $Q$ project onto $S_1^\perp$, and let the columns of $V$ span the second support $S_2$ in those coordinates. Put $W=QV$ and $N=Q-WW^+$, where $W^+$ is the Moore-Penrose inverse. The projector $WW^+$ has range inside the range of $Q$. A vector in the range of $Q$ is orthogonal to $QV$ exactly when it is orthogonal to $V$. Hence $N$ projects onto $S_1^\perp\cap S_2^\perp$. It annihilates both covariance increments under two-sided compression, even if their starting covariances differ. In dimension 650, supports of dimensions at most 550 and 4 leave at least 96 directions. A numerical evaluation still needs subspace accuracy and evidence that the selected signal matters.

For a measurable family of such projectors, split the work operator $K(t)$ into $N(t)K(t)N(t)$ and its exact remainder. Make the same split for the charge operator at the initial and final times. The remainders retain the cross terms. Trace cyclicity shows that both earlier increments contribute zero to compressed work and compressed charge. Possible overlap in the complementary readout remains to be assessed.

## Signed variation

Let $U_0(t)$ map the initial interaction coordinates to the physical coordinates, and write $C_*$ for the full initial covariance in the interaction picture. For a Hermitian interaction increment $dA$ in that picture, define $M(t)=\int_0^t dA(s)$. The coefficient linear in the perturbation is $X_1(t)=-i[M(t),C_*]$, with zero initial increment. Selecting an off-diagonal pair requires retaining both conjugate entries.

Put $\widehat K(t)=U_0(t)^\dagger N(t)K(t)N(t)U_0(t)$ and define $\widehat L(T)$ by the same conjugation for the endpoint charge. For a real weight $u$, the degree-one contribution to cumulative work minus final charge is
$$
B(s)=\int_s^T\widehat K(t)\,dt-\widehat L(T),
\qquad
F_1=-iu\int_0^T\operatorname{Tr}\!\left([C_*,B(s)]\,dA(s)\right).
$$
Substitute $X_1$ under each physical trace, interchange the order of integration, and cycle the factors. Bounded readouts over a finite interval and an interaction with finite total variation suffice. The projector enters only as a readout, so its derivative is unnecessary.

When $U_0=RD$ uses an orthonormal real eigenframe $R$ and diagonal phase propagation $D$, the interaction includes the connection $-iD^\dagger R^\dagger dR D$ and the conjugated remainder of the Hamiltonian. A jump in the chosen frame can be represented by a specified continuous unitary path at fixed physical time. Work and dynamical phase are frozen along that path. The degree-one coefficient depends on the chosen gauge and interpolation. The full ordered connection evolution implements the change of coordinates.

Application to an original physical source requires signed temporal integrals, controlled state and readout inputs, all complementary observables and the remaining orders. The physical response, its Ward identity and its continuum error remain open.
