"""S072: remove the common charged boson before computing neutral noise.

The exact reference is the unfolded free-fermion strip. Quantum finite-band
checks verify the source-response structure, not a 4D switching stress.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import expm, null_space
from scipy.integrate import quad, solve_ivp

ROOT=Path(__file__).resolve().parents[1]


def L_euclidean(t):
    return 2*mp.log(mp.coth(abs(t)/2))


def wick_cases():
    mp.mp.dps=75
    rows=[]
    for n in (1,2,3,5):
        for length in (2,3,4,6):
            times=[mp.mpf('.17')*i+mp.mpf('.037')*i*i for i in range(length)]
            signs=[1 if (i*i+length)%3 else -1 for i in range(length)]
            species=[(i*i+i+length)%n for i in range(length)]
            determinant=mp.mpc(1)
            for species_id in range(n):
                ids=[i for i in range(length) if species[i]==species_id]
                if not ids:continue
                matrix=mp.matrix(len(ids))
                for aa,i in enumerate(ids):
                    creation=0 if signs[i]==1 else mp.pi
                    for bb,j in enumerate(ids):
                        annihilation=mp.pi if signs[j]==1 else 0
                        matrix[aa,bb]=-1j/(4*mp.pi*mp.sin((annihilation-creation+1j*(times[j]-times[i]))/2))
                determinant*=mp.det(matrix)*mp.fprod(1j*signs[i] for i in ids)
            neutral_exponent=mp.mpf(0);charged_exponent=mp.mpf(0)
            for i in range(length):
                for j in range(i):
                    term=-signs[i]*signs[j]*L_euclidean(times[i]-times[j])
                    neutral_exponent+=term*((1 if species[i]==species[j] else 0)-mp.mpf(1)/n)
                    charged_exponent+=term/n
            reconstructed=(1/(4*mp.pi))**length*mp.exp(neutral_exponent+charged_exponent)
            error=abs(determinant-reconstructed)/max(abs(reconstructed),mp.mpf('1e-50'))
            assert error<mp.mpf('1e-60')
            rows.append(dict(n=n,insertions=length,relative_wick_error=str(error)))
    return rows


def power_coefficients(a,limit):
    # F_a(z)=((1+z)/(1-z))^(2a), (1-z^2)F_a'=4a F_a.
    coeff=[mp.mpf(1),4*a]
    for degree in range(2,limit+1):
        coeff.append((4*a*coeff[-1]+(degree-2)*coeff[-2])/degree)
    return coeff


def spectral_cases():
    rows=[];mp.mp.dps=75
    limit=240;z=mp.exp(-mp.mpf('.7')-mp.mpf('.31')*1j)
    L=2*(mp.log(1+z)-mp.log(1-z));b=1/(4*mp.pi)
    for n_int in (1,2,3,4,10,100,2000**4):
        n=mp.mpf(n_int);a=1-1/n;c=1/n
        fa=power_coefficients(a,limit);fc=power_coefficients(c,limit)
        xx=[];yy=[]
        for E in range(limit+1):
            vx=4*n*b*b*(fa[E]+(n-1)*fc[E]) if E>0 and E%2==0 else mp.mpf(0)
            vy=4*n*b*b*(fa[E]-(n-1)*fc[E]) if E%2 else mp.mpf(0)
            assert vx>=-mp.mpf('1e-58') and vy>=-mp.mpf('1e-58')
            if E>0:
                assert vx<=n/mp.pi**2*E+mp.mpf('1e-58')
                assert vy<=n/mp.pi**2*(E-mp.mpf(1)/E)+mp.mpf('1e-58')
            xx.append(vx);yy.append(vy)
        assert abs(yy[1])<mp.mpf('1e-58')
        if n_int<=2:assert max(abs(v) for v in yy)<mp.mpf('1e-58')
        closed_x=4*n*b*b*(mp.cosh(a*L)+(n-1)*mp.cosh(c*L)-n)
        closed_y=4*n*b*b*(mp.sinh(a*L)-(n-1)*mp.sinh(c*L))
        series_x=mp.fsum(xx[E]*z**E for E in range(limit+1))
        series_y=mp.fsum(yy[E]*z**E for E in range(limit+1))
        error=max(abs(series_x-closed_x),abs(series_y-closed_y))/max(1,n)
        assert error<mp.mpf('1e-55')
        rows.append(dict(n=n_int,checked_energies=limit,series_error_per_channel=str(error),
                         first_odd_weights_normalized=[float(yy[E]*mp.pi**2/n) for E in (1,3,5,7,9)],
                         first_even_weights_normalized=[float(xx[E]*mp.pi**2/n) for E in (2,4,6,8,10)]))
    return rows


def gaussian_matrix_checks():
    # Direct oscillator matrix products check non-time-ordered complex
    # phases, including the common-charge removal, independently of Wick.
    levels=35;ann=np.diag(np.sqrt(np.arange(1,levels)),1)
    rows=[];variance=.6;omega=1.3;b=1/(4*np.pi)
    for n in (2,3,5):
        basis=null_space(np.ones((1,n)))
        P=basis@basis.T
        assert np.max(abs(P-(np.eye(n)-np.ones((n,n))/n)))<1e-14
        times=[.1,.7,.3,1.1];indices=[0,1,0,n-1];signs=[1,-1,-1,1]
        result=complex(b**4*np.exp(sum(variance*P[i,i]/2 for i in indices)))
        for axis in range(n-1):
            product=np.eye(levels,dtype=complex)
            for t,index,sign in zip(times,indices,signs):
                coordinate=np.sqrt(variance)*basis[index,axis]*(ann*np.exp(-1j*omega*t)+ann.T*np.exp(1j*omega*t))
                product=product@expm(1j*sign*coordinate)
            result*=product[0,0]
        exponent=0j
        for i in range(4):
            for j in range(i+1,4):
                exponent-=signs[i]*signs[j]*P[indices[i],indices[j]]*variance*np.exp(-1j*omega*(times[i]-times[j]))
        expected=b**4*np.exp(exponent)
        error=abs(result-expected)/abs(expected)
        assert error<2e-13
        rows.append(dict(n=n,relative_four_vertex_error=error))
    return rows


def response_check(levels,coupling):
    n=3;variance=.6;omega=1.3;duration=1.1;b=1/(4*np.pi)
    ann=np.diag(np.sqrt(np.arange(1,levels)),1)
    basis=null_space(np.ones((1,n)))
    B=np.zeros((levels*levels,levels*levels),dtype=complex)
    for v in basis:
        B+=b*np.exp((1-1/n)*variance/2)*np.kron(
            expm(1j*np.sqrt(variance)*v[0]*(ann+ann.T)),
            expm(1j*np.sqrt(variance)*v[1]*(ann+ann.T)))
    X=B+B.conj().T;Y=1j*(B-B.conj().T)
    H0=omega*np.diag((np.arange(levels)[:,None]+np.arange(levels)[None,:]).ravel())
    initial=np.zeros(levels*levels,dtype=complex);initial[0]=1
    def ramp(t):
        v=t/duration;return 3*v*v-2*v*v*v
    def phase_shape(t):return .3+.7*t/duration
    phase=.0001
    observed=[]
    for sign in (-1,1):
        def rhs(t,state):
            alpha=sign*phase*phase_shape(t)
            H=H0-coupling*ramp(t)*(X*np.cos(alpha)+Y*np.sin(alpha))
            return -1j*(H@state)
        sol=solve_ivp(rhs,(0,duration),initial,method='DOP853',rtol=3e-12,atol=2e-14)
        assert sol.success
        state=sol.y[:,-1];alpha=sign*phase
        current=coupling*(Y*np.cos(alpha)-X*np.sin(alpha))
        observed.append(np.vdot(state,current@state).real)
    derivative=(observed[1]-observed[0])/(2*phase)
    def correlator(t,which):
        L=variance*np.exp(-1j*omega*t)
        if which=='X':return 4*n*b*b*(np.cosh((1-1/n)*L)+(n-1)*np.cosh(L/n)-n)
        return 4*n*b*b*(np.sinh((1-1/n)*L)-(n-1)*np.sinh(L/n))
    # H_pert=-f O implies chi=i theta<[O(t),O(0)]>=-2 Im C(t).
    local=quad(lambda s:-2*correlator(duration-s,'X').imag*coupling*ramp(s),0,duration,epsabs=1e-14)[0]
    nonlocal_y=quad(lambda s:-2*correlator(duration-s,'Y').imag*coupling*ramp(s)*phase_shape(s),0,duration,epsabs=1e-14)[0]
    first=-coupling*(2*n*b)
    complete=first-coupling*local+coupling*nonlocal_y
    incomplete=first+coupling*nonlocal_y
    error=abs(derivative-complete)
    return dict(levels=levels,coupling=coupling,exact_phase_derivative=derivative,
                first_order=first,with_both_second_order_terms=complete,
                error=error,without_local_response_error=abs(derivative-incomplete),
                local_second_order=-coupling*local,retarded_Y_second_order=coupling*nonlocal_y)


def main():
    wick=wick_cases();spectral=spectral_cases();gaussian=gaussian_matrix_checks()
    response=[response_check(10,g) for g in (.1,.05,.025)]
    refined=response_check(14,.1)
    assert abs(refined['exact_phase_derivative']-response[0]['exact_phase_derivative'])<2e-9
    for row in response:
        assert row['error']<row['without_local_response_error']/20
    assert response[1]['error']<response[0]['error']/5
    assert response[2]['error']<response[1]['error']/5
    result=dict(date='2026-10-13',scope='exact electrically neutral oscillator influence of the free LLL strip; nonlinear angular projection and full switching stress are not included',
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,
                fermion_wick_comparisons=wick,positive_spectral_cases=spectral,
                quantum_vertex_products=gaussian,quantum_response=response,
                quantum_response_refinement=refined,
                limitations=['The neutral fields retain the SU(n)_1 sector, including the not-yet-gauged angular current sector',
                             'The vertex normalization must be matched with the unitary charged factor and physical coupling',
                             'Second-order kernels are not an all-order bound at the collective physical coupling',
                             'Finite-band response checks do not define continuum stress counterterms'],passed=True)
    (ROOT/'evidence/neutral-transfer-current.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: 16 fermion determinants, 7 positive spectra, 3 quantum vertex products, neutral response and refinement')


if __name__=='__main__':
    main()
