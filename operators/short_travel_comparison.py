"""Short-travel comparison of the complete geodesic quantum form.

The physical Hessian is unchanged. Every removed long-travel contribution
belongs to the new error operator. This defines the exact scalar brackets
and certifies the full transfer cost; it does not evaluate a Hessian column.
"""
from flint import arb


class ShortTravelComparison:
    def __init__(self):
        self.radius = arb(1)/128

    def cutoff(self, distance):
        """C-infinity cutoff, one through a and zero from 2a.

        At an interval crossing an endpoint the entire [0,1] enclosure is
        returned. The exact cutoff is fixed by the formula, not by a ball's
        midpoint. This is a travel-time cutoff, not a physical boundary.
        """
        r = arb(distance)
        if not r >= 0:
            raise ValueError('Nonnegative travel distance required')
        t = r/self.radius
        if t <= 1:
            return arb(1)
        if t >= 2:
            return arb(0)
        if not (t > 1 and t < 2):
            return arb('.5', '.5')
        left = (-1/(2-t)).exp()
        right = (-1/(t-1)).exp()
        return left/(left+right)

    def brackets(self, distance, frequency):
        """Full, short and transferred spectral integrands.

        Each is multiplied by r sum_n d_n[-B_Mn(r)] and integrated over
        all r>0 in the exact definition. The original fourth finite part
        is additional and unchanged. No angular or physical thermal mode
        has been selected by this function.
        """
        r, w = arb(distance), arb(frequency)
        zeta = self.cutoff(r)
        x = r*w
        difference = 1-x.cos()
        local = -x*x/2+x**4/24
        return dict(full=difference+local,
                    short=zeta*difference+local,
                    transferred=(1-zeta)*difference)

    def full_transfer_L2_upper(self, C):
        """Complete physical-half base-L2 coefficient, including 2pi angle.

        Uses the retained entire-tower bound sum d_n[-B_Mn(r)]<=72C/r^6,
        spectral 0<=1-cos<=2 and all travel times r>=a.
        """
        C = arb(C)
        if not C > 0:
            raise ValueError('Positive sphere/orbit ratio required')
        return 72*arb.pi()*C/self.radius**4

    def physical_form_transfer_upper(self, C, B2, lam, F_lambda):
        """Current physical V->V* bound; original action units retained.

        This imports C J0/E0 <=1000 C^2/(lambda B^4) and
        l>=.99 F_lambda E0. It is NOT a strong-map or admission bound.
        """
        C, B2, lam, F_lambda = map(arb, (C, B2, lam, F_lambda))
        if not (C > 0 and B2 > 0 and lam > 0 and F_lambda > 0):
            raise ValueError('Positive physical parameters required')
        j0 = 1000*C*C/(lam*B2*B2)
        return 72*arb.pi()*j0/(self.radius**4*arb('.99')*F_lambda)
