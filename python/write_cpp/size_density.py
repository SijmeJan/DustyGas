import numpy as np
from scipy.integrate import quad

class SizeDensity():
    def __init__(self, p, Stokes_range):
        """Power-law size density"""
        self.p = p
        self.Stokes_range = Stokes_range
        self.s_range = np.log(np.asarray(Stokes_range))

    def sigma(self, Stokes):
        # Make sure we can handle both vector and scalar Stokes
        scalar_input = False
        if Stokes.ndim == 0:
            Stokes = Stokes[None]  # Makes 1D
            scalar_input = True
        else:
            original_shape = np.shape(Stokes)
            Stokes = np.ravel(Stokes)

        ret = self.p*np.power(Stokes, -self.p)/(np.power(self.Stokes_range[0], -self.p) - np.power(self.Stokes_range[1], -self.p))

        # Return value of original shape
        if scalar_input:
            return np.squeeze(ret)
        return np.reshape(ret, original_shape)

    def Jint(self, alpha, mu):
        f = lambda x: self.sigma(np.exp(x))*np.exp(alpha*x)/(1 + np.exp(2*x))

        ret = quad(f, self.s_range[0], self.s_range[1])[0]

        return ret*mu
