import numpy as np
import matplotlib.pyplot as plt
from .Circles_Lines import Circle

class Mobius:

def __init__(self, a=1.0, b=0.0, c=1.0, d=1.0, name=''):
    """Defines a Mobius transformation given four coefficients"""

    # Cast to clongdouble to ensure precision
    self.a = np.clongdouble(a)
    self.b = np.clongdouble(b)
    self.c = np.clongdouble(c)
    self.d = np.clongdouble(d)

    self.name = name


    def det(self):
        """Returns determinant of transformation"""
        return (self.a * self.d) - (self.b * self.c)

    def trace(self):
        """Returns trace of transformation"""
        return self.a + self.d

    def invert(self):
        """Returns inverse of transformation"""
        return Mobius(self.d, -self.b, -self.c, self.a, name='('+self.name+')^{-1}')

    def mat(self):
        """Returns matrix of transformation"""
        return np.array([[self.a, self.b], [self.c, self.d]], dtype=np.clongdouble)


    def get_isometric_circle(self):
        """Returns isometric circle of transformation.

        Note: If T(z) = (az + b) / (cz + d) is transformation, then complex equation of isometric circle is |cz + d| = 1, i.e. |z + d/c| = 1/c.
        """

        center = np.clongdouble(-self.d/self.c)
        radius = np.absolute(np.clongdouble(1.) / self.c, dtype = np.longdouble)

        circle = Circle(center, radius)

        return circle

    def plot_isometric_circle(self):
        """Plots isometric circle is matplotlib plot."""

        circle = self.get_isometric_circle()
        circle.plot()


    def scale(self, scalar):
        """Scales transformation as a function. If T(z) = (az +b)/(cz + d), then kT(z) = (kaz + kb)/(cz + d) where k is the scalar."""

        scalar = np.clongdouble(scalar)
        a = scalar * self.a
        b = scalar * self.b

        return Mobius(a, b, self.c, self.d)

    def normalize(self):
        """Normalizes transformation to determinant 1."""
        scalar = np.clongdouble(1. / np.sqrt(self.det(), dtype=np.clongdouble))
        a = self.a * scalar
        b = self.b * scalar
        c = self.c * scalar
        d = self.d * scalar
        name = self.name + '_norm'

        return Mobius(a, b, c, d)

    def fixed_points(self):
        """Returns fixed points of transformation."""
        fixed1 =  ((self.a-self.d) + np.sqrt(self.trace()**2 - 4*self.det())) / (2*self.c)
        fixed2 = ((self.a-self.d) - np.sqrt(self.trace()**2 - 4*self.det())) / (2*self.c)

        return fixed1, fixed2

    def set_name(self, name):
        """Sets name of transformation"""
        self.name = name


    def __call__(self, x, evalf=False):

        x = np.clongdouble(x)
        return ((self.a * x) + self.b) / ((self.c*x) + self.d)


    def __eq__(self, M):
        norm1 = self.normalize()
        norm2 = M.normalize()

        coeff_pairs = [(norm1.a, norm2.a), (norm1.b, norm2.b), (norm1.c, norm2.c), (norm1.d, norm2.d)]
        coeff_diffs = []
        coeff_sums = []

        for p in coeff_pairs:
            coeff_diffs.append(p[0] - p[1])
            coeff_sums.append(p[0] + p[1])

        coeff_diffs = np.array(coeff_diffs, dtype=np.clongdouble)
        coeff_sums = np.array(coeff_sums, dtype=np.clongdouble)

        if np.all(isclose(coeff_diffs, 0)) or np.all(isclose(coeff_sums, 0)):
            return True
        else:
            return False

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @staticmethod
    def compose(M1, M2):
        """Returns composition of two transformation via matrix multiplication formula."""
        a_ = M1.a * M2.a + M1.b * M2.c
        b_ = M1.a * M2.b + M1.b * M2.d
        c_ = M1.c * M2.a + M1.d * M2.c
        d_ = M1.c * M2.b + M1.d * M2.d

        name = M1.name + '-' + M2.name

        return Mobius(a_, b_, c_, d_, name)

    def compose_list(elements):
        """Returns composition of list of transformations in list order, i.e. if elements = [T1, T2, T3], then return T1-T2-T3, where - represents composition."""

        comp = Mobius(1,0,0,1)

        for i in range(len(elements)-1, -1, -1):
            comp = Mobius.compose(elements[i], comp)

        comp.name = comp.name[:-1] # Gets rid of dash at the end due to composing with identity

        return comp

    @staticmethod
    def geodesic_to_imag(center, radius, max_image=1j):
        """ Given geodesic in UHP (Euclidean circle with real center), return Mobius transformation sending geodesic to imaginary axis
        # max_image is where the top of Euclidean circle gets sent to on imaginary axis."""

        scalar = np.clongdouble(max_image / 1j)

        T = Mobius(-1, center - radius, 1, -(center + radius)).scale(scalar)
        T = T.normalize()

        return T

    @staticmethod
    def compare_list(G1, G2):
        G2 = list(G2)   # make a mutable copy
        try:
            for T in G1:
                G2.remove(T)
        except ValueError:
            return False

        return not G2
