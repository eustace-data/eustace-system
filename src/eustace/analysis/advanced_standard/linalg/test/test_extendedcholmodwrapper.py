"""Quick tests for extended cholmod wrapper features."""

from eustace.analysis.advanced_standard.linalg.extendedcholmodwrapper import ExtendedCholmodWrapper

import unittest
import numpy
import sksparse.cholmod
from scipy.sparse import csr_matrix

class TestExtendedCholmodWrapper(unittest.TestCase):

    def setUp(self):

        # -- Simulate problem

        # lower triangular 4x4 matrix
        L = csr_matrix(
            [ [ 2.5, 0.0, 0.0, 0.0 ],
              [ 0.3, 1.4, 0.0, 0.0 ],
              [ 0.0, 0.0, 3.8, 0.0 ],
              [ 0.2, 0.0, 0.0, 2.2 ] ] )

        # corresponding A = L . Lt
        self.A = L * L.T

        # simulate solution that had two right-hand-sides
        self.x = numpy.array([ [ 1.0,  900.0 ],
                               [ 2.0, -230.0 ],
                               [ 3.0,   22.2 ],
                               [ 4.0, -250.0 ] ])

        # simulate the right-hand sides
        self.b = self.A * self.x             # A.x  = b
        self.c = L * self.x                  # L.x  = c
        self.d = L.T * self.x                # Lt.x = d

        # -- Base class from cholmod library
        self.factor = sksparse.cholmod.cholesky(self.A)

        # Operations of extended class act on the fill-reducing permutation 
        # that was used during the solve
        self.permutation = self.factor.P()

    def test_base(self):

        # This should solve ok using base class
        # If this fails maybe the environment isn't configured correctly?
        numpy.testing.assert_almost_equal(self.x, self.factor.solve_A(self.b))

    def test_solve_A(self):

        wrapper = ExtendedCholmodWrapper.cholesky(self.A)
        numpy.testing.assert_almost_equal(self.x, wrapper.solve_A(self.b))

    def test_P(self):

        numpy.testing.assert_equal(ExtendedCholmodWrapper.cholesky(self.A).P(), self.factor.P())

    def test_apply_P(self):
     
        b = numpy.array([ [ 1.0 , 3.0 ,  2.0 , 7.0 ] ]).T
        numpy.testing.assert_almost_equal(
            ExtendedCholmodWrapper.cholesky(self.A).apply_P(b),
            self.factor.apply_P(b))

    def test_apply_Pt(self):
     
        b = numpy.array([ [ 1.0 , 3.0 , 2.0 , 7.0 ] ]).T
        numpy.testing.assert_almost_equal(
            ExtendedCholmodWrapper.cholesky(self.A).apply_Pt(b),
            self.factor.apply_Pt(b))

    def test_logdet(self):

        numpy.testing.assert_almost_equal(ExtendedCholmodWrapper.cholesky(self.A).logdet(), self.factor.logdet())

    def test_cholesky_inplace(self):

        # Initialise with identity instead of the intended A matrix
        wrapper = ExtendedCholmodWrapper.cholesky(csr_matrix(numpy.eye(4)))

        # Now set the matric using the inplace method
        wrapper.cholesky_inplace(self.A)

        # And should solve as before
        numpy.testing.assert_almost_equal(self.x, wrapper.solve_A(self.b))

    def test_solve_L_vector(self):

        wrapper = ExtendedCholmodWrapper.cholesky(self.A)

        # using individual vector
        v = self.c[self.permutation,0].ravel()
        numpy.testing.assert_almost_equal(self.x[self.permutation,0], wrapper.LLt_solve_L(v))

    def test_solve_Lt_vector(self):

        wrapper = ExtendedCholmodWrapper.cholesky(self.A)

        # using individual vector
        w = self.d[self.permutation,0].ravel()
        numpy.testing.assert_almost_equal(self.x[self.permutation,0], wrapper.LLt_solve_Lt(w))

    def test_solve_L_matrix(self):

        wrapper = ExtendedCholmodWrapper.cholesky(self.A)

        # using multiple right-hand sides
        numpy.testing.assert_almost_equal(self.x[self.permutation,:], wrapper.LLt_solve_L(self.c[self.permutation,:]))

        # we should see that this also solves just the same as using the L from original factor
        # (which already has permutation applied)
        x = numpy.matrix( [ 9.0, 8.0, 7.0, 5.0 ] ).T
        b =  self.factor.L() * x
        numpy.testing.assert_almost_equal(x,  wrapper.LLt_solve_L(b))

    def test_solve_Lt_matrix(self):

        wrapper = ExtendedCholmodWrapper.cholesky(self.A)

        # using multiple right-hand sides
        numpy.testing.assert_almost_equal(self.x[self.permutation,:], wrapper.LLt_solve_Lt(self.d[self.permutation,:]))
        
        # we should see that this also solves just the same as using the Lt from original factor
        # (which already has permutation applied)
        x = numpy.matrix( [ 9.0, 8.0, 7.0, 5.0 ] ).T
        b =  self.factor.L().T * x
        numpy.testing.assert_almost_equal(x,  wrapper.LLt_solve_Lt(b))
