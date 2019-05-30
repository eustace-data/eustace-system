"""Base class to wrap different linear algebra libraries for solving sparse positive-definite system."""

import numpy

class SolverLibraryWrapper(object):

    @staticmethod
    def as_matrix(b):
        """One-dimensional vectors require re-arranging as matrix objects for cholmod to work."""
      
        if isinstance(b, numpy.ndarray):

            if (b.ndim == 1) or ((b.ndim == 2) and (b.shape[1] == 1)):

                b = numpy.matrix(b.ravel()).T
                
        return b

    @staticmethod
    def return_like(exemplar, result):
        """Where appropriate this does the reverse of as_matrix and converts back to 1D."""

        if isinstance(exemplar, numpy.ndarray):

            if exemplar.ndim == 1:

                result = result.ravel()

        return result

    def solve_A(self, b):
        """Call underlying factor.solve_A method to solve A.x = b."""

        raise NotImplementedError

    def P(self):
        """Underlying permutation used."""

        raise NotImplementedError

    def solve_backward_substitution(self, b):
        """Returns :math:`x`, where :math:`L'x = b` for an LL' factorisation."""

        raise NotImplementedError
