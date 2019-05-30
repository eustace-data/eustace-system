"""Extend cholmod wrapper functionality."""

import numpy
import sksparse.cholmod
from wrapper import SolverLibraryWrapper

class ExtendedCholmodWrapper(SolverLibraryWrapper):
    """Extend behaviour of cholesky factor."""

    CHOLMOD_L = 4
    """Copied from definition in underlying C library."""

    CHOLMOD_Lt = 5
    """Copied from definition in underlying C library."""

    @staticmethod
    def cholesky(A, printstats=False):
        """Construct wrapper by making cholmod factor for matrix A."""

        factor = sksparse.cholmod.cholesky(A)
        return ExtendedCholmodWrapper(factor)
            
    def __init__(self, factor):
        """Construct from cholmod Factor object."""

        self.factor = factor

    def cholesky_inplace(self, A):
        """Call underlying factor.cholesky_inplace method."""

        return self.factor.cholesky_inplace(A)

    def solve_A(self, b):
        """Call underlying factor.solve_A method to solve A.x = b."""

        result = self.factor.solve_A(ExtendedCholmodWrapper.as_matrix(b))
        return ExtendedCholmodWrapper.return_like(b, result)

    def solve_backward_substitution(self, b):
        """Call underlying factor.LLt_solve_Lt to solve LT.x = b where A = L.LT"""

        b_permuted = self.apply_P(b)
        result_permuted = self.LLt_solve_Lt(b_permuted)
        return self.apply_Pt(result_permuted)

    def P(self):
        """Call underlying factor.P() to get permutation."""

        return self.factor.P()

    def logdet(self):
        """Call underlying factor.logdet."""

        return self.factor.logdet()

    def LLt_vector_instruction(self, b, instruction):
        """Issue solver instruction to operate on vector."""

        self.factor._ensure_L_or_LD_inplace(True)
        result = self.factor._solve(ExtendedCholmodWrapper.as_matrix(b), instruction)
        return ExtendedCholmodWrapper.return_like(b, result)

    def LLt_solve_L(self, b):
        """Returns :math:`x`, where :math:`Lx = b` for an LL' factorisation."""
        
        return self.LLt_vector_instruction(b, ExtendedCholmodWrapper.CHOLMOD_L)

    def LLt_solve_Lt(self, b):
        """Returns :math:`x`, where :math:`L'x = b` for an LL' factorisation."""

        return self.LLt_vector_instruction(b, ExtendedCholmodWrapper.CHOLMOD_Lt)

    def apply_P(self, b):
        """Call underyling factor.apply_P"""

        return ExtendedCholmodWrapper.return_like(b, self.factor.apply_P(ExtendedCholmodWrapper.as_matrix(b)))

    def apply_Pt(self, b):
        """Call underyling factor.apply_Pt"""

        return ExtendedCholmodWrapper.return_like(b, self.factor.apply_Pt(ExtendedCholmodWrapper.as_matrix(b)))
