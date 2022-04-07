import sympy as sym
from . import hf_nemethi as hf

def ribbon_obstruction_corr(L,M):
     """This checks whether a pair of 3-manifolds with given lists of correction terms passes the d-invariant obstruction to the existence of a ribbon cobordism between them.

     Args:
         L (list[int]): List of the correction terms of the first 3-manifold.
         M (list[int]): List of the correction terms of the second 3-manifold.

     Returns:
         bool: Whether or not the pair passes the obstruction.
     """     
     u = sym.Rational(len(M),len(L))
     if u.is_integer and sym.ntheory.primetest.is_square(u):
          for d in L:
               if M.count(d) != sym.sqrt(u)*L.count(d):
                    return False
          return True
     else:
          return False
    
def ribbon_obstruction_sfs(C_1,C_2):
     """This checks whether a pair of Seifert fibered spaces passes the d-invariant obstruction to the existence of a ribbon cobordism between them.

     Args:
         C_1 (list[int]): List of coefficients specifying the first Seifert fibered space.
         C_2 (list[int]): List of coefficients specifying the second Seifert fibered space.

     Returns:
         bool: Whether or not the pair passes the obstruction.
     """     
     L, M = hf.correction_terms(C_1), hf.correction_terms(C_2)
     return ribbon_obstruction_corr(L,M)

def ribbon_obstruction_lens(p,q,r,s):
     """This checks whether a pair of lens spaces passes the d-invariant obstruction to the existence of a ribbon cobordism between them.

     Args:
         p (int): First parameter of first lens space.
         q (int): Second parameter of first lens space.
         r (int): First parameter of second lens space.
         s (int): Second parameter of second lens space.

     Returns:
         bool: Whether or not the pair passes the obstruction.
     """     
     L, M = hf.correction_terms(hf.lens(p,q)), hf.correction_terms(hf.lens(r,s))
     return ribbon_obstruction_corr(L,M)

def ribbon_obstruction_prism(p,q,r,s):
     """This checks whether a pair of prism manifolds passes the d-invariant obstruction to the existence of a ribbon cobordism between them.

     Args:
         p (int): First parameter of first prism manifold.
         q (int): Second parameter of first prism manifold.
         r (int): First parameter of second prism manifold.
         s (int): Second parameter of second prism manifold.

     Returns:
         bool: Whether or not the pair passes the obstruction.
     """     
     L, M = hf.correction_terms(hf.prism(p,q)), hf.correction_terms(hf.prism(r,s))
     return ribbon_obstruction_corr(L,M)