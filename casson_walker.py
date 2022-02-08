import sympy as sym
import hf_nemethi as hf

##############################################################################################################################################################################################

def casson_walker(C):
    '''Computes the Casson-Walker invariant of a Seifert fibered space specified by a definite plumbing.
     
          Args:
               C(list of int): Coefficients specifying a definite plumbing.

          Returns:
               casson_walker(sym.Rational): The Casson-Walker invariant of the Seifert fibered space.
    '''
    for q in [sym.gcd(C[i],C[i+1]) for i in range(1,len(C),2)]:
        if q != 1:
            raise Exception('The framings on the exceptional fibers should be reduced fractions!')
    S = hf.linking_matrix(C)
    if not max(S.is_negative_definite, S.is_positive_definite):
        raise Exception('The corresponding plumbing is not definite!')
    epsilon = sym.sign(C[0])
    if epsilon == 1:
        C = [-abs(i) for i in C]
        for i in range(2,len(C),2):
            C[i] = abs(C[i])
    S = hf.linking_matrix(C)
    if sym.det(S) == 0:
        raise Exception('The corresponding plumbing should specify a rational homology sphere!')
    T, n = S.inv(), abs(sym.det(S))
    return epsilon*sym.Rational(1,24)*(S.trace() + 3*sym.shape(S)[0] + sum((2 - (sum(S.row(i)) - S[i,i]))*T[i,i] for i in range(sym.shape(S)[0])))

def dedekind_sum(m,n):
    '''Computes the Dedekind sum s(m,n).
     
          Args:
               m(int): The first parameter of s(m,n), a non-zero integer.
               n(int): The second parameter of s(m,n), a non-zero integer.

          Returns:
               dedekind_sum(sym.Rational): The Dedekind sum s(m,n).
    '''
    frac = lambda x : 0 if x%1 == 0 else x-sym.floor(x)-sym.Rational(1,2)
    return sum(frac(sym.Rational(i,n))*frac(sym.Rational(i*m,n)) for i in range(1, n))

def casson_walker_lens(p,q):
    '''Computes the Casson-Walker invariant of the lens space L(p,q), for p and q coprime.
     
          Args:
               p(int): The first parameter of L(p,q), a non-zero integer.
               q(int): The second parameter of L(p,q), a non-zero integer.

          Returns:
               casson_walker_lens(sym.Rational): The Casson-Walker invariant of L(p,q).
    '''
    if sym.gcd(p,q) != 1:
        raise Exception('The parameters of L(p,q) should be coprime!')
    if p*q == 0:
        raise Exception('The parameters should be non-zero integers!')
    epsilon = sym.sign(p)
    p, q = abs(p), q%abs(p)
    if epsilon == -1:
        q = p-q
    return sym.Rational(dedekind_sum(q,p), 2)

def casson_walker_prism(p,q):
    '''Computes the Casson-Walker invariant of the prism manifold P(p,q), for p and q coprime.
     
          Args:
               p(int): The first parameter of P(p,q), a non-zero integer.
               q(int): The second parameter of P(p,q), a non-zero integer.

          Returns:
               casson_walker_prism(sym.Rational): The Casson-Walker invariant of P(p,q).
    '''
    if sym.gcd(p,q) != 1:
        raise Exception('The parameters of P(p,q) should be coprime!')
    if p*q == 0:
        raise Exception('The parameters should be non-zero integers!')
    if q < 0:
        p, q = -p, -q
    return sym.Rational(1,2)*(sym.Rational(p,8*q)-dedekind_sum(p,q))

##############################################################################################################################################################################################

# Implementation of Rustamov's formula for the Casson-Walker invariant (where chi is an auxiliary function to compute the Euler characteristic of HF^+):

# import continued_fractions as cf

# def chi(dict):
#     even, odd = 0, 0
#     for key in dict.keys():
#         if key == 0:
#                 continue
#         even += key*len([i for i in dict[key] if i%2 == 0])
#         odd += key*len([i for i in dict[key] if i%2 == 1])
#     return even - odd

# def casson_walker(C):
#     '''Computes the Casson-Walker invariant of a Seifert fibered space specified by a definite plumbing as in hf.
     
#           Args:
#                C(list): Coefficients specifying a definite plumbing.

#           Returns:
#                casson_walker(sym.Rational): The Casson-Walker invariant of the Seifert fibered space.
#     '''
#     res = hf.spinc_to_HF(C)
#     n = len(res)
#     chi_hat = lambda s : chi(res[s]) - sym.Rational(res[s][0][0], 2)
#     return sym.Rational(1, n)*sum(chi_hat(s) for s in res.keys())

# Alternative formula for computing the Casson-Walker invariant of a lens space.

# def casson_walker_lens(p,q):
#     '''Computes the Casson-Walker invariant of the lens space L(p,q), for p and q coprime.
     
#           Args:
#                p(int): The first parameter of L(p,q), a non-zero integer.
#                q(int): The second parameter of L(p,q), a non-zero integer.

#           Returns:
#                casson_walker_lens(sym.Rational): The Casson-Walker invariant of L(p,q).
#     '''
#     if sym.gcd(p,q) != 1:
#         raise Exception('The parameters of L(p,q) should be coprime!')
#     if p*q == 0:
#         raise Exception('The parameters should be non-zero integers!')
#     if abs(p) == 1:
#         return 0
#     epsilon = sym.sign(p)
#     p, q = abs(p), q%abs(p)
#     if epsilon == -1:
#         q = p-q
#     L = cf.neg_continued_fraction_coefficients(p,q)
#     return sym.Rational(1,24)*(sym.Rational(q,p)+sym.Rational(sym.mod_inverse(q,p),p)+sum(a-3 for a in L))

##############################################################################################################################################################################################