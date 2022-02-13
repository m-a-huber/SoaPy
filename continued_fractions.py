import sympy as sym

def number_to_neg_cont_frac(p,q):
     '''This computes the coefficients of the (-)-continued fraction expansion of a rational number p/q, for p, q coprime.
     
          Args:
               p(int): The numerator of the rational number p/q.
               q(int): The denominator of the rational number p/q.

          Returns:
               L(list of int): The coefficients of the (-)-continued fraction expansion of a rational number p/q.
     '''
     L = []
     epsilon = int(sym.sign(p*q))
     if epsilon < 0:
          p, q = abs(p), abs(q)
     while q != 0:
          L.append(sym.ceiling(p/q))
          p, q = q, sym.ceiling(p/q)*q - p
     return list(map(int, [epsilon*l for l in L]))

def number_to_pos_cont_frac(p,q):
     '''This computes the coefficients of the (+)-continued fraction expansion of a rational number p/q, for p, q coprime.
     
          Args:
               p(int): The numerator of the rational number p/q.
               q(int): The denominator of the rational number p/q.

          Returns:
               L(list of int): The coefficients of the (+)-continued fraction expansion of a rational number p/q.
     '''
     L = []
     epsilon = sym.sign(p*q)
     if epsilon < 0:
          p, q = abs(p), abs(q)
     while q != 0:
          L.append(sym.floor(p/q))
          p, q = q, p - sym.floor(p/q)*q
     return list(map(int, [epsilon*l for l in L]))

def number_from_neg_cont_frac(*args):
     '''This computes (-)-continued fraction of an arbitrary list of coefficients.

          Args:
               *args(int): A variable number of integers which are the coefficients of a (-)-continued fraction decomposition.
          
          Returns:
               x(sym.Rational): The rational number p/q with corresponding coefficients in its (-)-continued fraction decomposition.
     '''
     x = args[-1]
     for i in args[-2::-1]:
          x = i - sym.Rational(1, x)
     return sym.Rational(x)

def number_from_pos_cont_frac(*args):
     '''This computes (+)-continued fraction of an arbitrary list of coefficients.

          Args:
               *args(int): A variable number of integers which are the coefficients of a (+)-continued fraction decomposition.
          
          Returns:
               x(sym.Rational): The rational number p/q with corresponding coefficients in its (+)-continued fraction decomposition.
     '''
     x = args[-1]
     for i in args[-2::-1]:
          x = i + sym.Rational(1, x)
     return sym.Rational(x)