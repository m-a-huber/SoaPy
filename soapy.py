import sympy as sym
import hf_nemethi as hf
import casson_walker as cw
import continued_fractions as cf
import ribbon_obstruction as rib

##############################################################################################################################################################################################

class SFS:
    '''This is a class for working with orientable Seifert fibered spaces (SFS) whose base orbifold is the 2-sphere. 
     
        Attributes:
            params (list of int): List of integer coefficients representing the SFS specified. These do not necessarily coincide with the input parameters, but rather are normalized in such a way that the corresponding integer surgery diagram is definite.
            central_weight (int): The weight of the central vertex of the normalized surgery description of the SFS specified.
            branch_weights (tuple of sym.Rational): Tuple containing the rational surgery coefficients of the exceptional fibers of the SFS specified.
            fractional_branch_weights (tuple of sym.Rational): Tuple containing the fractional parts of the branch weights.
            euler_number (sym.Rational): The orbifold Euler number of the normalized surgery description of the SFS specified.
            exceptional_fibers (int): The number of exceptional fibers of the normalized surgery description of the SFS specified.
    '''
    
    def __init__(self, *params):
        '''The constructor for the SFS class.
     
            Parameters:
                params (list of int): A variable number of integers specifying a SFS. The list must be of the format (e, a_1, b_1, ..., a_n, b_n), where e is the integer framing of the central unknot, and the a_i/b_i are the rational framings of the exceptional fibers (expressed as a reduced fraction). If the input list specifies a SFS with non-zero first Betti number, an error is raised.
        '''
        cond_1 = all(isinstance(el, int) for el in params)
        cond_2 = len(params)%2 == 1
        cond_3 = all([sym.gcd(params[i], params[i+1]) == 1 for i in range(1,len(params),2)])
        cond_4 = all(el!=0 for el in params[1:])
        cond_5 = hf.normalize(params)[0] != 0
        if not min(cond_1, cond_2, cond_3, cond_4):
            raise Exception('A SFS should be specified by an odd number of integers, with each exceptional fiber specified by a pair of non-zero coprime integers!')
        if not cond_5:
            raise Exception('The SFS specified is not a rational homology sphere!')
        self.params = tuple(hf.normalize(params)[1])
        self.central_weight = self.params[0]
        self.branch_weights = tuple(sym.Rational(self.params[i], self.params[i+1]) for i in range(1,len(self.params),2))
        self.fractional_branch_weights = tuple((1/w) - sym.floor(1/w) for w in self.branch_weights)
        self.euler_number = self.central_weight - sum([1/w for w in self.branch_weights])
        self.exceptional_fibers = len(self.branch_weights)

    def __repr__(self):
        branch_weights_string = ', '.join([str(q) for q in self.branch_weights])
        return 'Y({}; {})'.format(self.central_weight, branch_weights_string)
    
    def __neg__(self):
        params = tuple(-((-1)**i)*self.params[i] for i in range(len(self.params)))
        return SFS(*params)
    
    def __eq__(self, other):
        if isinstance(other, SFS):
            if not max(self.is_lens_space(), other.is_lens_space()):
                seifert_self = (self.euler_number, sorted(self.fractional_branch_weights))
                seifert_other = (other.euler_number, sorted(other.fractional_branch_weights))
                return seifert_self == seifert_other
            if self.is_lens_space() and other.is_lens_space():
                self, other = self.to_lens_space(), other.to_lens_space()
                return all([self.p == other.p, max(self.q == other.q, (self.q*other.q)%self.p == 1)])
        return False
    
    def __le__(self, other):
        if isinstance(other, SFS):
            if self == other:
                return True
            return rib.ribbon_obstruction_sfs(self.params, other.params)
        return False
    
    def __lt__(self, other):
        if isinstance(other, SFS):
            if self == other:
                return False
            return self <= other
        return False
    
    def seifert_invariants(self):
        '''Returns the Seifert invariants of the SFS specified.
     
            Parameters:
                None
            
            Returns:
                seifert_invariants (tuple): A tuple of the format (Euler number, (tuple of fractional branch weights)).
        '''
        return (self.euler_number, self.fractional_branch_weights)
    
    def linking_matrix(self):
        '''Returns the linking matrix of the SFS specified.
     
            Parameters:
                None
            
            Returns:
                linking_matrix (sym.Matrix): A SymPy-matrix with SymPy-integers as entries, representing the linking matrix of the integer plumbing corresponding to the SFS specified.
        '''
        return hf.linking_matrix(self.params)
    
    def first_homology(self):
        '''Returns the order of the first homology of the SFS specified.
     
            Parameters:
                None
            
            Returns:
                first_homology (int): The order of the first homology of the SFS specified.
        '''
        return int(abs(sym.det(self.linking_matrix())))
    
    def spinc_to_HF(self):
        '''Computes HF^+ in each spin^c-structure of the SFS specified. The Z[U]-module-structure of HF^+ is encoded as a dictionary of the format {'order of Z[U]-module-summand' : 'list of bottommost gradings of all Z[U]-module-summands of that order'}.
     
            Parameters:
                None

            Returns:
                spinc_to_HF(dict): A dictionary of the format {'spin^c-structure' : 'Z[U]-module-structure of HF^+'}.
        '''
        return hf.spinc_to_HF(self.params)
    
    def print_HF(self):
        '''Prints HF^+ of the SFS specified by a definite plumbing in a more legible manner.
     
            Parameters:
                None

            Returns:
                None: Prints HF^+ of the Seifert fibered space specified.
        '''
        print('HF^+({}):'.format(self))
        hf.print_HF(self.params)
        return
    
    def correction_terms(self):
        '''Returns a list of the corrections terms of the SFS specified.
     
            Parameters:
                None

            Returns:
                correction_terms (list of sym.Rational): List of all correction terms of the SFS specified.
        '''
        return tuple(sorted(hf.correction_terms(self.params)))
    
    def is_lspace(self):
        '''Checks whether or not the SFS specified is a Heegaard Floer L-space.
     
            Parameters:
                None

            Returns:
                bool: Whether or not the the SFS specified is an L-space.
        '''
        return hf.is_lspace(self.params)
    
    def casson_walker(self):
        '''Computes the Casson-Walker invariant of the SFS specified.
     
            Parameters:
                None

            Returns:
                casson_walker (sym.Rational): The Casson-Walker invariant of the SFS specified.
        '''
        if self.is_lens_space():
            self = self.to_lens_space()
            return cw.casson_walker_lens(self.p, self.q)
        if self.is_prism_mfld():
            self = self.to_prism_mfld()
            return cw.casson_walker_prism(self.p, self.q)
        return cw.casson_walker(self.params)

    def is_lens_space(self):
        '''Checks whether or not the SFS specified is homeomorphic to a lens space.
     
            Parameters:
                None

            Returns:
                bool: Whether or not the the SFS specified is homeomorphic to a lens space.
        '''
        return self.exceptional_fibers <= 2
    
    def to_lens_space(self):
        '''Transforms the SFS specified into the corresponding lens space, provided it is homeomorphic to one. If the SFS specified is not homeomorphic to any lens space, an error is raised.
     
            Parameters:
                None

            Returns:
                to_lens_space (soapy.Lens): Lens space homeomorphic to the SFS specified.
        '''
        if not self.is_lens_space():
            raise Exception('The SFS specified is not homeomorphic to a lens space!')
        coeffs = [sym.Rational(self.central_weight)]
        if self.exceptional_fibers >= 1:
            coeffs += cf.number_to_neg_cont_frac(self.branch_weights[0].p, self.branch_weights[0].q)
        if self.exceptional_fibers == 2:
            coeffs.reverse()
            coeffs += cf.number_to_neg_cont_frac(self.branch_weights[1].p, self.branch_weights[1].q)
        frac = -cf.neg_cont_frac_to_number(*coeffs)
        return Lens(frac.p, frac.q)
    
    def is_prism_mfld(self):
        '''Checks whether or not the SFS specified is homeomorphic to a prism manifold.
     
            Parameters:
                None

            Returns:
                bool: Whether or not the the SFS specified is homeomorphic to a prism manifold.
        '''
        return self.fractional_branch_weights.count(sym.Rational(1,2)) >= 2 and self.exceptional_fibers == 3
    
    def to_prism_mfld(self):
        '''Transforms the SFS specified into the corresponding prism manifold, provided it is homeomorphic to one. If the SFS specified is not homeomorphic to any prism manifold, an error is raised.
     
            Parameters:
                None

            Returns:
                to_prism_mfld (soapy.Prism): Prism manifold homeomorphic to the SFS specified.
        '''
        if not self.is_prism_mfld():
            raise Exception('The SFS specified is not homeomorphic to a prism manifold!')
        return Prism((self.euler_number).q, (self.euler_number).p)

three_sphere = SFS(1)
poincare_sphere = SFS(-2,-5,4,-3,2,-2,1)

##############################################################################################################################################################################################

class Lens(SFS):
    def __init__(self, p, q):
        cond_1 = all(isinstance(el, int) for el in (p,q))
        cond_2 = sym.gcd(p,q) == 1
        cond_3 = p*q != 0
        if not min(cond_1, cond_2, cond_3):
            raise Exception('A lens space should be specified by a pair of non-zero coprime integers!')
        super().__init__(*hf.lens(p,q))
        self.p = abs(p)
        self.q = (lambda x : 1 if x==0 else x)((lambda q : q%p if p>0 else (self.p-q)%self.p)(q))
    
    def __repr__(self):
        return 'L({}, {})'.format(self.p,self.q)
    
    def __neg__(self):
        return Lens(self.p, -self.q)

    def to_SFS(self):
        return SFS(*hf.lens(self.p,self.q))

    def linear_lattice(self, epsilon=-1):
        if epsilon not in  {-1, 1}:
            raise Exception('The second (optional) argument should be plus/minus 1!')
        if epsilon == -1:
            return cf.number_to_neg_cont_frac(-self.p, self.q)
        return cf.number_to_neg_cont_frac(self.p, self.p - self.q)

##############################################################################################################################################################################################

class Prism(SFS):
    def __init__(self, p, q):
        cond_1 = all(isinstance(el, int) for el in (p,q))
        cond_2 = abs(p) > 1
        cond_3 = sym.gcd(p,q) == 1
        if not min(cond_1, cond_2, cond_3):
            raise Exception('A prism manifold should be specified by a pair of coprime integers, with the first parameter greater than 1 in absolute value!')
        super().__init__(*hf.prism(p,q))
        self.p = abs(p)
        self.q = int(sym.sign(p)*q)
    
    def __repr__(self):
        return 'P({}, {})'.format(self.p,self.q)
    
    def __neg__(self):
        return Prism(self.p, -self.q)
    
    def to_SFS(self):
        return SFS(*hf.prism(self.p,self.q))

##############################################################################################################################################################################################

class Brieskorn(SFS):
    def __init__(self, *params):
        cond_1 = all(isinstance(el, int) for el in params)
        cond_2 = min(map(abs, params)) > 1
        cond_3 = len(set(map(sym.sign, params))) == 1
        cond_4 = sym.lcm_list(params) == sym.prod(params)
        if not min(cond_1, cond_2, cond_3, cond_4):
            raise Exception('A Brieskorn sphere should be specified by a list of pairwise coprime integers, each greater than 1 in absolute value, and all of the same sign!')
        super().__init__(*hf.brieskorn(*params))
        self.coeffs = params

    def __repr__(self):
        epsilon = sym.sign(self.coeffs[0])
        sign = (lambda x : '' if x > 0 else '-')(epsilon)
        coeffs_string = ','.join([str(abs(a)) for a in self.coeffs])
        return '{}Sigma({})'.format(sign, coeffs_string)

    def __neg__(self):
        return Brieskorn(*map(lambda x : -x, self.coeffs))
    
    def to_SFS(self):
        return SFS(*hf.brieskorn(*self.coeffs))

##############################################################################################################################################################################################