import sympy as sym
import hf_nemethi as hf
import casson_walker as cw
import continued_fractions as cf
import ribbon_obstruction as rib

##############################################################################################################################################################################################

class SFS:
    
    def __init__(self, *params):
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
        return (self.euler_number, self.fractional_branch_weights)
    
    def linking_matrix(self):
        return hf.linking_matrix(self.params)
    
    def first_homology(self):
        return abs(sym.det(self.linking_matrix()))
    
    def spinc_to_HF(self):
        return hf.spinc_to_HF(self.params)
    
    def print_HF(self):
        print('HF^+({}):'.format(self))
        hf.print_HF(self.params)
        return
    
    def correction_terms(self):
        return tuple(sorted(hf.correction_terms(self.params)))
    
    def is_lspace(self):
        return hf.is_lspace(self.params)
    
    def casson_walker(self):
        if self.is_lens_space():
            self = self.to_lens_space()
            return cw.casson_walker_lens(self.p, self.q)
        if self.is_prism_mfld():
            self = self.to_prism_mfld()
            return cw.casson_walker_prism(self.p, self.q)
        return cw.casson_walker(self.params)

    def is_lens_space(self):
        return self.exceptional_fibers <= 2
    
    def to_lens_space(self):
        if not self.is_lens_space():
            raise Exception('The SFS specified is not homeomorphic to a lens space!')
        coeffs = [sym.Rational(self.central_weight)]
        if self.exceptional_fibers >= 1:
            coeffs += cf.number_to_neg_cont_frac(self.branch_weights[0].p, self.branch_weights[0].q)
        if self.exceptional_fibers == 2:
            coeffs.reverse()
            coeffs += cf.number_to_neg_cont_frac(self.branch_weights[1].p, self.branch_weights[1].q)
        frac = -cf.number_from_neg_cont_frac(*coeffs)
        return Lens(frac.p, frac.q)
    
    def is_prism_mfld(self):
        return self.fractional_branch_weights.count(sym.Rational(1,2)) >= 2 and self.exceptional_fibers == 3
    
    def to_prism_mfld(self):
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
            return tuple(cf.number_to_neg_cont_frac(-self.p, self.q))
        return tuple(cf.number_to_neg_cont_frac(self.p, self.p - self.q))

def lens_from_linear_lattice(*params):
    p, q = cf.number_from_neg_cont_frac(*params).p, cf.number_from_neg_cont_frac(*params).q
    return -Lens(p,q)

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