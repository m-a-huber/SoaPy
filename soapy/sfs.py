import sympy as sym
from sympy.matrices.normalforms import invariant_factors
from soapy import hf_nemethi as hf
from soapy import casson_walker as cw
from soapy import continued_fractions as cf
from soapy import ribbon_obstruction as rib


class SFS:
    """This is a class for working with orientable Seifert fibered spaces (SFS)
    whose base orbifold is the 2-sphere.

    Attributes:
        params (list[int]): List of integer coefficients representing the SFS
            specified. These do not necessarily coincide with the input
            parameters, but rather are normalized in such a way that the
            corresponding integer surgery diagram is definite.
        central_weight (int): The weight of the central vertex of the
            normalized surgery description of the SFS specified.
        branch_weights (tuple(sym.Rational)): Tuple containing the rational
            surgery coefficients of the exceptional fibers of the SFS
            specified.
    """

    def __init__(self, *params):
        """The constructor for the SFS class.

        Args:
            *params (list[int]): A variable number of integers specifying a
                SFS. The list must be of the format
                (e, a_1, b_1, ..., a_n, b_n), where e is the integer framing
                of the central unknot, and the a_i/b_i are the rational
                framings of the exceptional fibers, each expressed as a
                reduced fraction.

        Raises:
            ValueError: If the input list is not of the required format.
            ValueError: If the SFS specified has positive first Betti number.
        """
        cond_1 = all(isinstance(el, int) for el in params)
        cond_2 = len(params) % 2 == 1
        cond_3 = all([
            sym.gcd(params[i], params[i+1]) == 1
            for i in range(1, len(params), 2)
        ])
        cond_4 = all(el != 0 for el in params[1:])
        if not min(cond_1, cond_2, cond_3, cond_4):
            raise ValueError(
                "A SFS should be specified by an odd number of integers, with "
                "each exceptional fiber specified by a pair of non-zero "
                "coprime integers!"
            )
        cond_5 = hf.normalize(params)[0] != 0
        if not cond_5:
            raise ValueError(
                "The SFS specified is not a rational homology sphere!"
            )
        self.params = tuple(hf.normalize(params)[1])
        self.central_weight = self.params[0]
        self.branch_weights = tuple(
            sym.Rational(self.params[i], self.params[i+1])
            for i in range(1, len(self.params), 2)
        )

    @classmethod
    def from_plumbing(cls, central_weight, *lists_of_coeffs):
        """Allows one to construct a soapy.SFS object from an integer plumbing
        description, all of whose weights are non-zero.

        Args:
            central_weight (int): The weight of the central vertex of the
                plumbing tree.
            *lists_of_coeffs (list[int]): A variable number of lists of
                weights of the branches, each read starting from the central
                vertex.

        Raises:
            ValueError: If any of the weights specified is zero.

        Returns:
            soapy.SFS: The SFS corresponding to the integer plumbing
                description specified.
        """
        cond_1 = isinstance(central_weight, int)
        cond_2 = all(isinstance(lst, list) for lst in lists_of_coeffs)
        cond_3 = all(
            all(isinstance(el, int) for el in lst)
            for lst in lists_of_coeffs
        )
        cond_4 = all(
            all(el != 0 for el in lst)
            for lst in lists_of_coeffs
        )
        if not min(cond_1, cond_2, cond_3, cond_4):
            raise ValueError(
                "The weights on the plumbing graph should be non-zero "
                "integers!"
            )
        branch_weights = tuple(
            cf.number_from_neg_cont_frac(*[i for i in lst])
            for lst in lists_of_coeffs
        )
        params = sum(
            ((w.p, w.q) for w in branch_weights),
            ()
        )
        return cls(central_weight, *params)

    def __repr__(self):
        branch_weights_string = ", ".join([
            str(q)
            for q in self.branch_weights
        ])
        return f"Y({self.central_weight}; {branch_weights_string})"

    def __neg__(self):
        """Returns the SFS with reversed orientation.
        """
        params = tuple(
            -((-1)**i)*self.params[i]
            for i in range(len(self.params))
        )
        return SFS(*params)

    def __eq__(self, other):
        """Checks if the two SFS are orientation-preservingly homeomorphic.
        """
        if isinstance(other, SFS):
            if not max(self.is_lens_space(), other.is_lens_space()):
                seifert_self = (
                    self.euler_number(),
                    sorted(self.fractional_branch_weights())
                )
                seifert_other = (
                    other.euler_number(),
                    sorted(other.fractional_branch_weights())
                )
                return seifert_self == seifert_other
            if self.is_lens_space() and other.is_lens_space():
                self, other = self.to_lens_space(), other.to_lens_space()
                return all([
                    self.p == other.p,
                    max(
                        self.q == other.q,
                        (self.q*other.q) % self.p == 1
                    )
                ])
        return False

    def __hash__(self):
        if not self.is_lens_space():
            seifert_tuple = (
                self.euler_number(),
                sorted(self.fractional_branch_weights())
            )
            return hash(seifert_tuple)
        else:
            L = self.to_lens_space()
            q = min(
                self.q,
                pow(self.q, -1, self.p)
            )
            return hash((L.p, q))

    def __le__(self, other):
        """Checks if the first SFS passes the d-invariant obstruction to
        admitting a ribbon rational homology cobordism to the second SFS.
        """
        if isinstance(other, SFS):
            if self == other:
                return True
            return rib.ribbon_obstruction_sfs(self.params, other.params)
        return False

    def __lt__(self, other):
        """Same as __le__, but returns False if the two SFS are, in fact,
        orientation-preservingly homeomorphic to one another.
        """
        if isinstance(other, SFS):
            if self == other:
                return False
            return self <= other
        return False

    def fractional_branch_weights(self):
        """Returns the fractional parts of the framings on the exceptional
        fibers.

        Returns:
            tuple(sym.Rational): Tuple containing the fractional parts of the
                branch weights.
        """
        return tuple(
            (1 / w) - sym.floor(1 / w)
            for w in self.branch_weights
        )

    def euler_number(self):
        """Returns the orbifold Euler number of the SFS specified.

        Returns:
            sym.Rational: The orbifold Euler number of the normalized surgery
                description of the SFS specified.
        """
        return self.central_weight - sum([1 / w for w in self.branch_weights])

    def seifert_invariants(self):
        """Returns the Seifert invariants of the SFS specified.

        Returns:
            tuple(sym.Rational, tuple(sym.Rational)): A tuple of the format
                (Euler number, (tuple of fractional branch weights)).
        """
        return (self.euler_number(), self.fractional_branch_weights())

    def number_of_exceptional_fibers(self):
        """Returns the number of exceptional fibers of the SFS specified.

        Returns:
            int: The number of exceptional fibers of the normalized surgery
                description of the SFS specified.
        """
        return len(self.branch_weights)

    def to_plumbing(self):
        """Returns the definite plumbing (equivalently: an integral surgery
        description) corresponding to the SFS specified.

        Returns:
            tuple(int, lists[int]): A tuple whose first elements is the
                central weight, followed by the lists of integer weights on
                the branches (read starting from the central vertex).
        """
        lists_of_coeffs = (
            cf.number_to_neg_cont_frac(w.p, w.q)
            for w in self.branch_weights
        )
        return ((self.central_weight,)
                + tuple(lst for lst in lists_of_coeffs))

    def linking_matrix(self):
        """Returns the linking matrix of the SFS specified.

        Returns:
            sym.Matrix: A SymPy-matrix with SymPy-integers as entries,
                representing the linking matrix of the integer plumbing
                corresponding to the SFS specified.
        """
        return hf.linking_matrix(self.params)

    def first_homology(self):
        """Returns the first homology of the SFS specified.

        Returns:
            tuple(int): The orders of the non-trivial cyclic summands of the
                first homology of the SFS specified.
        """
        return tuple(
            int(i)
            for i in map(abs, invariant_factors(self.linking_matrix()))
            if i != 1
        )

    def order_of_first_homology(self):
        """Returns the order of the first homology of the SFS specified.

        Returns:
            int: The order of the first homology of the SFS specified.
        """
        return int(abs(sym.det(self.linking_matrix())))

    def spinc_to_HF(self):
        """Computes HF^+ of the SFS specified in each spin^c-structure. The
        Z[U]-module-structure of HF^+ is encoded as a dictionary of the format
        {'order of Z[U]-module-summand': list of bottommost gradings of all
        Z[U]-module-summands of that order'}.

        Returns:
            dict: A dictionary of the format {'spin^c-structure':
                'Z[U]-module-structure of HF^+'}.
        """
        return hf.spinc_to_HF(self.params)

    def print_HF(self):
        """Prints HF^+ of the SFS specified in a more legible manner.

        Returns:
            None: Just prints HF^+ of the SFS specified.
        """
        print(f"HF^+({self}) = ")
        hf.print_HF(self.params)
        return

    def correction_terms(self):
        """Returns a list of the corrections terms of the SFS specified.

        Returns:
            list[sym.Rational]: List of all correction terms of the SFS
                specified.
        """
        return tuple(sorted(hf.correction_terms(self.params)))

    def is_lspace(self):
        """Checks whether or not the SFS specified is a Heegaard Floer L-space.

        Returns:
            bool: Whether or not the the SFS specified is an L-space.
        """
        return hf.is_lspace(self.params)

    def casson_walker(self):
        """Computes the Casson-Walker invariant of the SFS specified.

        Returns:
            sym.Rational: The Casson-Walker invariant of the SFS specified.
        """
        if self.is_lens_space():
            self = self.to_lens_space()
            return cw.casson_walker_lens(self.p, self.q)
        if self.is_prism_mfld():
            self = self.to_prism_mfld()
            return cw.casson_walker_prism(self.p, self.q)
        return cw.casson_walker(self.params)

    def is_lens_space(self):
        """Checks whether or not the SFS specified is homeomorphic to a lens
        space.

        Returns:
            bool: Whether or not the the SFS specified is homeomorphic to a
                lens space.
        """
        return self.number_of_exceptional_fibers() <= 2

    def to_lens_space(self):
        """Transforms the SFS specified into the corresponding lens space.

        Raises:
            ValueError: If the SFS specified is not homeomorphic to any lens
                space.

        Returns:
            soapy.Lens: Lens space homeomorphic to the SFS specified.
        """
        if not self.is_lens_space():
            raise ValueError(
                "The SFS specified is not homeomorphic to a lens space!"
            )
        coeffs = [sym.Rational(self.central_weight)]
        if self.number_of_exceptional_fibers() >= 1:
            coeffs += cf.number_to_neg_cont_frac(
                self.branch_weights[0].p,
                self.branch_weights[0].q
            )
        if self.number_of_exceptional_fibers() == 2:
            coeffs.reverse()
            coeffs += cf.number_to_neg_cont_frac(
                self.branch_weights[1].p,
                self.branch_weights[1].q
            )
        frac = -cf.number_from_neg_cont_frac(*coeffs)
        return Lens(frac.p, frac.q)

    def is_prism_mfld(self):
        """Checks whether or not the SFS specified is homeomorphic to a prism
        manifold.

        Returns:
            bool: Whether or not the the SFS specified is homeomorphic to a
                prism manifold.
        """
        return (self.fractional_branch_weights().count(sym.Rational(1, 2)) >= 2
                and self.number_of_exceptional_fibers() == 3)

    def to_prism_mfld(self):
        """Transforms the SFS specified into the corresponding prism manifold.

        Raises:
            ValueError: If the SFS specified is not homeomorphic to any prism
                manifold.

        Returns:
            soapy.Prism: Prism manifold homeomorphic to the SFS specified.
        """
        if not self.is_prism_mfld():
            raise ValueError(
                "The SFS specified is not homeomorphic to a prism manifold!"
            )
        return Prism((self.euler_number()).q, (self.euler_number()).p)

    # Define default variables that represent S^3 and
    # the Poincaré homology sphere, respectively.


three_sphere = SFS(1)
poincare_sphere = SFS(-2, -5, 4, -3, 2, -2, 1)


class Lens(SFS):
    """This is a subclass of SFS representing lens spaces.

    Attributes:
        p (int): The first parameter of the lens space specified, normalized
            to be greater than zero.
        q (int): The second parameter of the lens space specified, normalized
            so that p > q > 0.
    """

    def __init__(self, p, q):
        """The constructor for the Lens subclass. A lens space is specified by
        a pair of non-zero coprime integers p and q. The constructor
        normalizes the parameters to satisfy p > q > 0.

        Args:
            p (int): The first parameter of the lens space.
            q (int): The second parameter of the lens space.

        Raises:
            ValueError: If the parameters are not a pair of non-zero coprime
                integers.
        """
        cond_1 = all(isinstance(el, int) for el in (p, q))
        cond_2 = sym.gcd(p, q) == 1
        cond_3 = p*q != 0
        if not min(cond_1, cond_2, cond_3):
            raise ValueError(
                "A lens space should be specified by a pair of non-zero "
                "coprime integers!"
            )
        super().__init__(0, q, p)
        self.p = abs(p)
        self.q = (lambda x: 1 if x == 0 else x)(
            (lambda q: q % p if p > 0 else (self.p-q) % self.p)(q)
        )

    @classmethod
    def from_linear_lattice(cls, *params):
        """Allows one to construct a soapy.Lens object from a linear lattice
        specifying a lens space, all of whose weights are non-zero.

        Args:
            *params (int): A variable number of integer weights of the linear
                lattice, read starting from either end.

        Raises:
            ValueError: If any of the weights specified is zero.

        Returns:
            soapy.Lens: The lens space corresponding to the linear lattice
                specified.
        """
        cond_1 = all(isinstance(el, int) for el in params)
        cond_2 = all(el != 0 for el in params)
        if not min(cond_1, cond_2):
            raise ValueError(
                "The weights on the linear lattice should be non-zero "
                "integers!"
            )
        p, q = (
            cf.number_from_neg_cont_frac(*params).p,
            cf.number_from_neg_cont_frac(*params).q
        )
        return -cls(p, q)

    def __repr__(self):
        return f"L({self.p},{self.q})"

    def __neg__(self):
        """Returns the lens space with reversed orientation.
        """
        return Lens(self.p, -self.q)

    def to_SFS(self):
        """Transforms the lens space specified into a SFS.

        Returns:
            soapy.SFS: The SFS homeomorphic to the lens space specified.
        """
        return SFS(0, self.q, self.p)

    def to_linear_lattice(self, epsilon=-1):
        """Returns the weights of the linear lattice bounded by the lens space
        specified. By default, the negative definite linear lattice is
        returned, unless epsilon is set to 1.

        Args:
            epsilon (int, optional): The sign of definiteness of the linear
                lattice to be returned. Defaults to -1.

        Raises:
            ValueError: If epsilon is not 1 in absolute value.

        Returns:
            tuple(int): A tuple containing the weights of the linear lattice
            bounded by the lens space specified.
        """
        if epsilon not in {-1, 1}:
            raise ValueError(
                "The second (optional) argument should be plus/minus 1!"
            )
        if epsilon == -1:
            return tuple(cf.number_to_neg_cont_frac(-self.p, self.q))
        return tuple(cf.number_to_neg_cont_frac(self.p, self.p - self.q))


class Prism(SFS):
    """This is a subclass of SFS representing prism manifolds.

    Attributes:
        p (int): The first parameter of the prism manifold specified,
            normalized to be greater than 0.
        q (int): The second parameter of the prism manifold specified.
    """

    def __init__(self, p, q):
        """The constructor for the Prism subclass. A prism manifold is
        specified by a pair of non-zero coprime integers p and q. Note,
        however, that P(p,q) is, in fact, homeomorphic to a lens space if p
        equals 1 in absolute value.

        Args:
            p (int): The first parameter of the prism manifold.
            q (int): The second parameter of the prism manifold.

        Raises:
            ValueError: If the parameters are not a pair of non-zero coprime
                integers p and q, where p is greater than 1.
        """
        cond_1 = all(isinstance(el, int) for el in (p, q))
        cond_3 = sym.gcd(p, q) == 1
        if not min(cond_1, cond_3):
            raise ValueError(
                "A prism manifold should be specified by a pair of coprime "
                "integers!"
            )
        super().__init__(-1, -2, 1, -2, 1, -p, q)
        self.p = abs(p)
        self.q = int(sym.sign(p)*q)

    def __repr__(self):
        return f"P({self.p},{self.q})"

    def __neg__(self):
        """Returns the prism manifold with reversed orientation.
        """
        return Prism(self.p, -self.q)

    def to_SFS(self):
        """Transforms the prism manifold specified into a SFS.

        Returns:
            soapy.SFS: The SFS homeomorphic to the prism manifold specified.
        """
        return SFS(-1, -2, 1, -2, 1, -self.p, self.q)


class Brieskorn(SFS):
    """This is a subclass of SFS representing Brieskorn homology spheres.

    Attributes:
        coeffs (list[int]): The coefficients of the Brieskorn homology sphere
            specified.
    """

    def __init__(self, *coeffs):
        """The constructor for the Brieskorn subclass. A Brieskorn homology
        sphere is specified by a list of integers a_1, .., a_n that are
        pairwise coprime, each greater than 1 in absolute value, and all of
        the same sign. If all the coefficients are negative, this returns
        -Sigma(a_1, ...,a_n) (i.e. with orientation reversed).

        Args:
            *coeffs (int): The coefficients specifying a Brieskorn homology
                sphere.

        Raises:
            ValueError: If the coefficients are not a list of integers
                a_1, ..., a_n that are pairwise coprime, each greater than 1
                in absolute value, and all of the same sign.
        """
        cond_1 = all(isinstance(el, int) for el in coeffs)
        cond_2 = min(map(abs, coeffs)) > 1
        cond_3 = len(set(map(sym.sign, coeffs))) == 1
        cond_4 = sym.lcm_list(coeffs) == sym.prod(coeffs)
        if not min(cond_1, cond_2, cond_3, cond_4):
            raise ValueError(
                "A Brieskorn homology sphere should be specified by a list of "
                "pairwise coprime integers, each greater than 1 in absolute "
                "value, and all of the same sign!"
            )
        super().__init__(*hf.brieskorn(*coeffs))
        self.coeffs = coeffs

    def __repr__(self):
        epsilon = sym.sign(self.coeffs[0])
        sign = (lambda x: "" if x > 0 else "-")(epsilon)
        coeffs_string = ",".join([str(abs(a)) for a in self.coeffs])
        return f"{sign}Sigma({coeffs_string})"

    def __neg__(self):
        """Returns the Brieskorn homology sphere with reversed orientation.
        """
        return Brieskorn(*map(lambda x: -x, self.coeffs))

    def to_SFS(self):
        """Transforms the Brieskorn homology sphere specified into a SFS.

        Returns:
            soapy.SFS: The SFS homeomorphic to the Brieskorn homology sphere
                specified.
        """
        return SFS(*hf.brieskorn(*self.coeffs))
