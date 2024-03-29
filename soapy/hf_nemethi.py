import sympy as sym
from soapy import continued_fractions as cf

# Auxiliary functions that convert 3-manifold
# notation into the format required for input.


def normalize(C):
    """This returns the Euler number as well as the normalized input string
    (corresponding to a definite plumbing) for the SFS specified. If the SFS
    specified has vanishing Euler number, an empty list of coefficients is
    returned.

    Args:
        C (list of int): List of integer coefficients specifying a SFS.

    Returns:
        tuple(sym.Rational, list[int]): A tuple containing the Euler number
            and the coefficients of the input string of SFS specified,
            expressed as a definite plumbing.
    """
    central_weight = sym.Rational(C[0])
    branch_weights = tuple(
        sym.Rational(C[i], C[i+1])
        for i in range(1, len(C), 2)
    )
    fractional_branch_weights = [
        (1 / w) - sym.floor(1 / w)
        for w in branch_weights
    ]
    euler_number = central_weight - sum([1 / w for w in branch_weights])
    if euler_number == 0:
        return (0, [])
    new_branch_weights = [
        2 / (2 * w + sym.sign(euler_number) - 1)
        for w in fractional_branch_weights
        if w != 0
    ]
    new_central_weight = euler_number + sum([
        1 / w
        for w in new_branch_weights
    ])
    return (
        euler_number,
        ([int(new_central_weight)]
         + [i for w in new_branch_weights for i in [w.p, w.q]])
    )


# def lens(p, q):
#     """This returns the input string (corresponding to the negative definite
#     plumbing) for the lens space L(p, q), for p and q coprime.

#     Args:
#         p (int): The first parameter of L(p, q), a non-zero integer.
#         q (int): The second parameter of L(p, q), a non-zero integer.

#     Returns:
#         list[int]: The coefficients of the input string for the lens space
#             L(p, q).
#     """
#     if abs(p) == 1:
#         return [-1]
#     epsilon = sym.sign(p)
#     p, q = abs(p), q % abs(p)
#     if epsilon == -1:
#         q = p-q
#     if abs(q) == 1:
#         return list(map(int, [-sym.Rational(p, q)]))
#     return list(map(int, [
#         -sym.ceiling(sym.Rational(p, q)),
#         -q,
#         sym.ceiling(sym.Rational(p, q)) * q - p
#     ]))


# def prism(p, q):
#     """This returns the input string for the prism manifold P(p, q).

#     Args:
#         p (int): The first parameter of P(p, q), an integer greater than 1 in
#             absolute value.
#         q (int): The second parameter of P(p, q), a non-zero integer.

#     Returns:
#         list[int]: The coefficients of the input string for the prism
#             manifold P(p, q).
#     """
#     if p <= 0:
#         p, q = -p, -q
#     epsilon = -1
#     if q > 0:
#         epsilon, q = 1, -q
#     r = -(1 - sym.Rational(q, p))
#     c = -sym.ceiling(-r)
#     if p == 1:
#         L = [c, -2, 1, -2, 1]
#     else:
#         x = sym.Rational(p, p*c+p-q)
#         L = [c, -2, 1, -2, 1, x.p, x.q]
#     L = [((-epsilon)**i)*L[i] for i in range(len(L))]
#     L[0] = -epsilon*L[0]
#     return list(map(int, L))


def brieskorn(*A):
    """This returns the input string for the Brieskorn homology sphere
    Sigma(a_1, ..., a_n). If all of a_1, ..., a_n are negative, the input
    string for -Sigma(a_1, ..., a_n) is returned.

    Args:
        *args (int): The coefficients of Sigma(a_1, ..., a_n), each greater
            than 1 in abolute value and all of the same sign.

    Returns:
        list[int]: The coefficients of the input string for the Brieskorn
            homology sphere Sigma(a_1, ..., a_n).
    """
    epsilon = sym.sign(A[0])
    if epsilon == -1:
        A = list(map(abs, A))
    a_0 = sym.prod(A)
    A_hat = [
        sym.floor(sym.Rational(a_0, A[i]))
        for i in range(len(A))
    ]
    B = [
        sym.mod_inverse(-A_hat[i], A[i])
        for i in range(len(A))
    ]
    e0 = sym.floor(
        sym.Rational(1, a_0) * (-1 - sum([
            A_hat[i]*B[i]
            for i in range(len(A))
        ]))
    )
    L = [
        e0 if j == 0 else
        -A[sym.floor(sym.Rational(j, 2))] if j % 2 == 1 else
        B[sym.floor(sym.Rational(j, 2))-1] if j % 2 == 0 else
        "void"
        for j in range(2*len(A)+1)
    ]
    if epsilon == -1:
        L[0] = - L[0]
        for i in range(1, len(L), 2):
            L[i] = -L[i]
    return list(map(int, L))

# Functions that execute Némethi's algorithm.


def linking_matrix(C):
    L = [
        len(cf.number_to_neg_cont_frac(C[i], C[i+1]))
        for i in range(1, len(C), 2)
    ]
    Ell = sym.zeros(1+sum(L), 1+sum(L))
    Ell[0, 0] = C[0]
    for j in range(len(L)):
        for k in range(L[j]):
            if k == 0:
                Ell[0, 1+sum(L[:j])] = 1
                Ell[1+sum(L[:j]), 0] = 1
            else:
                Ell[k+sum(L[:j]), 1+k+sum(L[:j])] = 1
                Ell[1+k+sum(L[:j]), k+sum(L[:j])] = 1
            Ell[1+k+sum(L[:j]), 1+k+sum(L[:j])] = cf.number_to_neg_cont_frac(
                C[2*j+1], C[2*(j+1)]
            )[k]
    return Ell


def add_one(x, y):
    for i in range(len(x)):
        if x[i] < y[i]:
            x[i] += 1
            return x
        elif x[i] == y[i]:
            x[i] = 0
    return x


def checkall(a, a_0, e, e_0, Alpha, Beta):
    n = sym.ceiling(sym.Rational(1+a_0+len(Alpha), abs(e)))+1
    for i in range(1, n):
        cond = (1
                + a_0
                + (i*e_0)
                + sum([
                    sym.floor(
                        sym.Rational(i * Beta[ix] + a[ix], Alpha[ix])
                    )
                    for ix in range(len(Alpha))
                ])) > 0
        if cond:
            return False
    return True


def truncated_spinc(e, e_0, Alpha, Beta):
    S = []
    for j in range(len(Alpha)):
        for i in reversed(range(Alpha[j])):
            a = [0]*len(Beta)
            a.insert(j, i)
            if checkall(a, 0, e, e_0, Alpha, Beta):
                S.append(max(1, i))
                break
    return S


def spinc(e, e_0, Alpha, Beta, n):
    S = truncated_spinc(e, e_0, Alpha, Beta)
    p = sym.prod([i + 1 for i in S])
    a = [0]*len(Alpha)
    A = {tuple(a)}
    for i in range(p):
        a = add_one(a, S)
        A.add(tuple(a))
    B = set()
    for a in A:
        i = 0
        while len(B) < n and checkall(a, i, e, e_0, Alpha, Beta):
            b = (i,) + a
            B.add(tuple(b))
            i = i+1
    return B


def kill_repetitions(tau):
    tau_norep = [tau[0]]
    for i in range(1, len(tau)):
        if not tau[i] == tau[i-1]:
            tau_norep.append(tau[i])
    return tau_norep


def reduced_sequence(tau):
    tau = kill_repetitions(tau)
    red_tau = [tau[0]]
    for i in range(1, len(tau)):
        if i == len(tau)-1:
            red_tau.append(tau[i])
        elif tau[i] >= tau[i-1] and tau[i] >= tau[i+1]:
            red_tau.append(tau[i])
        elif tau[i] <= tau[i-1] and tau[i] <= tau[i+1]:
            red_tau.append(tau[i])
    return red_tau


def tau_delta(a, e, e_0, Alpha, Beta):
    P = [0]
    D = [
        1
        + a[0]
        + i*abs(e_0)
        + sum([
            sym.floor(sym.Rational(-i * Beta[j] + a[j+1], Alpha[j]))
            for j in range(len(Alpha))
        ])
        for i in range(sym.ceiling(
            sym.Rational(len(Alpha), abs(e))) + 1
        )
    ]
    for i in range(len(D)):
        P.append(P[-1] + D[i])
    return P


def shift(a, e_0, Alpha, Beta):
    e = e_0 + sum([
        sym.Rational(Beta[i], Alpha[i])
        for i in range(len(Alpha))
    ])
    omega = [
        sym.mod_inverse(Beta[i], Alpha[i])
        for i in range(len(Alpha))
    ]
    epsilon = sym.Rational(
        2 - len(Alpha) + sum([
            sym.Rational(1, Alpha[i])
            for i in range(len(Alpha))
        ]),
        e
    )
    a_tilde = a[0] + sum([
        sym.Rational(a[i+1], Alpha[i])
        for i in range(len(Alpha))
    ])
    x, y, z, v = (
        sym.Rational(sum(a), 2),
        sym.Rational(epsilon*a_tilde, 2),
        sym.Rational(a_tilde**2, 2*e),
        sym.Rational(0, 1)
    )
    for ix in range(len(Alpha)):
        v += sum([
            (sym.Rational(i * omega[ix], Alpha[ix])) % 1
            for i in range(a[ix+1]+1)
        ])
    return x + y + z - v


def vect_degree(v, Ell):
    M = (sym.Matrix(v).transpose()
         * sym.Matrix(Ell).inv()
         * sym.Matrix(v))
    return -(sym.Rational((M[0, 0]).p, (M[0, 0]).q) + len(v)) / 4


def degree_shift(S):
    K = [
        -(S[i, i] + 2)
        for i in range(sym.shape(S)[0])
    ]
    return vect_degree(K, S)

# Auxiliary functions that convert the graded
# roots and their correction terms into HF^+.


def tau_to_module(tau):
    """Computes the Z[U]-module corresponding to a reduced tau-sequence. In
    the resulting module, multiplication by U lowers the grading by 1.

    Args:
        tau (list[int]): A list of integers representing a reduced
            tau-sequence.

    Returns:
        dict: A dictionary of the format {'order of Z[U]-module-summand':
            'list of bottommost gradings of all Z[U]-module-summands of that
            order'}, encoding a Z[U]-module.
    """
    s = tau.index(min(tau))
    module = {0: [tau[s]]}
    for ix in reversed(range(s)):
        if tau[ix+1] - tau[ix] > 0:
            try:
                module[tau[ix+1] - tau[ix]].append(tau[ix])
            except IndexError:
                module[tau[ix+1] - tau[ix]] = [tau[ix]]
    for r in range(s+1, len(tau)):
        if tau[r-1] - tau[r] > 0:
            try:
                module[tau[r-1] - tau[r]].append(tau[r])
            except IndexError:
                module[tau[r-1] - tau[r]] = [tau[r]]
    return module


def module_corr_to_neg_HF(module, corr):
    """Computes the Z[U]-module corresponding to a given Z[U]-module endowed
    with a correction term. In the resulting module, multiplication by U
    lowers the grading by 2.

    Args:
        module (dict): A dictionary of the format {'order of
            Z[U]-module-summand': 'list of bottommost gradings of all
            Z[U]-module-summands of that order'}, encoding a Z[U]-module.
        corr (sym.Rational): The correction term corresponding to the module.

    Returns:
        dict: A dictionary of the format {'order of Z[U]-module-summand':
            'list of bottommost gradings of all Z[U]-module-summands of that
            order'}, encoding a Z[U]-module.
    """
    for key in module.keys():
        module[key] = [2 * i for i in module[key]]
    delta = corr - module[0][0]
    for key in module.keys():
        module[key] = [i + delta for i in module[key]]
    return module


def minus_HF(module):
    """Computes the Z[U]-module corresponding to HF^+(Y), given the
    Z[U]-module corresponding to HF^+(-Y).

    Args:
        module (dict): A dictionary of the format {'order of
            Z[U]-module-summand': 'list of bottommost gradings of all
            Z[U]-module-summands of that order'}, encoding a Z[U]-module.

    Returns:
        dict: A dictionary of the format {'order of Z[U]-module-summand':
            'list of bottommost gradings of all Z[U]-module-summands of that
            order'}, encoding a Z[U]-module.
    """
    for key in module.keys():
        module[key] = [-(i + key) for i in module[key]]
    return module

# Functions that return the correspondence between spin^c-structures and the
# reduced tau sequences, Z[U]-modules and HF^+, respectively, and the
# correction terms.


def spinc_to_tau_corr(C):
    """Computes the reduced tau-sequences and the correction terms
    corresponding to a spin^c-structure on the manifold specified by the
    plumbing (which has to be negative definite).

    Args:
        C (list[int]): The coefficients specifying a negative definite
            plumbing.

    Returns:
        dict: A dictionary of the format {'spin^c-structure':
        [reduced tau-sequence, correction term]}.
    """
    S = linking_matrix(C)
    n = abs(sym.det(S))
    Alpha = [
        abs(C[i])
        for i in range(1, len(C), 2)
    ]
    Beta = [
        abs(C[j])
        for j in range(2, len(C)+1, 2)
    ]
    e = C[0] + sum([
        sym.Rational(Beta[i], Alpha[i])
        for i in range(len(Beta))
    ])
    A = spinc(e, C[0], Alpha, Beta, n)
    def s(a): return reduced_sequence(tau_delta(list(a), e, C[0], Alpha, Beta))
    spinc_to_tau_corr = {
        a: [
            s(a)[:-1],
            2*min(s(a)) + degree_shift(S) - 2 * shift(a, C[0], Alpha, Beta)
        ]
        for a in A
    }
    return spinc_to_tau_corr


def spinc_to_module_corr(C):
    """Computes the Z[U]-modules (coming from a reduced tau-sequence) and the
    correction terms corresponding to a spin^c-structure on the manifold
    specified by the input.

    Args:
        C (list[int]): The coefficients specifying a negative definite
            plumbing.

    Returns:
        dict: A dictionary of the format {'spin^c-structure':
        [Z[U]-module, correction term]}.
    """
    spinc_to_module_corr = spinc_to_tau_corr(C)
    for key in spinc_to_module_corr.keys():
        spinc_to_module_corr[key] = [
            tau_to_module(spinc_to_module_corr[key][0]),
            spinc_to_module_corr[key][1]
        ]
    return spinc_to_module_corr


def spinc_to_HF(C):
    """Computes HF^+ in each spin^c-structure of the manifold specified by the
    input. The Z[U]-module-structure of HF^+ is encoded as a dictionary of the
    format {'order of Z[U]-module-summand': 'list of bottommost gradings of
    all Z[U]-module-summands of that order'}.

    Args:
        C (list[int]): The coefficients specifying a definite plumbing.

    Returns:
        dict: A dictionary of the format {'spin^c-structure':
        'Z[U]-module-structure of HF^+'}.
    """
    epsilon = sym.sign(C[0])
    if epsilon == 1:
        C = [-abs(i) for i in C]
        for i in range(2, len(C), 2):
            C[i] = abs(C[i])
    spinc_to_HF = spinc_to_module_corr(C)
    for key in spinc_to_HF.keys():
        spinc_to_HF[key] = module_corr_to_neg_HF(
            spinc_to_HF[key][0],
            spinc_to_HF[key][1]
        )
    if epsilon == -1:
        for key in spinc_to_HF.keys():
            spinc_to_HF[key] = minus_HF(spinc_to_HF[key])
    return spinc_to_HF

# Additional functions processing the output of spinc_to_HF.


def print_HF(C):
    """Prints HF^+ of the Seifert fibered space specified by a definite
    plumbing in a more legible manner.

    Args:
        C (list[int]): Coefficients specifying a definite plumbing.

    Returns:
        None: Prints HF^+ of the Seifert fibered space specified.
    """
    res = spinc_to_HF(C)
    for key in res.keys():
        if len(key) == 1:
            print(f"spin^c-structure: ({key[0]})")
        else:
            print(f"spin^c-structure: {key}")
        s = f"HF^+ = T_({res[key][0][0]})"
        for keykey in [i for i in res[key].keys() if i != 0]:
            for deg in res[key][keykey]:
                s += f" + Z({keykey})_({deg})"
        print(s)
    return


def correction_terms(C):
    """Returns a list of the corrections terms of the Seifert fibered space
    specified by a definite plumbing.

    Args:
        C (list[int]): Coefficients specifying a definite plumbing.

    Returns:
        list[sym.Rational]: List of all correction terms of the Seifert
            fibered space.
    """
    epsilon = sym.sign(C[0])
    if epsilon == 1:
        C = [-abs(i) for i in C]
        for i in range(2, len(C), 2):
            C[i] = abs(C[i])
    corr = spinc_to_tau_corr(C)
    return [
        epsilon*corr[key][1]
        for key in corr.keys()
    ]


def is_lspace(C):
    """Checks whether or not a Seifert fibered space specified by a definite
    plumbing is an L-space.

    Args:
        C (list[int]): Coefficients specifying a definite plumbing.

    Returns:
        bool: Whether or not the Seifert fibered space is an L-space.
    """
    res = spinc_to_HF(C)
    for key in res.keys():
        for keykey in res[key].keys():
            if keykey != 0:
                return False
    return True
