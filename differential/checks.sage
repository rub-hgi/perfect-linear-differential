from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.crypto.sboxes import SBox

from ciphers.cipher import AESLikeCipher


def get_linear_layers(cipher):
    """
    Returns the linear layer, its inverse and the number of s-boxes for
    either two rounds, or a superbox if the cipher is AES-like
    :param cipher: Cipher to get the linear layer and its inverse from
    :return: Tuple of linear layer, its inverse and the number of s-boxes
    """
    if isinstance(cipher, AESLikeCipher):
        return cipher.mc_binary_matrix, \
            cipher.mc_inverse_binary_matrix, \
            cipher.nbr_sboxes // cipher.nbr_superboxes
    return cipher.L, \
        cipher.L_inverse, \
        cipher.nbr_sboxes


def filter_output_differences(cipher):
    """
    Use the main theorem to filter possible output differences for a probability one differential
    :param cipher: The cipher for which we want to filter
    :return: Output difference candidates that could possibly lead to a probability one differential over two rounds
    """
    S = cipher.S
    S_inv = cipher.S_inverse
    m = S.input_size()
    L, L_inv, d = get_linear_layers(cipher)
    BCT = S.boomerang_connectivity_table()
    LS = S_inv.linear_structures()
    beta_candidates = [{i for i in range(2**m)} for _ in range(d)]

    def leading_to_max_bct_entries(subspace):
        beta_i_candidates = set()
        for delta in subspace:
            if delta.is_zero():
                continue
            beta_i_candidates |= {beta_i for beta_i, entry in enumerate(BCT[S.from_bits(delta)]) if entry == 2 ** m}
        return beta_i_candidates

    def leading_to_linear_structures(subspace):
        beta_i_candidates = {0}
        for gamma in subspace:
            if gamma.is_zero():
                continue
            beta_i_candidates |= {beta_i for gamma_, beta_i, c in LS if gamma_ == S.from_bits(gamma)}
        return beta_i_candidates

    # For each pair i,j, use filter out output differences that cannot lead to a probability one differential
    for i in range(d):
        for j in range(d):
            L_ij = L[i*m:(i+1)*m, j*m:(j+1)*m]
            if L_ij.rank() == 0:
                continue  # L_ij = 0 which does not gives us anything
            L_inv_j_neq_i = L_inv[j*m:(j+1)*m, [l for l in range(m*d) if l not in range(i*m,(i+1)*m)]]
            rank = (L_ij * L_inv_j_neq_i).rank()
            if rank == 0:  # Case 1 of the main theorem
                if L_ij.rank() > 2:
                    continue
                # Filter candidates for beta_i based on boomerang connectivity table
                beta_candidates[i] &= leading_to_max_bct_entries(L_ij.column_space())
                # Filter candidates for beta_i based on linear structures
                beta_candidates[i] &= leading_to_linear_structures(
                    L[i*m:(i+1)*m, [l for l in range(m*d) if l not in range(j*m, (j+1)*m)]].left_kernel()
                )
            elif rank == m:  # Case 3 of the main theorem
                # Update candidates using all beta_i such that S^{-1}(S(x) + beta_i) is affine
                beta_candidates[i] &= beta_leading_to_affine(S)
            else:  # 0 < rank < m, i.e. Case 2 of the main theorem
                beta_candidates[i] &= leading_to_max_bct_entries(GF(2)**m)
                for k in range(d):
                    if k == i:
                        continue
                    if (L_ij * L[j*m:(j+1)*m, k*m:(k+1)*m]).rank() > 0:
                        beta_candidates[k] &= leading_to_linear_structures(GF(2)**m)
            if len(beta_candidates[i]) == 1:  # Only 0 is left as a candidate
                if not 0 in beta_candidates[i]:  # Sanity check
                    raise Exception(
                        "Zero is not a possible output difference? Something must have gone terribly wrong!"
                    )
                break
    return beta_candidates

def beta_leading_to_affine(S, check_existence=False):
    """
    Finds all beta such that S^{-1}(S(x) + beta) is affine in x
    :param S: SBox
    :param check_existence: Whether only the existence of a non-zero beta should be checked
    :return: List of beta such that S^{-1}(S(x) + beta) is affine in x
    """
    betas = {0}
    m = S.input_size()
    BCT = S.boomerang_connectivity_table()
    S_inv = S.inverse()

    for beta in range(1, 2 ** m):
        # Check if S^{-1}(S(x) + beta_i) is affine in x
        is_affine = True
        # As this can be quite expensive, use Lemma 10 (2) as a pre-filter
        for entry in BCT[:, beta]:
            if 0 < entry[0] < 2 ** m:
                is_affine = False
                break
        if not is_affine:
            continue
        # Explicitly calculate the algebraic degree of S^{-1}(S(x) + beta_i)
        F = SBox([S_inv[S[l] ^^ beta] for l in range(2 ** m)])
        for l in range(m):
            if F.component_function(1 << l).algebraic_degree() > 1:
                is_affine = False
                break
        if is_affine:
            if check_existence:
                return True
            betas.add(beta)
    if check_existence:
        return False
    return betas

def check_corollaries(cipher):
    """
    Checks Corollaries 7 to p
    :param cipher: Cipher to be checked
    :return: Tuple indicating for each Corollary whether it rules out
    the existence of a probability one differential over two rounds
    """
    L, _, d = get_linear_layers(cipher)
    cached_is_bn_2 = is_branch_number_2(L, cipher.S.input_size(), d)
    return check_corollary_7(cipher, cached_is_bn_2), \
        check_corollary_8(cipher, cached_is_bn_2), \
        check_corollary_9(cipher)

def check_corollary_7(cipher, cached_is_bn_2=None):
    """
    Checks Corollary 7
    :param cipher: Cipher to be checked
    :param cached_is_bn_2: Cached result of is_branch_number_2 call
    :return: Whether the corollary can rule out the existence of a probability one differential over two rounds
    """
    m = cipher.S.input_size()
    if cached_is_bn_2 is None:
        cached_is_bn_2 = is_branch_number_2(cipher.L, m, cipher.nbr_sboxes)
    if cached_is_bn_2:
        return False  # Corollary 9 is not applicable for branch number 2

    L, _, d = get_linear_layers(cipher)
    for i in range(d):  # For any i, try to find a block L_ij of full rank
        found = False
        for j in range(d):
            if L[i*m:(i+1)*m, j*m:(j+1)*m].rank() == m:
                found = True
                break
        if not found:  # Corollary 9 not applicable
            return False
    # Check if S^{-1}(S(x) + beta_i) is affine for some beta_i
    return beta_leading_to_affine(cipher.S, check_existence=True)

def check_corollary_8(cipher, cached_is_bn_2=None):
    """
    Checks Corollary 8
    :param cipher: Cipher to be checked
    :param cached_is_bn_2: Cached result of is_branch_number_2 call
    :return: Whether the corollary can rule out the existence of a probability one differential over two rounds
    """
    m = cipher.S.input_size()
    if cached_is_bn_2 is None:
        cached_is_bn_2 = is_branch_number_2(cipher.L, m, cipher.nbr_sboxes)
    if cached_is_bn_2:
        return False  # Corollary 10 is not applicable for branch number 2
    return cipher.S.boomerang_uniformity() != 2**m or len(cipher.S.inverse().linear_structures()) == 0

def check_corollary_9(cipher):
    """
    Checks Corollary 9
    :param cipher: Cipher to be checked
    :return: Whether the corollary can rule out the existence of a probability one differential over two rounds
    """
    m = cipher.S.input_size()
    if not is_permutation_matrix(cipher.L):
        return False  # Corollary 11 is only applicable for bit-permutations
    if cipher.S.differential_uniformity() == 2**m:
        return False  # Corollary 11 is only applicable if differential uniformity is not maximal
    # Check if all blocks L_ij have at most rank 1
    L, _, d = get_linear_layers(cipher)
    for i in range(d):
        for j in range(d):
            if L[i*m:(i+1)*m, j*m:(j+1)*m].rank() > 1:
                return False
    # Cipher passed Corollary 11
    return True

def is_branch_number_2(L, m, d):
    """
    Naive algorithm for checking is L has (differential) branch number 2 (with respect to the s-boxes)
    :param L: Linear layer as matrix over GF(2)
    :param m: S-box dimension
    :param d: Number of s-boxes
    :return: Whether the (differential) branch number of L (with respect to the s-boxes) is 2
    """
    for i in range(d):
        for v in L[i*m:(i+1)*m].image():
            activity_vector = vector(GF(2), [0 if v[m*l:m*(l+1)].is_zero() else 1 for l in range(d)])
            if activity_vector.hamming_weight() == 1:
                return True
    return False

def is_permutation_matrix(M):
    """
    :param M: Matrix over GF(2)
    :return: Whether M is a bit-permutation
    """
    for r in M:
        if r.hamming_weight() != 1:
            return False
    for c in M.transpose():
        if c.hamming_weight() != 1:
            return False
    return True

