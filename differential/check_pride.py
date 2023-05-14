from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF

from ciphers.pride import Pride


def check_pride():
    cipher = Pride()
    # Define matrices A and B as in section B.2
    A = cipher.L[0:4,0:4]
    B = cipher.L[32:36, 32:36]
    # make sure that the primary diagonal only consists of A and B
    for i in range(cipher.nbr_sboxes // 2):
        assert cipher.L[i*4:(i+1)*4, i*4:(i+1)*4] == A
        j = i + (cipher.nbr_sboxes // 2)
        assert cipher.L[j*4:(j+1)*4, j*4:(j+1)*4] == B

    S = cipher.S
    S_inv = cipher.S_inverse
    m = cipher.S.input_size()
    F2_m = GF(2)**m

    # Try to find a compatible non-zero output difference for some weak keys
    for beta in {1,8}:  # We know that the non-zero output differences per s-box are either 0x1 or 0x8
        vec_b = vector(GF(2), S.to_bits(beta))
        for M in [A,B]:
            # For any choice of x, generate right side of the equations we have to check
            right_sides = [S_inv[S[M * x] + vec_b] + S_inv[S[0] + vec_b] for x in F2_m]
            # Try to find a weak key and a (non-zero) alpha such that the equation holds for all x
            for key in F2_m:
                for alpha in {1,8}:
                    vec_a = vector(GF(2), S.to_bits(alpha))
                    const_a = S[S_inv[key] + vec_a]
                    compatible = True
                    for x, y in zip(F2_m, right_sides):
                        if M* (S[S_inv[x + key] + vec_a] + const_a) != y:
                            compatible = False
                            break
                    if compatible:
                        # While we are technically not done here, as there are other equations that could still violate
                        # the compatibility of alpha and beta, this condition is never met for Pride.
                        # Hence, we can stop here
                        print(f"Found a potential probability on differential (over two rounds) for Pride")
                        print(f"alpha={alpha}, beta={beta}, key={key} are compatible for matrix {'A' if M == A else 'B'}")
                        return True
    print(f"No probability one differential (over two rounds) exists for Pride")
    return False
