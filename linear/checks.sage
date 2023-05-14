from ciphers.cipher import AESLikeCipher

def linear_check(cipher, r=2):
    if (r == 2) and isinstance(cipher, AESLikeCipher):
        # only check one supberbox
        nbr_sboxes = cipher.nbr_sboxes // cipher.nbr_superboxes
        S = lambda state: cipher.sbox_layer_faster(state, nbr_sboxes)
        H = lambda state: cipher.mc_binary_matrix*state
        dim = cipher.S.input_size()
        n = nbr_sboxes * dim
    elif r == 2:
        nbr_sboxes = cipher.nbr_sboxes
        S = cipher.sbox_layer_faster
        H = cipher.linear_layer
        dim = cipher.S.input_size()
        n = nbr_sboxes * dim
    elif (r == 3) and isinstance(cipher, AESLikeCipher):
        nbr_sboxes = cipher.nbr_superboxes
        S = cipher.sbox_layer_faster
        H = lambda state: cipher.sc(cipher.mc(cipher.sc(state)))
        dim = (cipher.nbr_sboxes // cipher.nbr_superboxes) * cipher.S.input_size()
        n = nbr_sboxes * dim
    elif r == 3:
        nbr_sboxes = cipher.nbr_sboxes
        H = lambda state: cipher.linear_layer(cipher.sbox_layer_faster(cipher.linear_layer(state)))
        S = cipher.sbox_layer # with symbolic evaluation
        dim = cipher.S.input_size()
        n = nbr_sboxes * dim
    elif (r == 4) and isinstance(cipher, AESLikeCipher):
        nbr_sboxes = cipher.nbr_superboxes
        S = cipher.sbox_layer
        H = lambda state: cipher.mc(cipher.sbox_layer_faster(cipher.sc(cipher.mc(cipher.sc(state)))))
        dim = (cipher.nbr_sboxes // cipher.nbr_superboxes) * cipher.S.input_size()
        n = nbr_sboxes * dim
    else:
        return "NotChecked"

    if (r == 2) or ((r==3) and isinstance(cipher, AESLikeCipher)):
        PR = BooleanPolynomialRing(n, names=[f"b{i}" for i in range(n)])
        # no key needed so we set it to the all zero vector
        # (so that we can use the same code for all cases)
        K = vector(GF(2), n)
    else:
        PR = BooleanPolynomialRing(2*n, names=[f"{c}{i}" for c in ["b", "k"] for i in range(n)])
        K = vector(PR, n, PR.gens()[n:2*n]) # key added before last sbox layer
    B = vector(PR, n, PR.gens()[:n]) # output maske
    d = [vector(GF(2), n,  i*dim*[0]+dim*[1]+(n-(i+1)*dim)*[0]) for i in range(nbr_sboxes)]

    # generate system of (non-)linear equations
    def random_equation():
        x = random_vector(GF(2), n)
        z = S(H(x) + K)
        if is_even(nbr_sboxes):
            z += S(H(vector(GF(2), n, [0]*n))+K)
        for i in range(nbr_sboxes):
            x_ = d[i].pairwise_product(x)
            z += S(H(x_)+K)
        return z
    EQs = [random_equation() for _ in range(3*n)]

    # solve system
    if r == 2 or ((r == 3) and isinstance(cipher, AESLikeCipher)):
        # system is linear
        M = matrix(GF(2), len(EQs), n, EQs)
        return M.right_kernel()

    # else: system is non-linear
    SYS = [B.inner_product(eq) for eq in EQs]


    I = Ideal(SYS)
    try:
        alarm(300)
        GB = I.groebner_basis()
    except AlarmInterrupt:
        print("Groebner basis computation stopped (timeout)")
        to_magma(SYS, str(cipher), r)
        return "NotChecked"
    except Exception as e:
        cancel_alarm()
        print("Error during Groebner basis computation")
        #print(e)
        to_magma(SYS, str(cipher), r)
        return "NotChecked"
    else:
        cancel_alarm()

    # unfortunately .variety() is buggy in sage
    # when the solution is just linear
    # (also see https://groups.google.com/g/sage-support/c/hLfSZNc6w7A)
    # as a workaround we check the groebner basis manually first
    if all([b in GB for b in B]):
        # all bits of beta must be zero
        return [vector(GF(2), n)]
    try:
        alarm(60)
        V = Ideal(GB).variety()
    except AlarmInterrupt:
        print("Variety computation stopped (timeout)")
        print(list(GB))
        to_magma(SYS, str(cipher), r)
        return "NotChecked"
    except Exception as e:
        cancel_alarm()
        print(list(GB))
        to_magma(SYS, str(cipher), r)
        print("Error during variety computation")
        #print(e)
        return "NotChecked"
    else:
        cancel_alarm()

    return [vector(GF(2), n, [sol[b] for b in B]) for sol in V]



# generate textfile that can be loaded into magma
# to do the groebner basis computation
def to_magma(SYS, name, r):
    fn = f"{name}-{r}-rounds-equations.magma"
    print(f"Generating magma file {fn}")
    PR = SYS[0].ring()
    with open(fn, "w") as dat:
        dat.write("clear;\n")
        s = "P<"
        for i, a in enumerate(PR.gens()):
            if i != len(PR.gens()) - 1:
                s += str(a) + ","
            else:
                s += str(a)
        s += f">:= PolynomialRing(GF(2),{len(PR.gens())});"
        dat.write(s)
        dat.write("\nS:=[")
        for a in PR.gens():
            dat.write(f"{a}+{a}^2,")
        for i, eq in enumerate(SYS):
            s = str(eq)+"," if i != len(SYS)-1 else str(eq)
            dat.write(s)
        dat.write("];\n")
        dat.write('SetVerbose("Groebner", 1);\n')
        dat.write('GB:=GroebnerBasis(S);\n')
        dat.write('V:=Variety(Ideal(GB));\n')
        dat.write('#V;\n')


def check_approximation(A, B, F, N=10):
    C = 0
    for _ in range(N):
        x = random_vector(GF(2), len(A))
        y = F(x)
        c = A.inner_product(x) + B.inner_product(y)
        C += int(c)
    print(f"Pr[Ax = BF(x) + 1] = {C/N}")




