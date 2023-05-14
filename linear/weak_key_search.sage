from ciphers.boomslang import Boomslang
from ciphers.craft import Craft
from ciphers.midori import Midori
from ciphers.skinny import Skinny

# search for approximations and corresponding weak keys
# as a starting point we take the results (i.e. groebner bases) of algorithm 1
# for the cipher and its inverse

def weak_key_search_two_rounds(cipher, GB, PR):
    dim = len(cipher.S)
    nbr_sboxes = 4
    n = dim*nbr_sboxes
    B = vector(PR, n, PR.gens()[n:2*n])

    single_S = lambda x: cipher.sbox_layer(x, nbr_sboxes=1)
    single_S_inverse = lambda x: cipher.sbox_layer_inverse(x, nbr_sboxes=1)
    S = lambda x: cipher.sbox_layer(x, nbr_sboxes)
    L = lambda x: cipher.mc_binary_matrix * x

    SYS = GB
    for i in range(nbr_sboxes):
        ai = vector(PR, dim, PR.gens()[:n][dim*i:dim*(i+1)])
        ki = vector(PR, dim, PR.gens()[2*n:3*n][dim*i:dim*(i+1)])
        for _ in range(30):
            x = vector(GF(2), n)
            x_ = vector(GF(2), n)
            while x == x_:
                xi, xi_ = random_vector(GF(2), dim), random_vector(GF(2), dim)
                x[dim*i:dim*(i+1)] = xi
                x_[dim*i:dim*(i+1)] = xi_
            y = single_S_inverse(xi + ki) + single_S_inverse(xi_ + ki)
            z = S(L(x)) + S(L(x_))
            SYS.append(ai.inner_product(y) + B.inner_product(z))
    I = Ideal(SYS)
    GB = I.groebner_basis()
    print(f"Groebner Basis for masks and weak keys of {cipher} Superbox:")
    print(list(GB))

def weak_key_search_three_rounds(cipher, GB, PR):
    dim = 4*len(cipher.S) # len(superbox)
    nbr_sboxes = cipher.nbr_superboxes
    n = cipher.nbr_sboxes * len(cipher.S)

    A = vector(PR, n, PR.gens()[:n])
    A = cipher.sc_inverse(A)
    B = vector(PR, n, PR.gens()[n:2*n])
    K0 = vector(PR, n, PR.gens()[2*n:3*n]) # directly after first sbox layer
    K1 = vector(PR, n, PR.gens()[3*n:4*n]) # directly after second sbox layer

    mc_inv = cipher.mc_binary_matrix.inverse()
    def single_S_inverse(x, k):
        h = cipher.sbox_layer_inverse(x, nbr_sboxes=4)
        h = mc_inv * h
        h = h + k
        h = cipher.sbox_layer_inverse(h, nbr_sboxes=4)
        return h

    S = lambda x: cipher.sbox_layer(x, nbr_sboxes)
    L = lambda x: cipher.sc(cipher.mc(cipher.sc(x)))
    SYS = GB
    for i in range(nbr_sboxes):
        print('.', end="", flush=True)
        ai = vector(PR, dim, PR.gens()[:n][dim*i:dim*(i+1)])
        ki = K0[dim*i:dim*(i+1)]
        ki_ = K1[dim*i:dim*(i+1)]
        for _ in range(5*n):
            x = vector(GF(2), n)
            x_ = vector(GF(2), n)
            while x == x_:
                xi, xi_ = random_vector(GF(2), dim), random_vector(GF(2), dim)
                x[dim*i:dim*(i+1)] = xi
                x_[dim*i:dim*(i+1)] = xi_
            y = single_S_inverse(xi + ki_, ki) + single_S_inverse(xi_ + ki_, ki)
            z = S(L(x)) + S(L(x_))
            SYS.append(ai.inner_product(y) + B.inner_product(z))
    I = Ideal(SYS)
    GB = I.groebner_basis(full_prot=True)
    print(f"Groebner Basis for masks and weak keys for 3-round {cipher}")
    print(list(GB))


def weak_key_search_four_rounds(cipher, GB, PR):
    dim = 4*len(cipher.S) # len(superbox)
    nbr_sboxes = cipher.nbr_superboxes
    n = cipher.nbr_sboxes * len(cipher.S)

    A = vector(PR, n, PR.gens()[:n])
    B = vector(PR, n, PR.gens()[n:2*n])
    K0 = vector(PR, n, PR.gens()[2*n:3*n]) # directly after first sbox layer
    K1 = vector(PR, n, PR.gens()[3*n:4*n]) # directly after second sbox layer
    K2 = vector(PR, n, PR.gens()[4*n:5*n]) # directly before last sbox layer

    mc_inv = cipher.mc_binary_matrix.inverse()
    def single_S_inverse(x, k):
        h = cipher.sbox_layer_inverse(x, nbr_sboxes=4)
        h = mc_inv * h
        h = h + k
        h = cipher.sbox_layer_inverse(h, nbr_sboxes=4)
        return h

    S = lambda x: cipher.sbox_layer(x, nbr_sboxes)
    H = lambda x: cipher.mc(cipher.sbox_layer(cipher.sc(cipher.mc(cipher.sc(x)))))
    SYS = GB
    for i in range(nbr_sboxes):
        print('.', end="", flush=True)
        ai = vector(PR, dim, PR.gens()[:n][dim*i:dim*(i+1)])
        ki = K0[dim*i:dim*(i+1)]
        ki_ = K1[dim*i:dim*(i+1)]
        for _ in range(5*n):
            x = vector(GF(2), n)
            x_ = vector(GF(2), n)
            while x == x_:
                xi, xi_ = random_vector(GF(2), dim), random_vector(GF(2), dim)
                x[dim*i:dim*(i+1)] = xi
                x_[dim*i:dim*(i+1)] = xi_
            y = single_S_inverse(xi + ki_, ki) + single_S_inverse(xi_ + ki_, ki)
            z = S(H(x)+K2) + S(H(x_)+K2)
            SYS.append(ai.inner_product(y) + B.inner_product(z))
    I = Ideal(SYS)
    GB = I.groebner_basis(full_prot=True)
    print(f"Groebner Basis for masks and weak keys for 4-round {cipher}")
    print(list(GB))



# Two rounds

# groebner bases GB are essentially the result of Algorithm 1 # (linear_check(cipher(), 2))
PR = BooleanPolynomialRing(n=48, names=sum([[f"{c}{i}" for i in range(16)] for c in ["a", "b", "k"]], []))
B = PR.gens()[16:32]
GB = list(B[::2])
weak_key_search_two_rounds(Boomslang(), GB, PR)

GB = list(B[:8])
weak_key_search_two_rounds(Craft(), GB, PR)

GB = [B[0], B[4], B[8], B[12],
      B[1] + B[3], B[1] + B[5], B[1]+B[7], B[1]+B[9], B[1]+B[11], B[1]+B[13], B[1]+B[15],
      B[2] + B[6], B[2] + B[10], B[2] + B[14]]
weak_key_search_two_rounds(Midori(), GB, PR)

GB = list(B[:4]) + list(B[8:])
weak_key_search_two_rounds(Skinny(64), GB, PR)

PR = BooleanPolynomialRing(n=96, names=sum([[f"{c}{i}" for i in range(32)] for c in ["a", "b", "k"]], []))
B = PR.gens()[32:64]
GB = list(B[:8]) + list(B[16:32])
weak_key_search_two_rounds(Skinny(128), GB, PR)



# Three rounds

PR = BooleanPolynomialRing(n=4*128, names=sum([[f"{c}{i}" for i in range(128)]
                                               for c in ["a", "b", "k", "k_"]], []))
A = PR.gens()[:128]
B = PR.gens()[128:256]
# we now that half of the bits must be zero and the others are the same in each superbox
GB = list(A[1::2]) + list(B[::2])
for i in range(8): # suberboxen
    for j in range(1, 8):
        GB.append(A[16*i] + A[16*i + 2*j])
        GB.append(B[16*i+1] + B[16*i+1 + 2*j])
weak_key_search_three_rounds(Boomslang(), GB, PR)

GB = list(B[:8]) + list(B[16:40]) + list(B[48:72]) + list(B[80:104]) + list(B[112:128])
GB += list(A[8:32]) + list(A[40:64]) + list(A[72:96]) + list(A[104:128])
#weak_key_search_three_rounds(Skinny(128), GB, PR) # takes too long


PR = BooleanPolynomialRing(n=4*64, names=sum([[f"{c}{i}" for i in range(64)] for c in ["a", "b", "k", "k_"]], []))
A = PR.gens()[:64]
B = PR.gens()[64:128]
GB = list(A[8:16]) + list(A[24:32]) + list(A[40:48]) + list(A[56:64]) + list(B[8:16]) + list(B[24:32]) + list(B[40:48]) + list(B[56:64])
weak_key_search_three_rounds(Craft(), GB, PR)

PR.inject_variables()
GB = [b0, b4, b8, b12, b16, b20, b24, b28, b32, b36, b40, b44, b48, b52, b56, b60, b1 + b3, b1 + b21, b1 + b23, b1 + b41, b1 + b43, b1 + b61, b1 + b63, b2 + b22, b2 + b42, b2 + b62, b5 + b7, b5 + b17, b5 + b19, b5 + b45, b5 + b47, b5 + b57, b5 + b59, b6 + b18, b6 + b46, b6 + b58, b9 + b11, b9 + b29, b9 + b31, b9 + b33, b9 + b35, b9 + b53, b9 + b55, b10 + b30, b10 + b34, b10 + b54, b13 + b15, b13 + b25, b13 + b27, b13 + b37, b13 + b39, b13 + b49, b13 + b51, b14 + b26, b14 + b38, b14 + b50, a0, a4, a8, a12, a16, a20, a24, a28, a32, a36, a40, a44, a48, a52, a56, a60, a1 + a3, a1 + a29, a1 + a31, a1 + a37, a1 + a39, a1 + a57, a1 + a59, a2 + a30, a2 + a38, a2 + a58, a5 + a7, a5 + a25, a5 + a27, a5 + a33, a5 + a35, a5 + a61, a5 + a63, a6 + a26, a6 + a34, a6 + a62, a9 + a11, a9 + a21, a9 + a23, a9 + a45, a9 + a47, a9 + a49, a9 + a51, a10 + a22, a10 + a46, a10 + a50, a13 + a15, a13 + a17, a13 + a19, a13 + a41, a13 + a43, a13 + a53, a13 + a55, a14 + a18, a14 + a42, a14 + a54]
weak_key_search_three_rounds(Midori(), GB, PR)

GB = list(B[:4]) + list(B[8:20]) + list(B[24:36]) + list(B[40:52]) + list(B[56:64])
GB += list(A[4:16]) + list(A[20:32]) + list(A[36:48]) + list(A[52:64])
weak_key_search_three_rounds(Skinny(64), GB, PR)


# four rounds
PR = BooleanPolynomialRing(n=5*64, names=sum([[f"{c}{i}" for i in range(64)]
                                              for c in ["a", "b", "k", "k_", "k__"]], []))
PR.inject_variables()
GB = [b0, b1, b2, b3, b4, b5, b6, b7, b8*k__8 + b10*k__8 + b10*k__11 + b11*k__11, b8*k__9, b8*k__10 + b10*k__8 + b10*k__11 + b11*k__11, b8*k__11, b9*k__8 + b10*k__8 + b10*k__11 + b11*k__8 + b11*k__11, b9*k__9 + b11*k__11, b9*k__10 + b10*k__8 + b10*k__11 + b11*k__10 + b11*k__11, b9*k__11 + b11*k__11, b10*b11*k__11 + b11*k__11, b10*k__8*k__10 + b10*k__8 + b10*k__10*k__11 + b10*k__11 + b11*k__10*k__11 + b11*k__11, b10*k__8*k__11 + b10*k__11 + b11*k__11, b10*k__9 + b10*k__11, b11*k__8*k__11, b11*k__9 + b11*k__11, b12*k__12 + b14*k__12 + b14*k__15 + b15*k__15, b12*k__13, b12*k__14 + b14*k__12 + b14*k__15 + b15*k__15, b12*k__15, b13*k__12 + b14*k__12 + b14*k__15 + b15*k__12 + b15*k__15, b13*k__13 + b15*k__15, b13*k__14 + b14*k__12 + b14*k__15 + b15*k__14 + b15*k__15, b13*k__15 + b15*k__15, b14*b15*k__15 + b15*k__15, b14*k__12*k__14 + b14*k__12 + b14*k__14*k__15 + b14*k__15 + b15*k__14*k__15 + b15*k__15, b14*k__12*k__15 + b14*k__15 + b15*k__15, b14*k__13 + b14*k__15, b15*k__12*k__15, b15*k__13 + b15*k__15, b16, b17, b18, b19, b20, b21, b22, b23, b24*k__24 + b26*k__24 + b26*k__27 + b27*k__27, b24*k__25, b24*k__26 + b26*k__24 + b26*k__27 + b27*k__27, b24*k__27, b25*k__24 + b26*k__24 + b26*k__27 + b27*k__24 + b27*k__27, b25*k__25 + b27*k__27, b25*k__26 + b26*k__24 + b26*k__27 + b27*k__26 + b27*k__27, b25*k__27 + b27*k__27, b26*b27*k__27 + b27*k__27, b26*k__24*k__26 + b26*k__24 + b26*k__26*k__27 + b26*k__27 + b27*k__26*k__27 + b27*k__27, b26*k__24*k__27 + b26*k__27 + b27*k__27, b26*k__25 + b26*k__27, b27*k__24*k__27, b27*k__25 + b27*k__27, b28*k__28 + b30*k__28 + b30*k__31 + b31*k__31, b28*k__29, b28*k__30 + b30*k__28 + b30*k__31 + b31*k__31, b28*k__31, b29*k__28 + b30*k__28 + b30*k__31 + b31*k__28 + b31*k__31, b29*k__29 + b31*k__31, b29*k__30 + b30*k__28 + b30*k__31 + b31*k__30 + b31*k__31, b29*k__31 + b31*k__31, b30*b31*k__31 + b31*k__31, b30*k__28*k__30 + b30*k__28 + b30*k__30*k__31 + b30*k__31 + b31*k__30*k__31 + b31*k__31, b30*k__28*k__31 + b30*k__31 + b31*k__31, b30*k__29 + b30*k__31, b31*k__28*k__31, b31*k__29 + b31*k__31, b32, b33, b34, b35, b36, b37, b38, b39, b40*k__40 + b42*k__40 + b42*k__43 + b43*k__43, b40*k__41, b40*k__42 + b42*k__40 + b42*k__43 + b43*k__43, b40*k__43, b41*k__40 + b42*k__40 + b42*k__43 + b43*k__40 + b43*k__43, b41*k__41 + b43*k__43, b41*k__42 + b42*k__40 + b42*k__43 + b43*k__42 + b43*k__43, b41*k__43 + b43*k__43, b42*b43*k__43 + b43*k__43, b42*k__40*k__42 + b42*k__40 + b42*k__42*k__43 + b42*k__43 + b43*k__42*k__43 + b43*k__43, b42*k__40*k__43 + b42*k__43 + b43*k__43, b42*k__41 + b42*k__43, b43*k__40*k__43, b43*k__41 + b43*k__43, b44*k__44 + b46*k__44 + b46*k__47 + b47*k__47, b44*k__45, b44*k__46 + b46*k__44 + b46*k__47 + b47*k__47, b44*k__47, b45*k__44 + b46*k__44 + b46*k__47 + b47*k__44 + b47*k__47, b45*k__45 + b47*k__47, b45*k__46 + b46*k__44 + b46*k__47 + b47*k__46 + b47*k__47, b45*k__47 + b47*k__47, b46*b47*k__47 + b47*k__47, b46*k__44*k__46 + b46*k__44 + b46*k__46*k__47 + b46*k__47 + b47*k__46*k__47 + b47*k__47, b46*k__44*k__47 + b46*k__47 + b47*k__47, b46*k__45 + b46*k__47, b47*k__44*k__47, b47*k__45 + b47*k__47, b48, b49, b50, b51, b52, b53, b54, b55, b56*k__56 + b58*k__56 + b58*k__59 + b59*k__59, b56*k__57, b56*k__58 + b58*k__56 + b58*k__59 + b59*k__59, b56*k__59, b57*k__56 + b58*k__56 + b58*k__59 + b59*k__56 + b59*k__59, b57*k__57 + b59*k__59, b57*k__58 + b58*k__56 + b58*k__59 + b59*k__58 + b59*k__59, b57*k__59 + b59*k__59, b58*b59*k__59 + b59*k__59, b58*k__56*k__58 + b58*k__56 + b58*k__58*k__59 + b58*k__59 + b59*k__58*k__59 + b59*k__59, b58*k__56*k__59 + b58*k__59 + b59*k__59, b58*k__57 + b58*k__59, b59*k__56*k__59, b59*k__57 + b59*k__59, b60*k__60 + b62*k__60 + b62*k__63 + b63*k__63, b60*k__61, b60*k__62 + b62*k__60 + b62*k__63 + b63*k__63, b60*k__63, b61*k__60 + b62*k__60 + b62*k__63 + b63*k__60 + b63*k__63, b61*k__61 + b63*k__63, b61*k__62 + b62*k__60 + b62*k__63 + b63*k__62 + b63*k__63, b61*k__63 + b63*k__63, b62*b63*k__63 + b63*k__63, b62*k__60*k__62 + b62*k__60 + b62*k__62*k__63 + b62*k__63 + b63*k__62*k__63 + b63*k__63, b62*k__60*k__63 + b62*k__63 + b63*k__63, b62*k__61 + b62*k__63, b63*k__60*k__63, b63*k__61 + b63*k__63] + [a0, a1, a2, a3, a4, a5, a6, a7, a8*k8 + a10*k8 + a10*k11 + a11*k11, a8*k9, a8*k10 + a10*k8 + a10*k11 + a11*k11, a8*k11, a9*k8 + a10*k8 + a10*k11 + a11*k8 + a11*k11, a9*k9 + a11*k11, a9*k10 + a10*k8 + a10*k11 + a11*k10 + a11*k11, a9*k11 + a11*k11, a10*a11*k11 + a11*k11, a10*k8*k10 + a10*k8 + a10*k10*k11 + a10*k11 + a11*k10*k11 + a11*k11, a10*k8*k11 + a10*k11 + a11*k11, a10*k9 + a10*k11, a11*k8*k11, a11*k9 + a11*k11, a12*k12 + a14*k12 + a14*k15 + a15*k15, a12*k13, a12*k14 + a14*k12 + a14*k15 + a15*k15, a12*k15, a13*k12 + a14*k12 + a14*k15 + a15*k12 + a15*k15, a13*k13 + a15*k15, a13*k14 + a14*k12 + a14*k15 + a15*k14 + a15*k15, a13*k15 + a15*k15, a14*a15*k15 + a15*k15, a14*k12*k14 + a14*k12 + a14*k14*k15 + a14*k15 + a15*k14*k15 + a15*k15, a14*k12*k15 + a14*k15 + a15*k15, a14*k13 + a14*k15, a15*k12*k15, a15*k13 + a15*k15, a16, a17, a18, a19, a20, a21, a22, a23, a24*k24 + a26*k24 + a26*k27 + a27*k27, a24*k25, a24*k26 + a26*k24 + a26*k27 + a27*k27, a24*k27, a25*k24 + a26*k24 + a26*k27 + a27*k24 + a27*k27, a25*k25 + a27*k27, a25*k26 + a26*k24 + a26*k27 + a27*k26 + a27*k27, a25*k27 + a27*k27, a26*a27*k27 + a27*k27, a26*k24*k26 + a26*k24 + a26*k26*k27 + a26*k27 + a27*k26*k27 + a27*k27, a26*k24*k27 + a26*k27 + a27*k27, a26*k25 + a26*k27, a27*k24*k27, a27*k25 + a27*k27, a28*k28 + a30*k28 + a30*k31 + a31*k31, a28*k29, a28*k30 + a30*k28 + a30*k31 + a31*k31, a28*k31, a29*k28 + a30*k28 + a30*k31 + a31*k28 + a31*k31, a29*k29 + a31*k31, a29*k30 + a30*k28 + a30*k31 + a31*k30 + a31*k31, a29*k31 + a31*k31, a30*a31*k31 + a31*k31, a30*k28*k30 + a30*k28 + a30*k30*k31 + a30*k31 + a31*k30*k31 + a31*k31, a30*k28*k31 + a30*k31 + a31*k31, a30*k29 + a30*k31, a31*k28*k31, a31*k29 + a31*k31, a32, a33, a34, a35, a36, a37, a38, a39, a40*k40 + a42*k40 + a42*k43 + a43*k43, a40*k41, a40*k42 + a42*k40 + a42*k43 + a43*k43, a40*k43, a41*k40 + a42*k40 + a42*k43 + a43*k40 + a43*k43, a41*k41 + a43*k43, a41*k42 + a42*k40 + a42*k43 + a43*k42 + a43*k43, a41*k43 + a43*k43, a42*a43*k43 + a43*k43, a42*k40*k42 + a42*k40 + a42*k42*k43 + a42*k43 + a43*k42*k43 + a43*k43, a42*k40*k43 + a42*k43 + a43*k43, a42*k41 + a42*k43, a43*k40*k43, a43*k41 + a43*k43, a44*k44 + a46*k44 + a46*k47 + a47*k47, a44*k45, a44*k46 + a46*k44 + a46*k47 + a47*k47, a44*k47, a45*k44 + a46*k44 + a46*k47 + a47*k44 + a47*k47, a45*k45 + a47*k47, a45*k46 + a46*k44 + a46*k47 + a47*k46 + a47*k47, a45*k47 + a47*k47, a46*a47*k47 + a47*k47, a46*k44*k46 + a46*k44 + a46*k46*k47 + a46*k47 + a47*k46*k47 + a47*k47, a46*k44*k47 + a46*k47 + a47*k47, a46*k45 + a46*k47, a47*k44*k47, a47*k45 + a47*k47, a48, a49, a50, a51, a52, a53, a54, a55, a56*k56 + a58*k56 + a58*k59 + a59*k59, a56*k57, a56*k58 + a58*k56 + a58*k59 + a59*k59, a56*k59, a57*k56 + a58*k56 + a58*k59 + a59*k56 + a59*k59, a57*k57 + a59*k59, a57*k58 + a58*k56 + a58*k59 + a59*k58 + a59*k59, a57*k59 + a59*k59, a58*a59*k59 + a59*k59, a58*k56*k58 + a58*k56 + a58*k58*k59 + a58*k59 + a59*k58*k59 + a59*k59, a58*k56*k59 + a58*k59 + a59*k59, a58*k57 + a58*k59, a59*k56*k59, a59*k57 + a59*k59, a60*k60 + a62*k60 + a62*k63 + a63*k63, a60*k61, a60*k62 + a62*k60 + a62*k63 + a63*k63, a60*k63, a61*k60 + a62*k60 + a62*k63 + a63*k60 + a63*k63, a61*k61 + a63*k63, a61*k62 + a62*k60 + a62*k63 + a63*k62 + a63*k63, a61*k63 + a63*k63, a62*a63*k63 + a63*k63, a62*k60*k62 + a62*k60 + a62*k62*k63 + a62*k63 + a63*k62*k63 + a63*k63, a62*k60*k63 + a62*k63 + a63*k63, a62*k61 + a62*k63, a63*k60*k63, a63*k61 + a63*k63]
weak_key_search_four_rounds(Craft(), GB, PR)


GB = [b0, b1, b2, b3, b4*k__4, b4*k__5, b5, b6*k__4, b6*k__5, b6*k__7 + b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20*k__20, b20*k__21, b21, b22*k__20, b22*k__21, b22*k__23 + b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36*k__36, b36*k__37, b37, b38*k__36, b38*k__37, b38*k__39 + b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b52*k__52, b52*k__53, b53, b54*k__52, b54*k__53, b54*k__55 + b54, b55, b56, b57, b58, b59, b60, b61, b62, b63]+[a0 + a1, a1*k0, a1*k1, a1*k3 + a1, a2*k0, a2*k1, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16 + a17, a17*k16, a17*k17, a17*k19 + a17, a18*k16, a18*k17, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32 + a33, a33*k32, a33*k33, a33*k35 + a33, a34*k32, a34*k33, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48 + a49, a49*k48, a49*k49, a49*k51 + a49, a50*k48, a50*k49, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63]
weak_key_search_four_rounds(Skinny(64), GB, PR)


# for boomslang, midori64 and mantis the two round masks work since they are iterative
A = vector(GF(2), ZZ(0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa).digits(2, padto=128)[::-1])
B = vector(GF(2), ZZ(0x55555555555555555555555555555555).digits(2, padto=128)[::-1])
S = Boomslang().sbox_layer_faster
L = Boomslang().linear_layer
F = lambda x: S(L(S(L(S(L(S(x)))))))
print("Checking 4-round Boomslang approximation")
check_approximation(A, B, F)

A = vector(GF(2), ZZ(0x5555555555555555).digits(2, padto=64)[::-1])
B = vector(GF(2), ZZ(0x5555555555555555).digits(2, padto=64)[::-1])
S = Midori().sbox_layer_faster
L = Midori().linear_layer
F = lambda x: S(L(S(L(S(L(S(x)))))))
print("Checking 4-round Midori approximation")
check_approximation(A, B, F)

A = vector(GF(2), ZZ(0x5555555555555555).digits(2, padto=64)[::-1])
B = vector(GF(2), ZZ(0x5555555555555555).digits(2, padto=64)[::-1])
S = Mantis().sbox_layer_faster
L = Mantis().linear_layer
F = lambda x: S(L(S(L(S(L(S(x)))))))
print("Checking 4-round Mantis approximation")
check_approximation(A, B, F)



# RESULTS
# Two Rounds
# Groebner Basis for masks and weak keys of Boomslang Superbox:
# [a0 + b15, a1, a2 + b15, a3, a4 + b15, a5, a6 + b15, a7, a8 + b15, a9, a10 + b15, a11, a12 + b15, a13, a14 + b15, a15, b0, b1 + b15, b2, b3 + b15, b4, b5 + b15, b6, b7 + b15, b8, b9 + b15, b10, b11 + b15, b12, b13 + b15, b14, b15*k0 + b15*k2, b15*k1 + b15*k3, b15*k4 + b15*k6, b15*k5 + b15*k7, b15*k8 + b15*k10, b15*k9 + b15*k11, b15*k12 + b15*k14, b15*k13 + b15*k15]

# Groebner Basis for masks and weak keys of Craft Superbox:
# [a0, a1, a2, a3, a4, a5, a6, a7, a8 + b8, a9 + b9 + b10*k8 + b10*k11 + b11*k11, a10 + b10, a11 + b10*k8 + b10*k11 + b11*k11 + b11, a12 + b12, a13 + b13 + b14*k12 + b14*k15 + b15*k15, a14 + b14, a15 + b14*k12 + b14*k15 + b15*k15 + b15, b0, b1, b2, b3, b4, b5, b6, b7, b8*k8 + b10*k8 + b10*k11 + b11*k11, b8*k9, b8*k10 + b10*k8 + b10*k11 + b11*k11, b8*k11, b9*k8 + b10*k8 + b10*k11 + b11*k8 + b11*k11, b9*k9 + b11*k11, b9*k10 + b10*k8 + b10*k11 + b11*k10 + b11*k11, b9*k11 + b11*k11, b10*b11*k11 + b11*k11, b10*k8*k10 + b10*k8 + b10*k10*k11 + b10*k11 + b11*k10*k11 + b11*k11, b10*k8*k11 + b10*k11 + b11*k11, b10*k9 + b10*k11, b11*k8*k11, b11*k9 + b11*k11, b12*k12 + b14*k12 + b14*k15 + b15*k15, b12*k13, b12*k14 + b14*k12 + b14*k15 + b15*k15, b12*k15, b13*k12 + b14*k12 + b14*k15 + b15*k12 + b15*k15, b13*k13 + b15*k15, b13*k14 + b14*k12 + b14*k15 + b15*k14 + b15*k15, b13*k15 + b15*k15, b14*b15*k15 + b15*k15, b14*k12*k14 + b14*k12 + b14*k14*k15 + b14*k15 + b15*k14*k15 + b15*k15, b14*k12*k15 + b14*k15 + b15*k15, b14*k13 + b14*k15, b15*k12*k15, b15*k13 + b15*k15]

# Groebner Basis for masks and weak keys of Midori64 Superbox:
# [a0, a1 + b15, a2 + b14, a3 + b15, a4, a5 + b15, a6 + b14, a7 + b15, a8, a9 + b15, a10 + b14, a11 + b15, a12, a13 + b15, a14 + b14, a15 + b15, b0, b1 + b15, b2 + b14, b3 + b15, b4, b5 + b15, b6 + b14, b7 + b15, b8, b9 + b15, b10 + b14, b11 + b15, b12, b13 + b15, b14*b15*k3 + b15*k3, b14*b15*k7 + b15*k7, b14*b15*k11 + b15*k11, b14*b15*k15 + b15*k15, b14*k0 + b14*k3 + b15*k3, b14*k1 + b14*k3, b14*k4 + b14*k7 + b15*k7, b14*k5 + b14*k7, b14*k8 + b14*k11 + b15*k11, b14*k9 + b14*k11, b14*k12 + b14*k15 + b15*k15, b14*k13 + b14*k15, b15*k0*k3, b15*k0*k7, b15*k0*k11, b15*k0*k15, b15*k1 + b15*k3, b15*k3*k4, b15*k3*k8, b15*k3*k12, b15*k4*k7, b15*k4*k11, b15*k4*k15, b15*k5 + b15*k7, b15*k7*k8, b15*k7*k12, b15*k8*k11, b15*k8*k15, b15*k9 + b15*k11, b15*k11*k12, b15*k12*k15, b15*k13 + b15*k15]

# Groebner Basis for masks and weak keys of SKINNY-64 Superbox:
# [a0 + b6, a1 + b6, a2 + b4 + b6*k2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, b0, b1, b2, b3, b4*k0, b4*k1, b5, b6*k0, b6*k1, b6*k3 + b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]

# Groebner Basis for masks and weak keys of SKINNY-128 Superbox:
#[a0, a1, a2, a3, a4 + b9 + b11*k2, a5, a6 + b11, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20 , a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, b0, b1, b2, b3, b4, b5, b6, b7, b8 + b11, b9*k0, b9*k1, b 10, b11*k0, b11*k1, b11*k3 + b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31]


# Three Rounds
# Groebner Basis for masks and weak keys for 3-round Boomslang
# [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a80, a81, a82, a83, a84, a85, a86, a87, a88, a89, a90, a91, a92, a93, a94, a95, a96, a97, a98, a99, a100, a101, a102, a103, a104, a105, a106, a107, a108, a109, a110, a111, a112, a113, a114, a115, a116, a117, a118, a119, a120, a121, a122, a123, a124, a125, a126, a127, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33 + b47, b34, b35 + b47, b36, b37 + b47, b38, b39 + b47, b40, b41 + b47, b42, b43 + b47, b44, b45 + b47, b46, b48, b49 + b63, b50, b51 + b63, b52, b53 + b63, b54, b55 + b63, b56, b57 + b63, b58, b59 + b63, b60, b61 + b63, b62, b64, b65 + b79, b66, b67 + b79, b68, b69 + b79, b70, b71 + b79, b72, b73 + b79, b74, b75 + b79, b76, b77 + b79, b78, b80, b81 + b95, b82, b83 + b95, b84, b85 + b95, b86, b87 + b95, b88, b89 + b95, b90, b91 + b95, b92, b93 + b95, b94, b96, b97 + b111, b98, b99 + b111, b100, b101 + b111, b102, b103 + b111, b104, b105 + b111, b106, b107 + b111, b108, b109 + b111, b110, b112, b113 + b127, b114, b115 + b127, b116, b117 + b127, b118, b119 + b127, b120, b121 + b127, b122, b123 + b127, b124, b125 + b127, b126]
# -> alpha must be all zero

# Groebner Basis for masks and weak keys for 3-round Craft
# [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b24, b25, b26, b27, b28, b29, b30, b31, b40, b41, b42, b43, b44, b45, b46, b47, b56, b57, b58, b59, b60, b61, b62, b63]
# -> alpha must be all zero

# Groebner Basis for masks and weak keys for 3-round Midori64
# [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36, b37, b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b52, b53, b54, b55, b56, b57, b58, b59, b60, b61, b62, b63]
# -> alpha must be all zero

# Groebner Basis for masks and weak keys for 3-round SKINNY-64
#[a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b56, b57, b58, b59, b60, b61, b62, b63]
# -> alpha must be all zero


# Four Rounds
# Craft
# [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24*k__24 + b26*k__24 + b26*k__27 + b27*k__27, b24*k__25, b24*k__26 + b26*k__24 + b26*k__27 + b27*k__27, b24*k__27, b25*k__24 + b26*k__24 + b26*k__27 + b27*k__24 + b27*k__27, b25*k__25 + b27*k__27, b25*k__26 + b26*k__24 + b26*k__27 + b27*k__26 + b27*k__27, b25*k__27 + b27*k__27, b26*b27*k__27 + b27*k__27, b26*k__24*k__26 + b26*k__24 + b26*k__26*k__27 + b26*k__27 + b27*k__26*k__27 + b27*k__27, b26*k__24*k__27 + b26*k__27 + b27*k__27, b26*k__25 + b26*k__27, b27*k__24*k__27, b27*k__25 + b27*k__27, b28*k__28 + b30*k__28 + b30*k__31 + b31*k__31, b28*k__29, b28*k__30 + b30*k__28 + b30*k__31 + b31*k__31, b28*k__31, b29*k__28 + b30*k__28 + b30*k__31 + b31*k__28 + b31*k__31, b29*k__29 + b31*k__31, b29*k__30 + b30*k__28 + b30*k__31 + b31*k__30 + b31*k__31, b29*k__31 + b31*k__31, b30*b31*k__31 + b31*k__31, b30*k__28*k__30 + b30*k__28 + b30*k__30*k__31 + b30*k__31 + b31*k__30*k__31 + b31*k__31, b30*k__28*k__31 + b30*k__31 + b31*k__31, b30*k__29 + b30*k__31, b31*k__28*k__31, b31*k__29 + b31*k__31, b32, b33, b34, b35, b36, b37, b38, b39, b40*k__40 + b42*k__40 + b42*k__43 + b43*k__43, b40*k__41, b40*k__42 + b42*k__40 + b42*k__43 + b43*k__43, b40*k__43, b41*k__40 + b42*k__40 + b42*k__43 + b43*k__40 + b43*k__43, b41*k__41 + b43*k__43, b41*k__42 + b42*k__40 + b42*k__43 + b43*k__42 + b43*k__43, b41*k__43 + b43*k__43, b42*b43*k__43 + b43*k__43, b42*k__40*k__42 + b42*k__40 + b42*k__42*k__43 + b42*k__43 + b43*k__42*k__43 + b43*k__43, b42*k__40*k__43 + b42*k__43 + b43*k__43, b42*k__41 + b42*k__43, b43*k__40*k__43, b43*k__41 + b43*k__43, b44*k__44 + b46*k__44 + b46*k__47 + b47*k__47, b44*k__45, b44*k__46 + b46*k__44 + b46*k__47 + b47*k__47, b44*k__47, b45*k__44 + b46*k__44 + b46*k__47 + b47*k__44 + b47*k__47, b45*k__45 + b47*k__47, b45*k__46 + b46*k__44 + b46*k__47 + b47*k__46 + b47*k__47, b45*k__47 + b47*k__47, b46*b47*k__47 + b47*k__47, b46*k__44*k__46 + b46*k__44 + b46*k__46*k__47 + b46*k__47 + b47*k__46*k__47 + b47*k__47, b46*k__44*k__47 + b46*k__47 + b47*k__47, b46*k__45 + b46*k__47, b47*k__44*k__47, b47*k__45 + b47*k__47, b48, b49, b50, b51, b52, b53, b54, b55, b56*k__56 + b58*k__56 + b58*k__59 + b59*k__59, b56*k__57, b56*k__58 + b58*k__56 + b58*k__59 + b59*k__59, b56*k__59, b57*k__56 + b58*k__56 + b58*k__59 + b59*k__56 + b59*k__59, b57*k__57 + b59*k__59, b57*k__58 + b58*k__56 + b58*k__59 + b59*k__58 + b59*k__59, b57*k__59 + b59*k__59, b58*b59*k__59 + b59*k__59, b58*k__56*k__58 + b58*k__56 + b58*k__58*k__59 + b58*k__59 + b59*k__58*k__59 + b59*k__59, b58*k__56*k__59 + b58*k__59 + b59*k__59, b58*k__57 + b58*k__59, b59*k__56*k__59, b59*k__57 + b59*k__59, b60*k__60 + b62*k__60 + b62*k__63 + b63*k__63, b60*k__61, b60*k__62 + b62*k__60 + b62*k__63 + b63*k__63, b60*k__63, b61*k__60 + b62*k__60 + b62*k__63 + b63*k__60 + b63*k__63, b61*k__61 + b63*k__63, b61*k__62 + b62*k__60 + b62*k__63 + b63*k__62 + b63*k__63, b61*k__63 + b63*k__63, b62*b63*k__63 + b63*k__63, b62*k__60*k__62 + b62*k__60 + b62*k__62*k__63 + b62*k__63 + b63*k__62*k__63 + b63*k__63, b62*k__60*k__63 + b62*k__63 + b63*k__63, b62*k__61 + b62*k__63, b63*k__60*k__63, b63*k__61 + b63*k__63]


# Groebner Basis for masks and weak keys for 4-round SKINNY-64
# [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, a61, a62, a63, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20*k__20, b20*k__21, b21, b22*k__20, b22*k__21, b22*k__23 + b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36*k__36, b36*k__37, b37, b38*k__36, b38*k__37, b38*k__39 + b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b52*k__52, b52*k__53, b53, b54*k__52, b54*k__53, b54*k__55 + b54, b55, b56, b57, b58, b59, b60, b61, b62, b63]
