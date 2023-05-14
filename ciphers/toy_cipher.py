from sage.crypto.sbox import SBox
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF

from ciphers.cipher import Cipher


class DifferentialProbOneExample(Cipher):
    def __init__(self, S, L, nbr_sboxes, name, alpha, beta):
        self.alpha = alpha
        self.beta = beta
        super().__init__(S, L, nbr_sboxes, name)

    def test_differential(self):
        SLS = lambda v: self.sbox_layer(self.linear_layer(self.sbox_layer(v)))
        for x in GF(2) ** (self.S.input_size() * self.nbr_sboxes):
            if self.beta != SLS(x) + SLS(x + self.alpha):
                return False
        return True


class Example1(DifferentialProbOneExample):
    def __init__(self):
        S = SBox([28, 27, 30, 9, 23, 21, 29, 4, 25, 0, 8, 10, 13, 19, 15, 17, 12, 18, 14, 16, 31, 22, 5, 7, 26, 24, 11, 2, 20, 3, 6, 1])
        L = Matrix(GF(2), [
            [1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0],
            [0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],
            [0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],
            [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1],
            [0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
            [0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1],
            [0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0],
            [1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1],
            [0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0]
        ])
        alpha = vector(GF(2), list(S.to_bits(0x1d)) * 3)
        super().__init__(S, L, 3, "Example 1", alpha, alpha)


class Example2(DifferentialProbOneExample):
    def __init__(self):
        S = SBox([0x8, 0x0, 0xa, 0x2, 0x3, 0x5, 0x4, 0x7, 0x6, 0x9, 0x1, 0xb, 0xc, 0xd, 0xe, 0xf])
        L = Matrix(GF(2), [
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0],
        ])
        alpha = vector(GF(2), [0, 0, 0, 0] + S.to_bits(0xc))
        beta = vector(GF(2), [0, 0, 0, 0] + S.to_bits(0x5))
        super().__init__(S, L, 2, "Example 2", alpha, beta)
