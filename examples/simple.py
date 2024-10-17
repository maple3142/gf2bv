import secrets
from gf2bv import LinearSystem
from gf2bv.pym4ri import solve


def magic(x, y):
    MASK64 = (1 << 64) - 1
    z1 = ((x ^ (y >> 22) ^ (x << 13)) & MASK64) >> 3
    z2 = ((y ^ (x >> 7) ^ (y << 5)) & MASK64) >> 3
    z3 = (x ^ y) & 0b101101
    return z1, z2, z3


def simple_linear():
    lin = LinearSystem((64, 64))
    xs, ys = lin.gens()
    z1s, z2s, z3s = magic(xs, ys)
    zeros = [z1s, z2s, z3s]
    assert all([e & 1 == 0 for e in lin.get_eqs(zeros)]), "the system is not linear"
    for sol in lin.solve(zeros):
        print(f"{sol = }")
        assert magic(*sol) == (0, 0, 0)


def simple_affine():
    inp = secrets.randbits(64), secrets.randbits(64)
    print(f"{inp = }")
    z1, z2, z3 = magic(*inp)

    lin = LinearSystem((64, 64))
    xs, ys = lin.gens()
    z1s, z2s, z3s = magic(xs, ys)
    zeros = [z1s ^ z1, z2s ^ z2, z3s ^ z3]
    for sol in lin.solve(zeros):
        print(f"{sol = }")
        assert magic(*sol) == (z1, z2, z3)

if __name__ == "__main__":
    simple_linear()
    simple_affine()
