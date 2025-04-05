import secrets

from gf2bv import BitVec, LinearSystem


def magic(x, y):
    MASK64 = (1 << 64) - 1
    z1 = ((x ^ (y >> 22) ^ (x << 13)) & MASK64) >> 3
    z2 = ((y ^ (x >> 7) ^ (y << 5)) & MASK64) >> 3
    z3 = (x ^ y) & 0b101101
    return z1, z2, z3


def solve(lin: LinearSystem, zeros: list[BitVec], expected: tuple[int, int, int]):
    # solving all solutions
    for sol in lin.solve_all(zeros):
        print(f"{sol = }")
        assert magic(*sol) == expected

    # solving only one solution
    sol = lin.solve_one(zeros)
    print(f"{sol = }")
    assert magic(*sol) == expected

    # we also check evaluating the bitvectors at the zeros gives 0
    for z in zeros:
        assert lin.evaluate(z, sol) == 0


def simple_linear():
    lin = LinearSystem((64, 64))
    xs, ys = lin.gens()
    z1s, z2s, z3s = magic(xs, ys)
    zeros = [z1s, z2s, z3s]
    assert all([e & 1 == 0 for e in lin.get_eqs(zeros)]), "the system is not linear"

    solve(lin, zeros, (0, 0, 0))


def simple_affine():
    inp = secrets.randbits(64), secrets.randbits(64)
    print(f"{inp = }")
    z1, z2, z3 = magic(*inp)

    lin = LinearSystem((64, 64))
    xs, ys = lin.gens()
    z1s, z2s, z3s = magic(xs, ys)
    zeros = [z1s ^ z1, z2s ^ z2, z3s ^ z3]

    solve(lin, zeros, (z1, z2, z3))


if __name__ == "__main__":
    simple_linear()
    simple_affine()
