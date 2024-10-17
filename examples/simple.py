import secrets
from gf2bv import LinearSystem


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

def mt19937_32():
    rand = random.Random(1234)
    st = tuple(rand.getstate()[1][:-1])
    out = [rand.getrandbits(32) for _ in range(624)]

    lin = LinearSystem(**{f"mt_{i}": 32 for i in range(624)})
    mt = lin.gens()

    rng = MT19937(mt)
    zeros = [rng.getrandbits(32) ^ o for o in out] + [mt[0] ^ 0x80000000]
    print("solving...")
    for sol in lin.solve(zeros):
        print("solved", sol[:10])
        assert sol == st
        rng = MT19937(sol)
        assert all(rng.getrandbits(32) == o for o in out)

        rand2 = random.Random(1234)
        rand2.setstate((3, tuple(list(sol) + [624]), None))
        assert all(rand2.getrandbits(32) == o for o in out)


def mt19937_1():
    rand = random.Random(1234)
    st = tuple(rand.getstate()[1][:-1])
    out = [rand.getrandbits(1) for _ in range(624 * 32)]

    lin = LinearSystem(**{f"mt_{i}": 32 for i in range(624)})
    mt = lin.gens()

    rng = MT19937(mt)
    zeros = [rng.getrandbits(1) ^ o for o in out] + [mt[0] ^ 0x80000000]
    print("solving...")
    for sol in lin.solve(zeros):
        # assert sol == st
        rng = MT19937(sol)
        assert all(rng.getrandbits(1) == o for o in out)

        rand2 = random.Random(1234)
        rand2.setstate((3, tuple(list(sol) + [624]), None))
        assert all(rand2.getrandbits(1) == o for o in out)

if __name__ == "__main__":
    simple_linear()
    simple_affine()
# mt19937_32()
# mt19937_1()
