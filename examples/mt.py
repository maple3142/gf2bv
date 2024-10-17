import random
from gf2bv import LinearSystem
from gf2bv.crypto.mt import MT19937


def mt19937(bs):
    rand = random.Random()
    st = tuple(rand.getstate()[1][:-1])

    effective_bs = ((bs - 1) & bs) or bs
    out = [rand.getrandbits(bs) for _ in range(624 * 32 // effective_bs)]

    lin = LinearSystem([32] * 624)
    mt = lin.gens()

    rng = MT19937(mt)
    zeros = [rng.getrandbits(bs) ^ o for o in out] + [mt[0] ^ 0x80000000]
    print("solving...")
    for sol in lin.solve_all(zeros):
        print("solved", sol[:10])
        assert sol == st

        rng = MT19937(sol)
        pyrand = rng.to_python_random()
        assert all(rng.getrandbits(bs) == o for o in out)
        assert all(pyrand.getrandbits(bs) == o for o in out)


if __name__ == "__main__":
    mt19937(32)
    mt19937(17)
    mt19937(9)
    mt19937(1)
