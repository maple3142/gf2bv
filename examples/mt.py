import random, time
from gf2bv import LinearSystem
from gf2bv.crypto.mt import MT19937
from time import perf_counter
from contextlib import contextmanager


@contextmanager
def timeit(task_name):
    start = perf_counter()
    try:
        yield
    finally:
        end = perf_counter()
        print(f"{task_name} took {end - start:.2f} seconds")


def mt19937(bs):
    print("bs:", bs)
    rand = random.Random(3142)
    st = tuple(rand.getstate()[1][:-1])

    effective_bs = ((bs - 1) & bs) or bs
    out = [rand.getrandbits(bs) for _ in range(624 * 32 // effective_bs)]

    lin = LinearSystem([32] * 624)
    mt = lin.gens()

    rng = MT19937(mt)
    with timeit("generate system"):
        zeros = [rng.getrandbits(bs) ^ o for o in out] + [mt[0] ^ 0x80000000]
    print("solving...")
    with timeit("solve_one"):
        sol = lin.solve_one(zeros)
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
