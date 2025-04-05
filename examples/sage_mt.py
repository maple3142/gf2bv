import random
from gf2bv import LinearSystem
from gf2bv.crypto.mt import MT19937
from sage.all import vector, GF
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


def sage_test():
    bs = 32
    rand = random.Random(1234)

    effective_bs = ((bs - 1) & bs) or bs
    out = [rand.getrandbits(bs) for _ in range(624 * 32 // effective_bs)]

    lin = LinearSystem([32] * 624)
    mt = lin.gens()

    rng = MT19937(mt)
    zeros = [rng.getrandbits(bs) ^ o for o in out] + [mt[0] ^ 0x80000000]
    with timeit("get_sage_mat"):
        A, b = lin.get_sage_mat(zeros)
    A.set_immutable()
    print("dim", A.dimensions())
    with timeit("solve_right"):
        s = A.solve_right(b)  # vector
    with timeit("solve_raw (one solution)"):
        ss = lin.solve_raw_one(zeros)  # our raw solution as a single int

    # and they are equal:
    assert vector(GF(2), f"{ss:019968b}"[::-1]) == s


if __name__ == "__main__":
    sage_test()
