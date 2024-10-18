from gf2bv import QuadraticSystem
from gf2bv.crypto.lfsr import GaloisLFSR, FibonacciLFSR
from tqdm import tqdm
import secrets, itertools


n, mask = 128, 0xD670201BAC7515352A273372B2A95B23
select = (13, 24, 35, 46, 57)


def combiner(x0, x1, x2, x3, x4):
    # this combining function is even: 50% of the time it will return 1
    return (x0 * x1) ^ (x0 * x1 * x3 * x4) ^ x0 ^ x1 ^ x2


def non_linear_output(lfsr):
    lfsr()
    x0, x1, x2, x3, x4 = [(lfsr.state >> i) & 1 for i in select]
    return combiner(x0, x1, x2, x3, x4)


def annihilator(x0, x1, x2, x3, x4):
    # you can find this using sage.crypto.boolean_function.BooleanFunction
    return (x0 * x1) ^ x0 ^ (x1 * x2) ^ x1 ^ x2 ^ 1


def sanity_check():
    for x0, x1, x2, x3, x4 in itertools.product([0, 1], repeat=5):
        if combiner(x0, x1, x2, x3, x4) == 1:
            assert annihilator(x0, x1, x2, x3, x4) == 0


def nlfsr_test(LFSR):
    # this example shows how QuadraticSystem can be used to solve a non-linear LFSR by linearization
    print(f"Testing {LFSR.__name__}")
    init = secrets.randbits(n)
    print(f"{init = :0{n}b}")
    lfsr = LFSR(n, mask, init)

    N = 2**14 + 1000
    out = [non_linear_output(lfsr) for _ in range(N)]

    qsys = QuadraticSystem([128])
    (x,) = qsys.gens()
    lfsr_sys = LFSR(128, mask, x)
    zeros = []
    for o in tqdm(out):
        lfsr_sys()
        if o == 1:
            x0, x1, x2, x3, x4 = [lfsr_sys.state.bits[i] for i in select]
            # this is same as the annihilator function
            z = qsys.mul_bit(x0, x1) ^ x0 ^ qsys.mul_bit(x1, x2) ^ x1 ^ x2 ^ 1
            zeros.append(z)
    print(f"{len(zeros) = }")
    sols = list(qsys.solve_all(zeros))
    for (sol,) in sols:
        print(f"{sol = :0{n}b}")
        assert sol == init

    (sol,) = qsys.solve_one(zeros)
    print(f"{sol = :0{n}b}")
    assert sol == init


if __name__ == "__main__":
    sanity_check()
    nlfsr_test(GaloisLFSR)
    nlfsr_test(FibonacciLFSR)
