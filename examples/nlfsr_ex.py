from gf2bv import QuadraticSystem, DimensionTooLargeError
from gf2bv.crypto.lfsr import GaloisLFSR, FibonacciLFSR
from tqdm import trange
from pathlib import Path
import secrets, pickle, traceback, itertools, gzip

from nlfsr import n, mask, select, non_linear_output


def nlfsr_ex_test(LFSR):
    print(f"Testing {LFSR.__name__}")

    # not 2**14 + 1000 because I want to make it possible to result in DimensionTooLargeError
    N = 2**14

    # just to show that quadratic system also works when multiple sizes are given
    qsys = QuadraticSystem([64, 64])
    lo, hi = qsys.gens()
    x = lo.concat(hi)

    # note that for a given LFSR, its system is the same
    # so we can precompute the system and reuse it
    cache_file_name = Path(__file__).parent / f"cache_{LFSR.__name__}.pkl.gz"
    try:
        with gzip.open(cache_file_name, "rb") as f:
            maybe_zeros = pickle.load(f)
        assert len(maybe_zeros) == N
        print("cache found, reusing...")
    except:
        print("cache not found, generating...")
        lfsr_sys = LFSR(128, mask, x)
        maybe_zeros = []
        for o in trange(N):
            lfsr_sys()
            x0, x1, x2, x3, x4 = [lfsr_sys.state.bits[i] for i in select]
            # this is same as the annihilator function
            z = qsys.mul_bit(x0, x1) ^ x0 ^ qsys.mul_bit(x1, x2) ^ x1 ^ x2 ^ 1
            # maybe_zeros is like zeros, it is also the real zeros if the output is 1
            maybe_zeros.append(z)
        with gzip.open(cache_file_name, "wb") as f:
            # yes, it can be trivially serialized using pickle
            pickle.dump(maybe_zeros, f)
            # dumping the system is also supported, but you don't need to do that

    # then generate the outputs
    init = secrets.randbits(n)
    print(f"{init = :0{n}b}")
    lfsr = LFSR(n, mask, init)
    out = [non_linear_output(lfsr) for _ in range(N)]

    # finally, use the outputs to generate the real zeros
    zeros = [z for z, o in zip(maybe_zeros, out) if o == 1]
    print(f"{len(zeros) = }")

    # try to solve it
    try:
        (sol,) = qsys.solve_one(zeros)
        print(f"{sol = :0{n}b}")
        assert sol == init
        print(
            "Lucky, the number of zeros is sufficient to not result in DimensionTooLargeError"
        )
    except DimensionTooLargeError as err:
        traceback.print_exc()
        print("=" * 40)
        # Unluckily, the solution space is too large
        # So we have to resort to bruteforcing some known bits
        # note that we are dealing we a quadratic system by linearization
        # so asserting a bit equality is more than just adding `x.bits[i] ^ bit` to the zeros
        # QuadraticSystem has a helper method to do this called `bit_assert`
        # for example, we can try to bruteforce 2 bits in x like this:
        for b0, b1 in itertools.product([0, 1], repeat=2):
            sol_tuple = qsys.solve_one(
                zeros
                # we can assert a single bit:
                + qsys.bit_assert(x.bits[0], b0)
                # and some linear combinations of bits:
                + qsys.bit_assert(x.bits[1] ^ x.bits[2] ^ x.bits[87], b1)
            )
            print(b0, b1, sol_tuple)
            if sol_tuple:
                lo, hi = sol_tuple
                s = lo | (hi << 64)
                assert s == init
                # check if the solution does match our guess
                assert s & 1 == b0
                assert ((s >> 1) & 1) ^ ((s >> 2) & 1) ^ ((s >> 87) & 1) == b1


if __name__ == "__main__":
    nlfsr_ex_test(GaloisLFSR)
    nlfsr_ex_test(FibonacciLFSR)
