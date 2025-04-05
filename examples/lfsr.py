import secrets

from gf2bv import LinearSystem
from gf2bv.crypto.lfsr import FibonacciLFSR, GaloisLFSR


def lfsr_test(LFSR, n: int, mask: int):
    print(f"Testing {LFSR.__name__}")
    init_st = secrets.randbits(n)
    print(f"{init_st = :#x}")
    lfsr = LFSR(n, mask, init_st)
    out = [lfsr() for _ in range(256)]

    lin = LinearSystem([n])
    (sym,) = lin.gens()
    lfsr2 = LFSR(n, mask, sym)
    zeros = [lfsr2() ^ o for o in out]
    for (sol,) in lin.solve_all(zeros):
        print(f"{sol = :#x}")
        assert sol == init_st
    print()


if __name__ == "__main__":
    lfsr_test(GaloisLFSR, 128, 0x5C2B76970103D4EEFCD4A2C681CC400D)
    lfsr_test(FibonacciLFSR, 128, 0x6D6AC812F52A212D5A0B9F3117801FD5)
