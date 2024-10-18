import secrets
from gf2bv import LinearSystem
from gf2bv.crypto.xoshiro import Xoshiro256starstar


def xoshiro256starstar():
    xos = Xoshiro256starstar.generate()
    print(f"{xos.s = }")
    out = [xos() for _ in range(10)]

    lin = LinearSystem([64] * 4)
    xos2 = Xoshiro256starstar(lin.gens())
    zeros = [xos2.step() ^ Xoshiro256starstar.untemper(o) for o in out]
    for sol in lin.solve_all(zeros):
        print(f"{sol = }")
        xos2 = Xoshiro256starstar(sol)
        assert all(xos2() == o for o in out)


if __name__ == "__main__":
    xoshiro256starstar()
