from .. import BitVec
import secrets

MASK64 = (1 << 64) - 1


def rotl64(x, n):
    if isinstance(x, BitVec):
        return x.rotl(n)
    return ((x << n) | (x >> (64 - n))) & MASK64


class Xoshiro256starstar:
    def __init__(self, s: list[int]):
        if len(s) != 4:
            raise ValueError("invalid state")
        self.s = s

    @staticmethod
    def generate():
        return Xoshiro256starstar([secrets.randbits(64) for _ in range(4)])

    @staticmethod
    def tamper(s1):
        return rotl64(s1 * 5 & MASK64, 7) * 9 & MASK64

    inv9 = pow(9, -1, 1 << 64)
    inv5 = pow(5, -1, 1 << 64)

    @staticmethod
    def untamper(s1):
        return (
            rotl64(s1 * Xoshiro256starstar.inv9 & MASK64, 64 - 7)
            * Xoshiro256starstar.inv5
            & MASK64
        )

    def step(self):
        s0, s1, s2, s3 = self.s
        result = s1
        t = (s1 << 17) & MASK64
        s2 ^= s0
        s3 ^= s1
        s1 ^= s2
        s0 ^= s3
        s2 ^= t
        s3 = rotl64(s3, 45)
        self.s = [s0, s1, s2, s3]
        return result

    def __call__(self):
        return self.tamper(self.step())
