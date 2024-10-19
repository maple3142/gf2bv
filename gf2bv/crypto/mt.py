import random
from .. import BitVec


class MersenneTwister:
    def __init__(self, mt, w, n, m, r, a, u, d, s, b, t, c, l):
        w1 = (1 << w) - 1
        if len(mt) != n or min(r, u, s, t, l) > w and max(a, b, c, d) > w1:
            raise ValueError("invalid parameters")

        self.mt = list(mt)
        self.w = w
        self.n = n
        self.m = m
        self.r = r
        self.a = a
        self.u = u
        self.d = d
        self.s = s
        self.b = b
        self.t = t
        self.c = c
        self.l = l

        self.w1 = w1
        self.lmsk = w1 & ((1 << r) - 1)
        self.umsk = w1 ^ self.lmsk
        self.mti = 624

    def twist(self):
        for i in range(self.n):
            y = (self.mt[i] & self.umsk) ^ (self.mt[(i + 1) % self.n] & self.lmsk)
            sel = (
                y.broadcast(0, 32) & self.a
                if isinstance(y, BitVec)
                else (y & 1) * self.a
            )
            self.mt[i] = self.mt[(i + self.m) % self.n] ^ (y >> 1) ^ sel

    def temper(self, y):
        y ^= (y >> self.u) & self.d
        y ^= (y << self.s) & self.w1 & self.b
        y ^= (y << self.t) & self.w1 & self.c
        y ^= y >> self.l
        return y

    def __call__(self):
        if self.mti >= self.n:
            self.twist()
            self.mti = 0
        y = self.mt[self.mti]
        self.mti += 1
        return self.temper(y)

    def getrandbits(self, k=None):
        """Uses the CPython's implementation of random.getrandbits()"""
        if k is None:
            k = self.w
        if k == 0:
            return 0
        if k <= self.w:
            return self.__call__() >> (self.w - k)
        x = 0
        for i in range(0, k, self.w):
            r = self.__call__()
            if i + self.w > k:
                r = r >> (self.w - (k - i))
            x |= r << i
        return x


class MT19937(MersenneTwister):
    """32-bit Mersenne Twister by Matsumoto and Nishimura, 1998"""

    def __init__(self, mt):
        super().__init__(
            mt,
            32,
            624,
            397,
            31,
            0x9908B0DF,
            11,
            0xFFFFFFFF,
            7,
            0x9D2C5680,
            15,
            0xEFC60000,
            18,
        )

    def to_python_random(self):
        r = random.Random(0)
        r.setstate((3, (*self.mt, self.mti), None))
        return r
