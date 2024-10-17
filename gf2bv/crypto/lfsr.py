import secrets
from .. import BitVec


class GaloisLFSR:
    def __init__(self, n: int, mask: int, state: int):
        M = (1 << n) - 1
        self.mask = mask & M
        self.state = state & M

    def __call__(self):
        bit = self.state & 1
        self.state >>= 1
        sel = (
            bit[0].dup(len(self.state)) & self.mask
            if isinstance(bit, BitVec)
            else bit * self.mask
        )
        self.state ^= sel
        return bit


class FibonacciLFSR:
    def __init__(self, n: int, mask: int, state: int):
        self.n = n
        M = (1 << n) - 1
        self.mask = mask & M
        self.state = state & M

    def __call__(self):
        b = self.state & 1
        if isinstance(self.state, BitVec):
            o = (self.state & self.mask).sum()
            self.state = (self.state >> 1) ^ o.zeroext(self.n - 1) << (self.n - 1)
        else:
            self.state = (self.state >> 1) | (
                ((self.state & self.mask).bit_count() & 1) << (self.n - 1)
            )
        return b
