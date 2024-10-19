from __future__ import annotations
from typing import Any, Union, Optional, Literal, TypeVar, TYPE_CHECKING
from collections.abc import Sequence
from operator import xor
from functools import reduce
from ._internal import m4ri_solve, to_bits, mul_bit_quad, AffineSpace

TSolveMode = TypeVar("TSolveMode", Literal[0], Literal[1])


class BitVec:
    __slots__ = ("bits",)

    def __init__(self, bits: list[int]):
        # this represent the symbolic bits of the BitVec
        # in little-endian order (lsb first)
        self.bits = bits

    def __copy__(self):
        return BitVec(self.bits[:])

    def __len__(self):
        return len(self.bits)

    def __getitem__(self, key: int | slice):
        if isinstance(key, slice):
            return BitVec(self.bits[key])
        return self.bits[key]

    def __xor__(self, other: BitVec | int):
        if not isinstance(other, BitVec):
            bs = to_bits(len(self.bits), other)
            return BitVec(list(map(xor, self.bits, bs)))
        else:
            if len(self.bits) != len(other.bits):
                raise ValueError("Cannot mix bitvecs of different lengths")
        return BitVec(list(map(xor, self.bits, other.bits)))

    __rxor__ = __xor__
    __pow__ = __xor__  # alias to __xor__, for convenience in sage

    def __rshift__(self, n: int):
        return BitVec(self.bits[n:] + [0] * n)

    def __lshift__(self, n: int):
        return BitVec([0] * n + self.bits[:-n])

    def __and__(self, mask: int):
        bs = to_bits(len(self.bits), mask)
        return BitVec([0 if not b else a for a, b in zip(self.bits, bs)])

    __rand__ = __and__

    def __or__(self, mask: int):
        bs = to_bits(len(self.bits), mask)
        return BitVec([1 if b else a for a, b in zip(self.bits, bs)])

    __ror__ = __or__

    def rotr(self, n: int):
        return BitVec(self.bits[n:] + self.bits[:n])

    def rotl(self, n: int):
        return BitVec(self.bits[-n:] + self.bits[:-n])

    def sum(self):
        return BitVec([reduce(xor, self.bits)])

    def zeroext(self, n: int):
        return BitVec(self.bits + [0] * n)

    def signext(self, n: int):
        return BitVec(self.bits + [self.bits[-1]] * n)

    def broadcast(self, i: int, n: int):
        return BitVec([self.bits[i]] * n)

    def dup(self, n: int):
        return BitVec(self.bits * n)

    def concat(self, other: BitVec):
        return BitVec(self.bits + other.bits)


Zeros = Sequence[BitVec | int]


class DimensionTooLargeError(Exception):
    def __init__(self, message: str, space: AffineSpace):
        super().__init__(message)
        self.space = space


class LinearSystem:
    def __init__(self, sizes: list[int]):
        self._sizes = sizes[:]
        self._cols = sum(sizes)

        # lsb used to represent constant term (affine part), the first one in the basis
        self._basis = [1 << i for i in range(1 + self._cols)]

        _vars: list[BitVec] = []
        i = 1  # lsb used to represent constant term (affine part)
        for size in self._sizes:
            _vars.append(BitVec(self._basis[i : i + size]))
            i += size
        self._vars = tuple(_vars)

    def gens(self):
        return self._vars

    def __reduce__(self):
        return (self.__class__, (self._sizes,))

    def get_sage_mat(self, zeros: Zeros, tqdm=lambda x, desc: x):
        """
        Convert the system of equations to Sage, return a matrix A and a vector b such that Ax = b
        """
        from sage.all import GF, vector, MatrixSpace, matrix

        F2 = GF(2)
        eqs = self.get_eqs(zeros)
        affine = []
        for i in range(len(eqs)):
            affine.append(eqs[i] & 1)
            eqs[i] >>= 1
        affine = vector(F2, affine)
        cols = self._cols
        rows = len(eqs)
        mat = matrix(F2, rows, cols)

        # this is slow, idk how to optimize this
        i = 0
        for v in tqdm(eqs, desc="Converting equations"):
            bs = to_bits(cols, v)
            for j in range(cols):
                if bs[j]:
                    mat[i, j] = 1
            i += 1
        return mat, affine

    def get_eqs(self, zeros: Zeros) -> list[int]:
        """
        Convert zeros into a list of equations accepted by pym4ri.solve
        """
        # https://stackoverflow.com/questions/716477/join-list-of-lists-in-python/21034265#21034265
        # this seems to be the fastest way to flatten a list of lists
        eqs: list[int] = []
        list(
            map(
                eqs.extend,
                (bv.bits if isinstance(bv, BitVec) else [bv] for bv in zeros),
            )
        )
        return list(filter(None, eqs))  # remove literal zeros

    def solve_raw(self, zeros: Zeros, mode: TSolveMode):
        eqs = self.get_eqs(zeros)
        if 1 in eqs:
            # no solution
            return
        cols = self._cols
        if cols > len(eqs):
            # pym4ri.solve requires rows >= cols, pad with zeros
            eqs += [0] * (cols - len(eqs))
        # all mode: may return None if no solution, otherwise return an iterator
        # one solution mode: return the solution directly if it exists, otherwise return None
        return m4ri_solve(eqs, cols, mode)

    def _convert_sol(self, s: int) -> tuple[int, ...]:
        sol = []
        for size in self._sizes:
            sol.append(s & ((1 << size) - 1))
            s >>= size
        assert s == 0, "Invalid solution"
        return tuple(sol)

    def convert_sol(self, s: int) -> Optional[tuple[int, ...]]:
        return self._convert_sol(s)

    def solve_space(self, zeros: Zeros) -> Optional[AffineSpace]:
        return self.solve_raw(zeros, 1)

    def solve_all(self, zeros: Zeros, max_dimension: int = 16):
        space = self.solve_space(zeros)
        if space is None:
            return
        if space.dimension > max_dimension:
            raise DimensionTooLargeError(
                f"Solution space (dim {space.dimension}) is too large, try increase max_dimension ({max_dimension}) if you want (there will be 2**dim solutions)",
                space=space,
            )
        for s in space:
            ret = self.convert_sol(s)
            if ret is not None:
                yield ret

    def solve_one(self, zeros: Zeros):
        sol = self.solve_raw(zeros, 0)
        if sol is None:
            return
        return self.convert_sol(sol)


class QuadraticSystem(LinearSystem):
    def __init__(self, sizes: list[int]):
        n = sum(sizes)
        quad_terms = n * (n - 1) // 2
        super().__init__(sizes + [quad_terms])
        self._lin_size = n
        self._const_lin_mask = (1 << (1 + n)) - 1
        self._quad_size = quad_terms

        # for pickle
        self._quad_init_sizes = sizes

    def gens(self):
        return super().gens()[:-1]

    def __reduce__(self):
        return (self.__class__, (self._quad_init_sizes,))

    def mul_bit_slow(self, a: int, b: int):
        # TODO: find a way to optimize this
        n = self._lin_size
        # constant term and linear terms (x^2 = x in GF(2))
        v = (a & self._const_lin_mask) & b
        # quadratic terms
        # a >>= 1
        # b >>= 1
        # abits = [0] * n
        # bbits = [0] * n
        # for i in range(n):
        #     abits[i] = a & 1
        #     bbits[i] = b & 1
        #     a >>= 1
        #     b >>= 1
        abits = to_bits(n, a >> 1)
        bbits = to_bits(n, b >> 1)
        # assert a == b == 0, "a, b is not linear terms"
        mi = 1 + n
        for i in range(n):
            for j in range(i):
                r = (abits[i] & bbits[j]) ^ (abits[j] & bbits[i])
                if r:
                    v |= self._basis[mi]  # same as (1 << mi)
                mi += 1
        # assert (v >> (1 + self._lin_size + self._quad_size)) == 0, "Overflow"
        return v

    def mul_bit(self, a: int, b: int) -> int:
        # constant term and linear terms (x^2 = x in GF(2))
        v = (a & self._const_lin_mask) & b
        # quadratic terms
        return mul_bit_quad(self._lin_size, a >> 1, b >> 1, v, self._basis)

    def bit_assert(self, a: int, v: int):
        # this is to assert `a` bit is exactly equal to `v`
        assert v in (0, 1), "Invalid bit"
        assert a not in (0, 1), "a should not be a constant"
        assert a >> self._lin_size == 0, "Not a linear term"
        zeros = [a ^ v]  # assert simple linear term
        # and it a linearized system, we can see that:
        # (bits[x] + bits[y] + ...) * bits[?] = 0 if v == 0
        # (bits[x] + bits[y] + ...) * bits[?] = bits[?] if v == 1
        # and `a` represents (bits[x] + bits[y] + ...), `b` represents bits[?]
        for i in range(1, 1 + self._lin_size):
            b = self._basis[i]  # linear term
            if a == b:
                continue
            if v == 0:
                zeros.append(self.mul_bit(a, b))
            else:
                zeros.append(self.mul_bit(a, b) ^ b)
        return zeros

    def _check_lin_match_quad(self, lin: int, quad: int):
        n = self._lin_size
        lin_bits = [0] * n
        for i in range(n):
            lin_bits[i] = lin & 1
            lin >>= 1
        assert lin == 0, "Invalid linear part"
        for i in range(n):
            for j in range(i):
                r = lin_bits[i] & lin_bits[j]
                if r != (quad & 1):
                    return False
                quad >>= 1
        assert quad == 0, "Invalid quadratic part"
        return True

    def convert_sol(self, s: int) -> Optional[tuple[int, ...]]:
        lin = s & ((1 << self._lin_size) - 1)
        s >>= self._lin_size
        quad = s & ((1 << self._quad_size) - 1)
        s >>= self._quad_size
        assert s == 0, "Invalid solution"
        if self._check_lin_match_quad(lin, quad):
            return super()._convert_sol(lin)[:-1]

    def solve_one(self, zeros: Zeros):
        # we can't use the LinearSystem.solve_one because the returned solution might not pass convert_sol
        for sol in self.solve_all(zeros):
            return sol
