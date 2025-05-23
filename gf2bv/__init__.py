from __future__ import annotations

from collections.abc import Sequence
from functools import reduce
from operator import xor
from typing import Literal, Optional, TypeVar

from ._internal import (
    AffineSpace,
    eqs_to_sage_mat_helper,
    m4ri_solve,
    mul_bit_quad,
    to_bits,
    tuple_where,
    xor_tuple,
)

TSolveMode = TypeVar("TSolveMode", Literal[0], Literal[1])


class BitVec:
    __slots__ = ("_bits",)

    def __init__(self, bits: tuple[int, ...]):
        # this represent the symbolic bits of the BitVec
        # in little-endian order (lsb first)
        self._bits = bits

    def __len__(self):
        return len(self._bits)

    def __getitem__(self, key: int | slice):
        if isinstance(key, slice):
            return BitVec(self._bits[key])
        return BitVec(
            (self._bits[key],)
        )  # still wrap it to prevent misuse like bv[0] ^ bv

    def __xor__(self, other: BitVec | int):
        if not isinstance(other, BitVec):
            bs = to_bits(len(self._bits), other)
            return BitVec(xor_tuple(self._bits, bs))
        else:
            if len(self._bits) != len(other._bits):
                raise ValueError("Cannot mix bitvecs of different lengths")
        return BitVec(xor_tuple(self._bits, other._bits))

    __rxor__ = __xor__
    __pow__ = __xor__  # alias to __xor__, for convenience in sage

    def __rshift__(self, n: int):
        if n == 0:
            return self
        return BitVec(self._bits[n:] + (0,) * n)

    def __lshift__(self, n: int):
        if n == 0:
            return self
        return BitVec((0,) * n + self._bits[:-n])

    def lshift_ext(self, n: int):
        return BitVec((0,) * n + self._bits)

    def __and__(self, mask: int):
        bs = to_bits(len(self._bits), mask)
        if all(bs):
            # if all bits are set, it does not change anything
            return self
        return BitVec(tuple_where(bs, self._bits, 0))

    __rand__ = __and__

    def __or__(self, mask: BitVec | int):
        if isinstance(mask, BitVec):
            if len(self._bits) > len(mask._bits):
                self, mask = mask, self
            ar = [0] * len(mask._bits)
            for i in range(len(self._bits)):
                if self._bits[i] not in (0, 1) and mask._bits[i] not in (0, 1):
                    raise ValueError(
                        "Cannot compute logical or using bitvecs with non-zero bits"
                    )
                if self._bits[i] == 1 or mask._bits[i] == 1:
                    ar[i] = 1
                elif self._bits[i] == 0:
                    ar[i] = mask._bits[i]
                elif mask._bits[i] == 0:
                    ar[i] = self._bits[i]
            ar[len(self._bits) :] = mask._bits[len(self._bits) :]
            return BitVec(tuple(ar))
        bs = to_bits(len(self._bits), mask)
        if all(bs):
            # if all bits are set, it becomes all ones
            return BitVec(bs)
        return BitVec(tuple_where(bs, 1, self._bits))

    __ror__ = __or__

    def __mod__(self, n: int):
        if n & (n - 1) != 0:
            raise ValueError("modulo non-power-of-2 is not a linear operation")
        return self & (n - 1)

    def rotr(self, n: int):
        return BitVec(self._bits[n:] + self._bits[:n])

    def rotl(self, n: int):
        return BitVec(self._bits[-n:] + self._bits[:-n])

    def sum(self):
        return BitVec((reduce(xor, self._bits),))

    def zeroext(self, n: int):
        return BitVec(self._bits + (0,) * n)

    def signext(self, n: int):
        return BitVec(self._bits + (self._bits[-1],) * n)

    def broadcast(self, i: int, n: int):
        return BitVec((self._bits[i],) * n)

    def dup(self, n: int):
        return BitVec(self._bits * n)

    def concat(self, other: BitVec):
        return BitVec(self._bits + other._bits)

    def evaluate(self, s: int):
        """
        Evaluate the BitVec using the given raw solution
        """
        r1 = (s << 1) | 1
        bs = ((b & r1).bit_count() & 1 for b in reversed(self._bits))
        return int("".join(map(str, bs)), 2)


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
            _vars.append(BitVec(tuple(self._basis[i : i + size])))
            i += size
        self._vars = tuple(_vars)

    def gens(self):
        return self._vars

    def __reduce__(self):
        return (self.__class__, (self._sizes,))

    def get_sage_mat_slow(self, zeros: Zeros, *, tqdm=lambda x, desc: x):
        """
        Convert the system of equations to Sage, return a matrix A and a vector b such that Ax = b
        """
        from sage.all import GF, matrix, vector

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
                    mat[i, j] = bs[j]
            i += 1
        return mat, affine

    def get_sage_mat(self, zeros: Zeros):
        """
        Convert the system of equations to Sage, return a matrix A and a vector b such that Ax = b
        """
        import struct

        from sage.all import GF, vector
        from sage.matrix.matrix_mod2_dense import unpickle_matrix_mod2_dense_v2

        F2 = GF(2)
        eqs = self.get_eqs(zeros)
        cols = self._cols
        rows = len(eqs)
        buf, affine = eqs_to_sage_mat_helper(eqs, cols)
        affine = vector(F2, affine)
        # convert it to signed bytes and deallocate the original buffer
        buf = struct.unpack(f">{len(buf)}b", buf)  # type: ignore
        mat = unpickle_matrix_mod2_dense_v2(rows, cols, buf, len(buf), False)
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
                (bv._bits if isinstance(bv, BitVec) else [bv] for bv in zeros),
            )
        )
        return list(filter(None, eqs))  # remove literal zeros

    def _solve_internal(self, zeros: Zeros, mode: TSolveMode):
        eqs = self.get_eqs(zeros)
        if 1 in eqs:
            # no solution
            return
        cols = self._cols
        if cols > len(eqs):
            # m4ri_solve requires rows >= cols, pad with zeros
            eqs += [0] * (cols - len(eqs))
        # all mode (mode == 1): may return None if no solution, otherwise return an iterator
        # one solution mode (mode == 0): return the solution directly if it exists, otherwise return None
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

    def solve_raw_one(self, zeros: Zeros) -> Optional[AffineSpace]:
        return self._solve_internal(zeros, 0)

    def solve_raw_space(self, zeros: Zeros) -> Optional[AffineSpace]:
        return self._solve_internal(zeros, 1)

    def solve_all(self, zeros: Zeros, *, max_dimension: int = 16):
        space = self.solve_raw_space(zeros)
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
        sol = self._solve_internal(zeros, 0)
        if sol is None:
            return
        return self.convert_sol(sol)

    def evaluate(self, bv: BitVec, sol: tuple[int, ...]) -> int:
        """
        Evaluate the BitVec using the given solution
        """
        s = 0
        for v, sz in zip(reversed(sol), reversed(self._sizes)):
            s <<= sz
            s |= v
        return bv.evaluate(s)


class QuadraticSystem(LinearSystem):
    def __init__(self, sizes: list[int]):
        n = sum(sizes)
        quad_terms = n * (n - 1) // 2
        super().__init__(sizes + [quad_terms])
        self._quad_sizes = sizes[:]
        self._lin_size = n
        self._const_lin_mask = (1 << (1 + n)) - 1
        self._quad_size = quad_terms

    def gens(self):
        return super().gens()[:-1]

    def __reduce__(self):
        return (self.__class__, (self._quad_sizes,))

    def _mul_bit_slow(self, a: int, b: int):
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

    def _mul_bit(self, a: int, b: int) -> int:
        # constant term and linear terms (x^2 = x in GF(2))
        v = (a & self._const_lin_mask) & b
        # quadratic terms
        return mul_bit_quad(self._lin_size, a >> 1, b >> 1, v, self._basis)

    def mul_bit(self, a: BitVec, b: BitVec) -> BitVec:
        if len(a) != 1 or len(b) != 1:
            raise ValueError("The inputs should be single bits")
        return BitVec((self._mul_bit(a._bits[0], b._bits[0]),))

    def _bit_assert(self, a: int, v: int):
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
                zeros.append(self._mul_bit(a, b))
            else:
                zeros.append(self._mul_bit(a, b) ^ b)
        return zeros

    def bit_assert(self, a: BitVec, v: int) -> Zeros:
        if len(a) != 1:
            raise ValueError("The input should be a single bit")
        return self._bit_assert(a._bits[0], v)

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

    def evaluate(self, bv: BitVec, sol: tuple[int, ...]) -> int:
        """
        Evaluate the BitVec using the given solution
        """
        s = 0
        for v, sz in zip(reversed(sol), reversed(self._quad_sizes)):
            s <<= sz
            s |= v
        return bv.evaluate(s)
