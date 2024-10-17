from typing import Union
from .pym4ri import solve


class BitVec:
    def __init__(self, sys: "LinearSystem", st: list[int]):
        self.sys = sys
        self.st = st

    def _cast_int(self, n: int):
        if n < 0:
            raise ValueError("Negative value")
        st = []
        for i in range(len(self.st)):
            if n & 1:
                st.append(1)
            else:
                st.append(0)
            n >>= 1
        return BitVec(self.sys, st)

    def _check_other(self, other: "BitVec"):
        if self.sys is not other.sys:
            raise ValueError("Cannot mix bitvecs from different systems")
        if len(self.st) != len(other.st):
            raise ValueError("Cannot mix bitvecs of different lengths")

    def __xor__(self, other: Union["BitVec", int]):
        if not isinstance(other, BitVec):
            return self ^ self._cast_int(other)
        self._check_other(other)
        return BitVec(self.sys, [a ^ b for a, b in zip(self.st, other.st)])

    __rxor__ = __xor__

    def __pow__(self, other: Union["BitVec", int]):
        # Alias to __xor__, for convenience in sage
        return self ^ other

    def __rshift__(self, n: int):
        return BitVec(self.sys, self.st[n:] + [0] * n)

    def __lshift__(self, n: int):
        return BitVec(self.sys, [0] * n + self.st[:-n])

    def __and__(self, mask: int):
        if mask < 0:
            raise ValueError("Negative value")
        vecs = self.st[:]
        for i in range(len(self.st)):
            if mask & 1 == 0:
                vecs[i] = 0
            mask >>= 1
        return BitVec(self.sys, vecs)

    __rand__ = __and__

    def rotr(self, n: int):
        return BitVec(self.sys, self.st[n:] + self.st[:n])

    def rotl(self, n: int):
        return BitVec(self.sys, self.st[-n:] + self.st[:-n])

    def copy(self, n: int, m: int):
        """
        Construct a new bitvec by copying the n-th element m times
        """
        return BitVec(self.sys, [self.st[n]] * m)


class LinearSystem:
    def __init__(self, sizes: list[int]):
        self._sizes = sizes[:]
        self._vars = []
        i = 1  # lsb used to represent constant term (affine part)
        for size in self._sizes:
            # vecs = list(mat[i : i + size])
            st = [1 << j for j in range(i, i + size)]
            i += size
            self._vars.append(BitVec(self, st))
        self._vars = tuple(self._vars)

    def gens(self):
        return self._vars

    def get_sage_mat(self, zeros: list[BitVec], tqdm=lambda x, desc: x):
        """
        Convert the system of equations to Sage, return a matrix A and a vector b such that Ax = b
        """
        from sage.all import GF, vector, MatrixSpace
        from sage.matrix.matrix_mod2_dense import Matrix_mod2_dense

        F2 = GF(2)
        eqs = self.get_eqs(zeros)
        affine = []
        for i in range(len(eqs)):
            affine.append(eqs[i] & 1)
            eqs[i] >>= 1
        affine = vector(F2, affine)
        cols = sum(self._sizes)
        rows = len(eqs)
        mat = Matrix_mod2_dense(MatrixSpace(F2, rows, cols))
        i = 0
        for v in tqdm(eqs, desc="Converting equations"):
            for j in range(cols):
                if v & 1:
                    mat[i, j] = 1
                v >>= 1
            i += 1
        return mat, affine

    def get_eqs(self, zeros: list[BitVec]):
        """
        Convert zeros into a list of equations accepted by pym4ri.solve
        """
        # https://stackoverflow.com/questions/716477/join-list-of-lists-in-python/21034265#21034265
        # this seems to be the fastest way to flatten a list of lists
        eqs: list[int] = []
        list(map(eqs.extend, (bv.st for bv in zeros)))
        eqs = list(filter(None, eqs))  # remove literal zeros
        return eqs

    def solve_raw(self, zeros: list[BitVec]):
        eqs = self.get_eqs(zeros)
        if 1 in eqs:
            # no solution
            return
        cols = sum(self._sizes)
        if cols > len(eqs):
            # pym4ri.solve requires rows >= cols, pad with zeros
            eqs += [0] * (cols - len(eqs))
        # may return None if no solution, otherwise return an iterator
        return solve(eqs, cols, True)

    def convert_sol(self, s: int):
        sol = []
        for size in self._sizes:
            sol.append(s & ((1 << size) - 1))
            s >>= size
        assert s == 0
        return tuple(sol)

    def solve(self, zero_bvs: list[BitVec]):
        it = self.solve_raw(zero_bvs)
        if it is None:
            return
        for s in it:
            yield self.convert_sol(s)
