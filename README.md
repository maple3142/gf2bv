# gf2bv - Solving linear systems over GF(2) by manipulating bitvectors

`gf2bv` allows you to build and solve linear systems over GF(2) by manipulating bitvectors, based on [M4RI](https://github.com/malb/m4ri).

## Installation

First, you need to have `m4ri` installed on your system, try some package names like `mr4i`, `libm4ri-dev` or something else depending on [your system](https://repology.org/project/libm4ri/versions).

Clone the repository and install the package with pip:

```bash
pip install .
```

Requires Python 3.11 or later.

## Usage

Define a linear system using `LinearSystem` and get the symbolic bitvectors with `gens()`, then use the them to build the equations you want to solve. Equations are represented as a list of symbolic bitvectors that evaluate to zero named `zeros` (You can choose other names if you want), then pass them to `solve_one` or `solve_all` to get the solutions.

### Find all solutions to a simple linear system

```python
from gf2bv import LinearSystem

lin = LinearSystem([1, 1, 1, 1])
a, b, c, d = lin.gens()

"""
This is the system we want to solve:
a + b + c = 1
b + d = 0
a + c = 1
"""

zeros = [a ^ b ^ c ^ 1, b ^ d, a ^ c ^ 1]
for sol in lin.solve_all(zeros):
    print(sol)
```

### Break a linear hash function

```python
from gf2bv import LinearSystem
import secrets

MASK64 = (1 << 64) - 1


def just_some_hash(x, y):
    z1 = x ^ (y >> 22)
    z2 = y ^ (x << 13) & MASK64
    z1 ^= (z2 >> 5) & 0xDEADBEEFDEADBEEF
    z2 ^= (z1 << 7) | 0x1337314213373142
    return z1, z2


x, y = secrets.randbits(64), secrets.randbits(64)
z1, z2 = just_some_hash(x, y)

lin = LinearSystem([64, 64])
xx, yy = lin.gens()
zz1, zz2 = just_some_hash(xx, yy)

zeros = [zz1 ^ z1, zz2 ^ z2]
sol = lin.solve_one(zeros)
assert sol == (x, y)
assert just_some_hash(*sol) == (z1, z2)
print(sol)
```

### More examples

Check the [examples](examples) directory for more examples.

## Prior work

This package is pretty similar to [xorsat](https://github.com/Lydxn/xorsat), but I keep most of the code in Python for flexibility and ease of use.
