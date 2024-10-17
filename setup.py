from setuptools import setup, Extension

module = Extension(
    "gf2bv.pym4ri",
    sources=["gf2bv/pym4ri.c"],
    libraries=["m4ri"],
    extra_compile_args=["-O3", "-march=native", "-mtune=native"],
)

setup(
    name="gf2bv",
    version="0.1.0",
    description="Solving linear systems over GF(2) by manipulating bitvectors",
    packages=["gf2bv", "gf2bv.crypto"],
    ext_modules=[module],
)
