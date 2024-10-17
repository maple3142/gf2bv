from setuptools import setup, Extension


setup(
    name="gf2bv",
    packages=["gf2bv", "gf2bv.crypto"],
    ext_modules=[
        Extension(
            "gf2bv.pym4ri",
            sources=["gf2bv/pym4ri.c"],
            libraries=["m4ri"],
            extra_compile_args=["-O3", "-march=native", "-mtune=native"],
        )
    ],
)
