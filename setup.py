from setuptools import setup, Extension, find_packages


setup(
    name="gf2bv",
    packages=find_packages(),
    ext_modules=[
        Extension(
            "gf2bv._internal",
            sources=["gf2bv/_internal.c"],
            libraries=["m4ri"],
            extra_compile_args=["-O3", "-march=native", "-mtune=native"],
        )
    ],
)
