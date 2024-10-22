from setuptools import setup, Extension, find_packages
import os, tarfile, subprocess
from urllib.request import urlopen
from pathlib import Path


def download_and_build_m4ri():
    release = "20240729"
    workdir = Path(f"m4ri-release-{release}")
    libm4ri_a = workdir / ".libs" / "libm4ri.a"
    if libm4ri_a.exists():
        print(f"Using cached m4ri {release}")
    else:
        with urlopen(
            f"https://github.com/malb/m4ri/archive/refs/tags/release-{release}.tar.gz"
        ) as source:
            with tarfile.open(fileobj=source, mode="r|gz") as tar:
                tar.extractall()
        if not workdir.exists():
            raise FileNotFoundError(f"Failed to extract {workdir}")
        subprocess.run("autoreconf --install", shell=True, cwd=workdir, check=True)
        subprocess.run(
            f'./configure CFLAGS="-fPIC -O3 -march=native -mtune=native" --enable-openmp --enable-thread-safe',
            shell=True,
            cwd=workdir,
            check=True,
        )
        subprocess.run("make -j", shell=True, cwd=workdir, check=True)
        if not libm4ri_a.exists():
            raise FileNotFoundError(f"Failed to build {libm4ri_a}")
    return Extension(
        "gf2bv._internal",
        sources=["gf2bv/_internal.c"],
        libraries=["gomp"],
        extra_compile_args=["-O3", "-march=native", "-mtune=native"],
        library_dirs=[str(workdir)],
        extra_objects=[str(libm4ri_a)],
    )


if os.environ.get("GF2BV_BUILD_M4RI", "0") != "0":
    ext = download_and_build_m4ri()
else:
    ext = Extension(
        "gf2bv._internal",
        sources=["gf2bv/_internal.c"],
        libraries=["m4ri"],
        extra_compile_args=["-O3", "-march=native", "-mtune=native"],
    )

setup(
    name="gf2bv",
    packages=find_packages(),
    ext_modules=[ext],
)
