import hashlib
import io
import os
import subprocess
import sys
import tarfile
from pathlib import Path
from urllib.request import urlopen

from setuptools import Extension, find_packages, setup


def download_and_build_m4ri():
    release = "20260122"
    checksum = bytes.fromhex(
        "68196ed43bc2f20f8cf84433ff5d7161bf71011f2f429ce3c70b372740f4d4cf"
    )
    workdir = Path(f"m4ri-{release}")
    libm4ri_a = workdir / ".libs" / "libm4ri.a"
    if libm4ri_a.exists():
        print(f"Using cached m4ri {release}")
    else:
        with urlopen(
            f"https://github.com/malb/m4ri/archive/refs/tags/{release}.tar.gz"
        ) as source:
            fobj = io.BytesIO(source.read())
            if hashlib.sha256(fobj.getbuffer()).digest() != checksum:
                raise ValueError("Checksum mismatch for downloaded m4ri")
            with tarfile.open(fileobj=fobj, mode="r|gz") as tar:
                tar.extractall()
        if not workdir.exists():
            raise FileNotFoundError(f"Failed to extract {workdir}")
        subprocess.run("autoreconf --install", shell=True, cwd=workdir, check=True)
        subprocess.run(
            './configure CFLAGS="-fPIC -O3 -march=native -mtune=native" --enable-openmp --enable-thread-safe',
            shell=True,
            cwd=workdir,
            check=True,
        )
        subprocess.run("make -j", shell=True, cwd=workdir, check=True)
        if not libm4ri_a.exists():
            raise FileNotFoundError(f"Failed to build {libm4ri_a}")
    extra_compile_args = ["-O3", "-march=native", "-mtune=native"]
    extra_link_args = []
    if sys.platform.startswith("darwin"):
        # you may want to use the following line if you installed libomp via Homebrew
        # export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
        extra_compile_args += ["-Xpreprocessor", "-fopenmp"]
        extra_link_args += ["-lomp"]
    elif not sys.platform.startswith("win"):
        extra_compile_args += ["-fopenmp"]
        extra_link_args += ["-fopenmp"]
    else:
        raise NotImplementedError("Windows is not supported yet")
    return Extension(
        "gf2bv._internal",
        sources=["gf2bv/_internal.c"],
        extra_link_args=extra_link_args,
        extra_compile_args=extra_compile_args,
        include_dirs=[str(workdir)],
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
