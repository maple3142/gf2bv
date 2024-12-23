# Maintainer: maple3142 <kirby741852963@gmail.com>
pkgname=python-gf2bv
pkgver=0.0.1
pkgrel=1
pkgdesc="Solving linear systems over GF(2) by manipulating bitvectors"
arch=('x86_64')
url="https://github.com/maple3142/gf2bv"
license=('GPL')
depends=('m4ri')
makedepends=(python-build python-installer python-wheel)
source=("src::git+file://$startdir")
md5sums=('SKIP')
_name=${pkgname#python-}

build() {
    cd src
    python -m build --wheel --no-isolation
}

package() {
    cd src
    python -m installer --destdir="$pkgdir" dist/*.whl
}
