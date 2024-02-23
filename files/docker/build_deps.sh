#!/usr/bin/env bash

set -euo pipefail

# Build binary dependencies for PAV:
# htslib
# samtools
# tabix
# minimap2
# LRA
# ucsc (command-line tools)


PREFIX=/opt/pav

HTSLIB_VERSION=1.19

SAMTOOLS_VERSION=${HTSLIB_VERSION}

MINIMAP_VERSION=2.26

LRA_VERSION=1.3.7.2


### htslib ###

wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2

tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2

pushd htslib-${HTSLIB_VERSION}

./configure --prefix=${PREFIX} CPPFLAGS="-I${PREFIX}/include" LDFLAGS="-L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib"

make -j 4 && make install

popd

rm -r htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2

rm bin/htsfile


### samtools ###

wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2

pushd samtools-${SAMTOOLS_VERSION}

./configure --prefix=${PREFIX} CPPFLAGS="-I${PREFIX}/include" LDFLAGS="-L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib"

make -j 4 && make install

popd

rm -r samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

rm bin/{ace2sam,wgsim} bin/{plot-,maq2sam}* bin/*.pl


### minimap2 ###

wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP_VERSION}/minimap2-${MINIMAP_VERSION}.tar.bz2

tar -xjf minimap2-${MINIMAP_VERSION}.tar.bz2

make -j 4 -C minimap2-${MINIMAP_VERSION}

install minimap2-${MINIMAP_VERSION}/minimap2 ${PREFIX}/bin/

rm -r minimap2-${MINIMAP_VERSION} minimap2-${MINIMAP_VERSION}.tar.bz2


### LRA ###

git clone --recursive https://github.com/ChaissonLab/LRA.git -b v${LRA_VERSION}

pushd LRA

make -j 4 CFLAGS="-I${PREFIX}/include -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib -g"

install lra ${PREFIX}/bin/

popd

rm -r LRA


### UCSC ###

curl http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed --output bin/bedToBigBed
chmod ugo+x bin/bedToBigBed


### Cleanup ###

rm -r share

rm lib/libhts.a
