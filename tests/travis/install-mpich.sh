#/bin/sh

curl -L http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz | tar zxf -
cd mpich-3.3.2
./configure --disable-fortran
make && make install
