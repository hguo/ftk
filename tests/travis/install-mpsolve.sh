#/bin/sh

curl -L https://github.com/robol/MPSolve/archive/3.2.1.tar.gz | tar zxf - 
cd MPSolve-3.2.1
./autogen.sh
./configure --prefix=/tmp/MPSolve-3.2.1
make && make install
