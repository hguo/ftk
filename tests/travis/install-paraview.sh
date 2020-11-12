#/bin/sh

curl -L https://github.com/Kitware/ParaView/archive/v5.8.1.tar.gz | tar zxf -
cd ParaView-5.8.1
mkdir build && cd build

cmake .. -DPARAVIEW_BUILD_EDITION=CORE -DCMAKE_INSTALL_PREFIX=$HOME/paraview-5.8.1
make && make install
