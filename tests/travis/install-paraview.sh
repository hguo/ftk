#/bin/sh

curl https://www.paraview.org/paraview-downloads/download.php\?submit\=Download\&version\=v5.8\&type\=source\&os\=Sources\&downloadFile\=ParaView-v5.8.1.tar.gz -o paraview.tar.gz
tar zxf paraview.tar.gz
# curl -L https://github.com/Kitware/ParaView/archive/v5.8.1.tar.gz | tar zxf -
cd ParaView-v5.8.1
mkdir build && cd build

cmake .. -DPARAVIEW_BUILD_EDITION=CORE -DCMAKE_INSTALL_PREFIX=$HOME/paraview-5.8.1
make && make install
