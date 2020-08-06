#/bin/sh

curl -L https://vtk.org/files/release/9.0/VTK-9.0.1.tar.gz | tar zxf -
cd VTK-9.0.1
mkdir build && cd build
cmake .. 
#-DVTK_GROUP_ENABLE_Imaging=WANT \
#  -DVTK_GROUP_ENABLE_MPI=DONT_WANT \
#  -DVTK_GROUP_ENABLE_Qt=DONT_WANT \
#  -DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
#  -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT \
#  -DVTK_GROUP_ENABLE_Views=DONT_WANT \
#  -DVTK_GROUP_ENABLE_Web=DONT_WANT
make && make install
