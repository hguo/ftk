#/bin/sh

curl -L https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_0/source/hdf5-1.12.0.tar.gz | tar zxf  -
cd hdf5-1.12.0
./configure 
make && make install

cd ..
curl -L https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz | tar zxf -
cd netcdf-c-4.7.4
./configure 
make && make install
