#/bin/bash

INSTALLER=cuda-repo-${UBUNTU_VERSION}_${CUDA}_amd64.deb
wget http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/x86_64/${INSTALLER}
sudo dpkg -i ${INSTALLER}
wget https://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
sudo apt update -qq
sudo apt install -y cuda-core-${CUDA_SHORT_VER} cuda-cudart-dev-${CUDA_SHORT_VER} cuda-cufft-dev-${CUDA_SHORT_VER}
sudo apt clean
CUDA_HOME=/usr/local/cuda-${CUDA_SHORT}
LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}
PATH=${CUDA_HOME}/bin:${PATH}
