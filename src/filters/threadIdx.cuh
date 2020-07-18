#ifndef __THREADINDEX_CUH
#define __THREADINDEX_CUH

#define CUDA_SAFE_CALL(call) {\
  cudaError_t err = call;\
  if (cudaSuccess != err) {\
    fprintf(stderr, "[CUDA Error] %s, in file '%s', line%i.\n", cudaGetErrorString(err), __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }\
}

__device__ inline int getGlobalIdx_1D_1D()
{
	return blockIdx.x *blockDim.x + threadIdx.x;
}

__device__ inline int getGlobalIdx_1D_2D()
{
	return blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
}

__device__ inline int getGlobalIdx_1D_3D()
{
	return blockIdx.x * blockDim.x * blockDim.y * blockDim.z 
	 + threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
}

__device__ inline int getGlobalIdx_2D_1D()
{
	int blockId   = blockIdx.y * gridDim.x + blockIdx.x;			 	
	int threadId = blockId * blockDim.x + threadIdx.x; 
	return threadId;
}

 __device__ inline int getGlobalIdx_2D_2D()
{
	int blockId = blockIdx.x + blockIdx.y * gridDim.x; 
	int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
	return threadId;
}

__device__ inline int getGlobalIdx_2D_3D()
{
	int blockId = blockIdx.x 
			 + blockIdx.y * gridDim.x; 
	int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
			   + (threadIdx.z * (blockDim.x * blockDim.y))
			   + (threadIdx.y * blockDim.x)
			   + threadIdx.x;
	return threadId;
} 

__device__ inline int getGlobalIdx_3D_1D()
{
	int blockId = blockIdx.x 
			 + blockIdx.y * gridDim.x 
			 + gridDim.x * gridDim.y * blockIdx.z; 
	int threadId = blockId * blockDim.x + threadIdx.x;
	return threadId;
} 

__device__ inline int getGlobalIdx_3D_2D()
{
	int blockId = blockIdx.x 
		         + blockIdx.y * gridDim.x 
			 + gridDim.x * gridDim.y * blockIdx.z; 
	int threadId = blockId * (blockDim.x * blockDim.y)
			  + (threadIdx.y * blockDim.x)
			  + threadIdx.x;
	return threadId;
}

__device__ inline int getGlobalIdx_3D_3D()
{
	int blockId = blockIdx.x 
			 + blockIdx.y * gridDim.x 
			 + gridDim.x * gridDim.y * blockIdx.z; 
	int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
			  + (threadIdx.z * (blockDim.x * blockDim.y))
			  + (threadIdx.y * blockDim.x)
			  + threadIdx.x;
	return threadId;
}

#endif
