#ifndef _FTK_LOCATOR2_HH
#define _FTK_LOCATOR2_HH

#include <ftk/mesh/bvh2d.hh>

template <typename I, typename F>
__device__ __host__
inline bool bvh2_inside_quad(const bvh2d_node_t<I, F> &q, F x, F y)
{
#if __CUDA_ARCH__
  return __fmul_rz(x-q.Ax, x-q.Bx) <= 0 && __fmul_rz(y-q.Ay, y-q.By) <= 0;
#else
  return (x - q.Ax) * (x - q.Bx) <= 0 && (y - q.Ay) * (y - q.By) <= 0;
#endif
  // return x >= q.Ax && x < q.Bx && y >= q.Ay && y < q.By;
}

template <typename I, typename F>
__device__ __host__
inline bool bvh2_inside_triangle(const bvh2d_node_t<I, F> &q, F x, F y, F lambda[3], const F *invdet)
{
#if 0
  lambda.x = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) /
          ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  lambda.y = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) /
         ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
#endif
  const F d = invdet[q.triangleId];
  lambda[0] = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) * d; 
  lambda[1] = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) * d;
  lambda[2] = 1.0 - lambda[0] - lambda[1];
  // fprintf(stderr, "barycentric=%f, %f, %f\n", lambda.x, lambda.y, lambda.z);
  return lambda[0] >= 0 && lambda[1] >= 0 && lambda[2] >= 0 && 
         lambda[0] <= 1 && lambda[1] <= 1 && lambda[2] <= 1;
}

template <typename I, typename F>
__device__ __host__
inline I bvh2_locate_point_recursive(const bvh2d_node_t<I, F> *q, const bvh2d_node_t<I, F> *nodes, F x, F y, F lambda[3], const F *invdet)
{
  if (q->triangleId >= 0) { //leaf node
    bool succ = bvh2_inside_triangle(*q, x, y, lambda, invdet);
    if (succ) return q->triangleId;
  } else if (bvh2_inside_quad(*q, x, y)) {
    for (int j=0; j<4; j++) {
      if (q->childrenIds[j] > 0) {
        I result = bvh2_locate_point_recursive(&nodes[q->childrenIds[j]], nodes, x, y, lambda, invdet);
        if (result >= 0) return result;
      }
    }
  }
  return -1;
}

template <typename I, typename F>
__device__ __host__
inline I bvh2_locate_point(const bvh2d_node_t<I, F> *nodes, F x, F y, F lambda[3], const F *invdet, int root=0)
{
  // float lambda.x, lambda.y, lambda.z;
  static const int maxStackSize = 64;
  int stack[maxStackSize];
  int stackPos = 0;
  stack[stackPos++] = root; // push root

  while (stackPos > 0) {
    const int i = stack[--stackPos]; // pop
    const bvh2d_node_t<I, F> &q = nodes[i];

    // fprintf(stderr, "D_checking node %d, %f, %f, %f, %f\n", i, q.Ax, q.Ay, q.Bx, q.By);
    // fprintf(stderr, "D_checking node %d\n", i);

    if (q.triangleId >= 0) { // leaf node
      bool succ = bvh2_inside_triangle(q, x, y, lambda, invdet);
      if (succ) return i; // q.triangleId;
    } else if (bvh2_inside_quad(q, x, y)) { // non-leaf node
      for (int j=0; j<4; j++) {
        if (q.childrenIds[j] > 0)
          stack[stackPos++] = q.childrenIds[j];
      }
    }
  }
  return -1;
}

template <typename I, typename F>
__device__ __host__
inline I bvh2_locate_point_tri(const bvh2d_node_t<I, F> *nodes, F x, F y, F lambda[3], const F *invdet, int root=0)
{
  I nid = bvh2_locate_point(nodes, x, y, lambda, invdet, root);
  if (nid < 0) return -1;
  else return nodes[nid].triangleId;
}

template <typename I, typename F>
__device__ __host__
inline int bvh2_locate_point_coherent(const bvh2d_node_t<I, F> *bvh, int last_nid, F x, F y, F lambda[3], const F *invdet, const I *neighbors)
{
  // check if last_nid is valid
  if (last_nid<0) return bvh2_locate_point(bvh, x, y, lambda, invdet);

  // check if in the same triangle
  if (bvh2_inside_triangle(bvh[last_nid], x, y, lambda, invdet)) return last_nid;

  // check if neighbor triangles have the point
#if 0
  for (int i=0; i<3; i++) {
    int triangleId = bvh[last_nid].triangleId;
    int neighborQuadId = neighbors[triangleId*3];
    if (neighborQuadId<0) continue;
    else if (bvh2_inside_triangle(bvh[neighborQuadId], x, y, lambda, invdet)) return neighborQuadId;
  }
#endif

  // traverse from parents
  // int nid = bvh2_locate_point(bvh, x, y, lambda, invdet, bvh[bvh[last_nid].parentId].parentId);
  // int nid = bvh2_locate_point(bvh, x, y, lambda, invdet, bvh[last_nid].parentId);
  // if (nid >= 0) return nid;

  // TODO: check if in triangle neighbors of last_nid

  // fallback
  return bvh2_locate_point(bvh, x, y, lambda, invdet);
}

// template <typename I, typename F>
// __device__ __host__
// inline F bvh2_sample(int i0, int i1, int i2, const F lambda[3], const F *data) {
//   return lambda.x * data[i0] + lambda.y * data[i1] + lambda.z * data[i2];
// }

// template <typename I, typename F>
// __device__ __host__
// inline float2 bvh2_sample2(int i0, int i1, int i2, const F lambda[3], const F *data) {
//   return make_float2(lambda.x * data[i0*2] + lambda.y * data[i1*2] + lambda.z * data[i2*2],
//       lambda.x * data[i0*2+1] + lambda.y * data[i1*2+1] + lambda.z * data[i2*2+1]);
// }

// template <typename I, typename F>
// __device__ __host__
// inline F bvh2_sample(bvh2d_node_t<I, F>* bvh, int nid, const F lambda[3], const F *data) {
//   const bvh2d_node_t<I, F> &q = bvh[nid];
//   return lambda.x * data[q.i0] + lambda.y * data[q.i1] + lambda.z * data[q.i2];
// }

// template <typename I, typename F>
// __device__ __host__
// inline float2 bvh2_sample2(bvh2d_node_t<I, F>* bvh, int nid, float3 lambda, float *data) {
//   const bvh2d_node_t<I, F> &q = bvh[nid];
//   return bvh2_sample2(q.i0, q.i1, q.i2, lambda, data);
//   // return make_float2(lambda.x * data[q.i0*2] + lambda.y * data[q.i1*2] + lambda.z * data[q.i2*2],
//   //     lambda.x * data[q.i0*2+1] + lambda.y * data[q.i1*2+1] + lambda.z * data[q.i2*2+1]);
// }

#endif
