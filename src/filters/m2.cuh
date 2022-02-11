#ifndef _FTK_MESH2_CUH
#define _FTK_MESH2_CUH

template <typename I>
__device__
inline void m2_get_tri(I i, I tri[3], const I m2tris[])
{
  for (int k = 0; k < 3; k ++)
    tri[k] = m2tris[3*i+k];
}

template <typename I>
__device__
inline void m2_get_edge(I i, I edge[2], const I m2edges[])
{
  edge[0] = m2edges[2*i];
  edge[1] = m2edges[2*i+1];
}

#endif
