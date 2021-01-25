#ifndef _FTK_M2E_CUH
#define _FTK_M2E_CUH

#include "m2.cuh"

template <typename I>
__device__
inline I mod(I v, I m)
{
  I x = v % m;
  if (v < 0) x += m;
  return x;
}

template <typename I, typename F>
__device__
inline void me2_get_coords(
    I v3, F x[], 
    const I m2n0, const F m2coords[])
{
  const I v2 = mod(v3, m2n0);
  x[0] = m2coords[v2*2];
  x[1] = m2coords[v2*2+1];
  x[2] = v3 / m2n0;
}

template <typename I>
__device__
int me2_edge_type(I i, const I m2n1)
{
  if (i < m2n1) return 0;
  else if (i < 2 * m2n1) return 1;
  else return 2;
}

template <typename I>
__device__
int me2_tri_type(I i, const I m2n1, const I m2n2)
{
  if (i < m2n2) return 0;
  else if (i < 2 * m2n2) return 1;
  else if (i < 3 * m2n2) return 2; 
  else if (i < 3 * m2n2 + m2n1) return 3;
  else return 4;
}

template <typename I>
__device__
void me2_get_edge(I k, I verts[2], const I m2n0, const I m2n1, 
    const I me2edges[])
{
  const I me2n1 = 2 * m2n1 + m2n0;
  const I i = mod(k, me2n1), t = floor(double(k) / me2n1);
  const int type = me2_edge_type(i, m2n1);
  const I offset = t * m2n0;

  if (type < 2) {
    I m2edge[2];
    m2_get_edge(i % m2n1, m2edge, me2edges);
    if (type == 0) {
      verts[0] = m2edge[0];
      verts[1] = m2edge[1];
    } else { // type == 1
      verts[0] = m2edge[0];
      verts[1] = m2edge[1] + m2n0;
    }
  } else {
    verts[0] = i - 2 * m2n1;
    verts[1] = i - 2 * m2n1 + m2n0;
  }
  
  for (int j = 0; j < 2; j ++)
    verts[j] += offset;
}

template <typename I>
__device__
void me2_get_tri(I k, I verts[3], 
    const I m2n0, const I m2n1, const I m2n2,
    const I m2edges[], const I m2tris[])
{
  // call m2_get_tri, m2_get_edge
  const I me2n2 = 3 * m2n2 + 2 * m2n1;
  const I i = mod(k, me2n2), t = floor(double(k) / me2n2);
  const int type = me2_tri_type(i, m2n1, m2n2);
  const I offset = t * m2n0;

  if (type < 3) {
    I m2tri[3];
    m2_get_tri(i % m2n2, m2tri, m2tris);
    if (type == 0) {
      verts[0] = m2tri[0]; 
      verts[1] = m2tri[1];
      verts[2] = m2tri[2];
    } else if (type == 1) {
      verts[0] = m2tri[0];
      verts[1] = m2tri[1];
      verts[2] = m2tri[2] + m2n0;
    } else if (type == 2) {
      verts[0] = m2tri[0];
      verts[1] = m2tri[1] + m2n0;
      verts[2] = m2tri[2] + m2n0;
    }
  } else {
    I m2edge[2];
    m2_get_edge((i - 3 * m2n2) % m2n1, m2edge, m2edges);
    if (type == 3) {
      verts[0] = m2edge[0];
      verts[1] = m2edge[1]; 
      verts[2] = m2edge[1] + m2n0;
    } else { // type 4
      verts[0] = m2edge[0];
      verts[1] = m2edge[0] + m2n0;
      verts[2] = m2edge[1] + m2n0;
    }
  }

  for (int j = 0; j < 3; j ++)
    verts[j] += offset;
}

#endif
