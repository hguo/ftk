#ifndef _FTK_MX4_CUH
#define _FTK_MX4_CUH

#include "mx3.cuh"

template <typename I>
__device__
int mx4_tri_type(const I i, const I mx3n1, const I mx3n2)
{
  // const I i = mod(k, mx4n2); // , t = floor((double)k / mx4n2);
  if (i < mx3n2) return 0;
  else if (i < 2 * mx3n2) return 1;
  else if (i < 3 * mx3n2) return 2;
  else if (i < 3 * mx3n2 + mx3n1) return 3;
  else return 4;
}

template <typename I>
__device__
void mx4_get_tri(I k, I verts[3], 
    const I np,
    const I m2n0, const I m2n1, const I m2n2, 
    const I m2edges[], const I m2tris[])
{
  const I mx3n0 = m2n0 * np; 
  const I mx3n1 = (2 * m2n1 + m2n0) * np; 
  const I mx3n2 = (3 * m2n2 + 2 * m2n1) * np;

  const I mx4n2 = 3 * mx3n2 + 2 * mx3n1; 

  const I i = mod(k, mx4n2), t = floor(double(k) / mx4n2);
  const int type = mx4_tri_type(i, mx3n1, mx3n2);
  const I offset = t * mx3n0;

  // const I n2 = 3 * m2n2 + 2 * m2n1;
  // const I i = mod(k, n2), t = std::floor(double(k) / n2);
  // const int type = mx4_tri_type(i, m2n1, m2n2);

  if (type < 3) {
    I mx3tri[3];
    mx3_get_tri(i % mx3n2, mx3tri, np, m2n0, m2n1, m2n2, m2edges, m2tris);
    if (type == 0) {
      verts[0] = mx3tri[0]; 
      verts[1] = mx3tri[1];
      verts[2] = mx3tri[2];
    } else if (type == 1) {
      verts[0] = mx3tri[0];
      verts[1] = mx3tri[1];
      verts[2] = mx3tri[2] + mx3n0;
    } else if (type == 2) {
      verts[0] = mx3tri[0];
      verts[1] = mx3tri[1] + mx3n0;
      verts[2] = mx3tri[2] + mx3n0;
    } else 
      assert(false);
  } else {
    I mx3edge[2];
    mx3_get_edge((i - 3 * mx3n2) % mx3n1, mx3edge, np, m2n0, m2n1, m2edges);
    if (type == 3) {
      verts[0] = mx3edge[0];
      verts[1] = mx3edge[1]; 
      verts[2] = mx3edge[1] + mx3n0;
    } else { // type 4
      verts[0] = mx3edge[0];
      verts[1] = mx3edge[0] + mx3n0;
      verts[2] = mx3edge[1] + mx3n0;
    }
  }

  for (int j = 0; j < 3; j ++)
    verts[j] += offset;
}

#endif
