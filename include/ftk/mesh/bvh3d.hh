#ifndef _BVH3D_HH
#define _BVH3D_HH

template <typename I=int, typename F=double>
struct bvh3d_node_t {
  // tree
  I parentId;
  I childrenIds[8];

  // bounds
  F Ax, Ay, Az, Bx, By, Bz;

  // tet
  I tetId; // -1 if the node if not leaf
  I i0, i1, i2, i3;
  F x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
};


#endif
