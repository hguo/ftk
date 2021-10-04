#ifndef _BVH2D_HH
#define _BVH2D_HH

template <typename I=int, typename F=double>
struct bvh2d_node_t {
  // tree
  I parentId;
  I childrenIds[4];

  // bounds
  F Ax, Ay, Bx, By;

  // triangle
  I triangleId; // -1 if the node if not leaf
  I i0, i1, i2;
  F x0, y0, x1, y1, x2, y2;
};


#endif
