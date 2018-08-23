#ifndef _FTK_MESH
#define _FTK_MESH

#include <set>

namespace ftk {

template <class IdType>
std::set<size_t> Get4Neighbors2DRegular(IdType w, IdType h, IdType id) {
  auto nid2nidx = [w](IdType id, IdType &i, IdType &j) {
    j = id / w;
    i = id - j * w;
  };

  auto nidx2nid = [w](IdType i, IdType j) {
    return i + w * j;
  };

  IdType i, j;
  nid2nidx(id, i, j);

  std::set<IdType> neighbors;
  if (i > 0) neighbors.insert(nidx2nid(i-1, j));
  if (i < w-1) neighbors.insert(nidx2nid(i+1, j));
  if (j > 0) neighbors.insert(nidx2nid(i, j-1));
  if (j < h-1) neighbors.insert(nidx2nid(i, j+1));

  return neighbors;
}

template <class IdType>
std::set<size_t> Get8Neighbors2DRegular(IdType w, IdType h, IdType id) {
  auto nid2nidx = [w](IdType id, IdType &i, IdType &j) {
    j = id / w;
    i = id - j * w;
  };

  auto nidx2nid = [w](IdType i, IdType j) {
    return i + w * j;
  };
  
  IdType i, j;
  nid2nidx(id, i, j);

  IdType i0 = (i > 0)   ? (i - 1) : 0, 
         i1 = (i < w-1) ? (i + 1) : w - 1, 
         j0 = (j > 0)   ? (j - 1) : 0,
         j1 = (j < h-1) ? (j + 1) : h - 1;

  std::set<IdType> neighbors;
  neighbors.insert(nidx2nid(i0, j0));
  neighbors.insert(nidx2nid(i,  j0));
  neighbors.insert(nidx2nid(i1, j0));
  neighbors.insert(nidx2nid(i0, j ));
  
  neighbors.insert(nidx2nid(i1, j ));
  neighbors.insert(nidx2nid(i0, j1));
  neighbors.insert(nidx2nid(i,  j1));
  neighbors.insert(nidx2nid(i1, j1));

  return neighbors;
}

template <class IdType>
std::set<size_t> Get6Neighbors3DRegular(IdType w, IdType h, IdType d, IdType id) 
{
  const IdType s = w * h;

  auto nid2nidx = [w, s](IdType id, IdType &i, IdType &j, IdType &k) {
    k = id / s;
    j = (id - k*s) / w;
    i = id - k*s - j*w;
  };

  auto nidx2nid = [w, s](IdType i, IdType j, IdType k) {
    return i + j*w + k*s;
  };
 
  IdType i, j, k;
  nid2nidx(id, i, j, k);
  IdType i0 = (i > 0)   ? (i - 1) : 0, 
         i1 = (i < w-1) ? (i + 1) : w - 1, 
         j0 = (j > 0)   ? (j - 1) : 0,
         j1 = (j < h-1) ? (j + 1) : h - 1, 
         k0 = (k > 0)   ? (k - 1) : 0,
         k1 = (k < d-1) ? (k + 1) : d - 1;
  
  std::set<IdType> nb;
  nb.insert(nidx2nid(i,  j,  k0));
  nb.insert(nidx2nid(i,  j,  k1));
  nb.insert(nidx2nid(i0, j,  k ));
  nb.insert(nidx2nid(i1, j,  k ));
  nb.insert(nidx2nid(i,  j0, k ));
  nb.insert(nidx2nid(i,  j1, k ));

  return nb;
}

template <class IdType>
std::set<size_t> Get18Neighbors3DRegular(IdType W, IdType H, IdType D, IdType id);

template <class IdType>
std::set<size_t> Get26Neighbors3DRegular(IdType w, IdType h, IdType d, IdType id)
{
  const IdType s = w * h;

  auto nid2nidx = [w, s](IdType id, IdType &i, IdType &j, IdType &k) {
    k = id / s;
    j = (id - k*s) / w;
    i = id - k*s - j*w;
  };

  auto nidx2nid = [w, s](IdType i, IdType j, IdType k) {
    return i + j*w + k*s;
  };
 
  IdType i, j, k;
  nid2nidx(id, i, j, k);
  IdType i0 = (i > 0)   ? (i - 1) : 0, 
         i1 = (i < w-1) ? (i + 1) : w - 1, 
         j0 = (j > 0)   ? (j - 1) : 0,
         j1 = (j < h-1) ? (j + 1) : h - 1, 
         k0 = (k > 0)   ? (k - 1) : 0,
         k1 = (k < d-1) ? (k + 1) : d - 1;
  
  std::set<IdType> nb;
  nb.insert(nidx2nid(i0, j0, k0));
  nb.insert(nidx2nid(i,  j0, k0));
  nb.insert(nidx2nid(i1, j0, k0));
  
  nb.insert(nidx2nid(i0, j,  k0));
  nb.insert(nidx2nid(i,  j,  k0));
  nb.insert(nidx2nid(i1, j,  k0));
  
  nb.insert(nidx2nid(i0, j1, k0));
  nb.insert(nidx2nid(i,  j1, k0));
  nb.insert(nidx2nid(i1, j1, k0));
  
  nb.insert(nidx2nid(i0, j0, k ));
  nb.insert(nidx2nid(i,  j0, k ));
  nb.insert(nidx2nid(i1, j0, k ));
  
  nb.insert(nidx2nid(i0, j,  k ));
  // nb.insert(nidx2nid(i,  j,  k ));
  nb.insert(nidx2nid(i1, j,  k ));
  
  nb.insert(nidx2nid(i0, j1, k ));
  nb.insert(nidx2nid(i,  j1, k ));
  nb.insert(nidx2nid(i1, j1, k ));
  
  nb.insert(nidx2nid(i0, j0, k1));
  nb.insert(nidx2nid(i,  j0, k1));
  nb.insert(nidx2nid(i1, j0, k1));
  
  nb.insert(nidx2nid(i0, j,  k1));
  nb.insert(nidx2nid(i,  j,  k1));
  nb.insert(nidx2nid(i1, j,  k1));
  
  nb.insert(nidx2nid(i0, j1, k1));
  nb.insert(nidx2nid(i,  j1, k1));
  nb.insert(nidx2nid(i1, j1, k1));
  
  return nb;
}

}

#endif
