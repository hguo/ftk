#ifndef _FTK_MESH
#define _FTK_MESH

#include <set>

namespace ftk {

template <class IdType>
std::set<size_t> Get4Neighbors2DRegular(IdType w, IdType h, IdType id) {
  auto nid2nidx = [w, h](IdType id, IdType &i, IdType &j) {
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
std::set<size_t> Get8Neighbors2DRegular(IdType w, IdType h, IdType id);

template <class IdType>
std::set<size_t> Get6Neighbors3DRegular(IdType W, IdType H, IdType D, IdType id);

template <class IdType>
std::set<size_t> Get18Neighbors3DRegular(IdType W, IdType H, IdType D, IdType id);

template <class IdType>
std::set<size_t> Get26Neighbors3DRegular(IdType W, IdType H, IdType D, IdType id);

}

#endif
