#ifndef _FTK_REGULAR2D_H
#define _FTK_REGULAR2D_H

#include <set>

namespace ftk {

template <class IdType=size_t>
struct regular_mesh_2d {
  static void nid2nidx(IdType w, IdType h, IdType id, IdType &i, IdType &j) {
    j = id / w;
    i = id - j * w;
  };
 
  static IdType nidx2nid(IdType w, IdType h, IdType i, IdType j) {
    return i + w * j;
  };

  static std::set<IdType> get4neighbors(IdType w, IdType h, IdType id) {
    IdType i, j;
    nid2nidx(w, h, id, i, j);

    std::set<IdType> nb;
    if (i > 0) nb.insert(nidx2nid(w, h, i-1, j));
    if (i < w-1) nb.insert(nidx2nid(w, h, i+1, j));
    if (j > 0) nb.insert(nidx2nid(w, h, i, j-1));
    if (j < h-1) nb.insert(nidx2nid(w, h, i, j+1));

    return nb;
  }
  
  static std::set<IdType> get6neighbors(IdType w, IdType h, IdType id) {
    IdType i, j;
    nid2nidx(w, h, id, i, j);
    
    IdType i0 = (i > 0)   ? (i - 1) : 0, 
           i1 = (i < w-1) ? (i + 1) : w - 1, 
           j0 = (j > 0)   ? (j - 1) : 0,
           j1 = (j < h-1) ? (j + 1) : h - 1;

    std::set<IdType> nb;
    nb.insert(nidx2nid(w, h, i0, j ));
    nb.insert(nidx2nid(w, h, i0, j1));
    nb.insert(nidx2nid(w, h, i , j1));
    nb.insert(nidx2nid(w, h, i , j0));
    nb.insert(nidx2nid(w, h, i1, j ));
    nb.insert(nidx2nid(w, h, i1, j0));

    return nb;
  }
  
  static std::set<IdType> get8neighbors(IdType w, IdType h, IdType id) {
    IdType i, j;
    nid2nidx(w, h, id, i, j);
    
    IdType i0 = (i > 0)   ? (i - 1) : 0, 
           i1 = (i < w-1) ? (i + 1) : w - 1, 
           j0 = (j > 0)   ? (j - 1) : 0,
           j1 = (j < h-1) ? (j + 1) : h - 1;

    std::set<IdType> nb;
    nb.insert(nidx2nid(w, h, i0, j0));
    nb.insert(nidx2nid(w, h, i,  j0));
    nb.insert(nidx2nid(w, h, i1, j0));
    nb.insert(nidx2nid(w, h, i0, j ));
    
    nb.insert(nidx2nid(w, h, i1, j ));
    nb.insert(nidx2nid(w, h, i0, j1));
    nb.insert(nidx2nid(w, h, i,  j1));
    nb.insert(nidx2nid(w, h, i1, j1));

    return nb;
  }
};

}

#endif
