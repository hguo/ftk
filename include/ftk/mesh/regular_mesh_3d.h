#ifndef _FTK_REGULAR3D_H
#define _FTK_REGULAR3D_H

#include <set>

namespace ftk {

template <class IdType=size_t>
struct regular_mesh_3d {
  static void nid2nidx(IdType wh /* W*H */, IdType w, IdType id, IdType &i, IdType &j, IdType &k) {
    k = id / wh;
    j = (id - k*wh) / w;
    i = id - k*wh - j*w;
  };
 
  static IdType nidx2nid(IdType wh, IdType w, IdType i, IdType j, IdType k) {
    return i + j*w + k*wh;
  };

  static std::set<IdType> get6neighbors(IdType w, IdType h, IdType d, IdType id) {
    const IdType wh = w*h; 

    IdType i, j, k;
    nid2nidx(wh, w, id, i, j, k);
    IdType i0 = (i > 0)   ? (i - 1) : 0, 
           i1 = (i < w-1) ? (i + 1) : w - 1, 
           j0 = (j > 0)   ? (j - 1) : 0,
           j1 = (j < h-1) ? (j + 1) : h - 1, 
           k0 = (k > 0)   ? (k - 1) : 0,
           k1 = (k < d-1) ? (k + 1) : d - 1;
    
    std::set<IdType> nb;
    nb.insert(nidx2nid(wh, w, i,  j,  k0));
    nb.insert(nidx2nid(wh, w, i,  j,  k1));
    nb.insert(nidx2nid(wh, w, i0, j,  k ));
    nb.insert(nidx2nid(wh, w, i1, j,  k ));
    nb.insert(nidx2nid(wh, w, i,  j0, k ));
    nb.insert(nidx2nid(wh, w, i,  j1, k ));

    return nb;
  }
  
  static std::set<IdType> get18neighbors(IdType w, IdType h, IdType d, IdType id);
  
  static std::set<IdType> get26neighbors(IdType w, IdType h, IdType d, IdType id) {
    const IdType wh = w*h; 

    IdType i, j, k;
    nid2nidx(wh, w, id, i, j, k);
    IdType i0 = (i > 0)   ? (i - 1) : 0, 
           i1 = (i < w-1) ? (i + 1) : w - 1, 
           j0 = (j > 0)   ? (j - 1) : 0,
           j1 = (j < h-1) ? (j + 1) : h - 1, 
           k0 = (k > 0)   ? (k - 1) : 0,
           k1 = (k < d-1) ? (k + 1) : d - 1;
    
    std::set<IdType> nb;
    nb.insert(nidx2nid(wh, w, i0, j0, k0));
    nb.insert(nidx2nid(wh, w, i,  j0, k0));
    nb.insert(nidx2nid(wh, w, i1, j0, k0));
    
    nb.insert(nidx2nid(wh, w, i0, j,  k0));
    nb.insert(nidx2nid(wh, w, i,  j,  k0));
    nb.insert(nidx2nid(wh, w, i1, j,  k0));
    
    nb.insert(nidx2nid(wh, w, i0, j1, k0));
    nb.insert(nidx2nid(wh, w, i,  j1, k0));
    nb.insert(nidx2nid(wh, w, i1, j1, k0));
    
    nb.insert(nidx2nid(wh, w, i0, j0, k ));
    nb.insert(nidx2nid(wh, w, i,  j0, k ));
    nb.insert(nidx2nid(wh, w, i1, j0, k ));
    
    nb.insert(nidx2nid(wh, w, i0, j,  k ));
    // nb.insert(nidx2nid(i,  j,  k ));
    nb.insert(nidx2nid(wh, w, i1, j,  k ));
    
    nb.insert(nidx2nid(wh, w, i0, j1, k ));
    nb.insert(nidx2nid(wh, w, i,  j1, k ));
    nb.insert(nidx2nid(wh, w, i1, j1, k ));
    
    nb.insert(nidx2nid(wh, w, i0, j0, k1));
    nb.insert(nidx2nid(wh, w, i,  j0, k1));
    nb.insert(nidx2nid(wh, w, i1, j0, k1));
    
    nb.insert(nidx2nid(wh, w, i0, j,  k1));
    nb.insert(nidx2nid(wh, w, i,  j,  k1));
    nb.insert(nidx2nid(wh, w, i1, j,  k1));
    
    nb.insert(nidx2nid(wh, w, i0, j1, k1));
    nb.insert(nidx2nid(wh, w, i,  j1, k1));
    nb.insert(nidx2nid(wh, w, i1, j1, k1));

    return nb;
  }
};

}

#endif
