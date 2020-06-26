#ifndef _FTK_CCL2D_HH
#define _FTK_CCL2D_HH

#include <ftk/ndarray.hh>

namespace ftk {

// https://en.wikipedia.org/wiki/Hoshen%E2%80%93Kopelman_algorithm
template <typename LabelIdType>
LabelIdType hoshen_kopelman_2d(ndarray<LabelIdType>& matrix)
{
  LabelIdType largest_label = 0;

  if (matrix.nd() != 2) return 0;
  std::vector<LabelIdType> labels(matrix.nelem() / 2);

  auto make_set = [&]() {
    labels[0] ++;
    labels[labels[0]] = labels[0];
    return labels[0];
  };

  auto find = [&](LabelIdType x) {
    LabelIdType y = x;
    while (labels[y] != y) 
      y = labels[y];
    while (labels[x] != x) {
      LabelIdType z = labels[x];
      labels[x] = y;
      x = z;
    }
    return y;
  };

  auto unite = [&](LabelIdType x, LabelIdType y) {
    return labels[find(x)] = find(y);
  };

  // first pass
  for (auto j = 0; j < matrix.dim(1); j ++) {
    for (auto i = 0; i < matrix.dim(0); i ++) {
      if (matrix(i, j)) {
        LabelIdType left = (i == 0) ? 0 : matrix(i-1, j);
        LabelIdType up = (j == 0) ? 0 : matrix(i, j-1);
        switch (!!left + !!up) {
        case 0: 
          matrix(i, j) = make_set(); 
          break;
        case 1:
          matrix(i, j) = std::max(up, left);
          break;
        case 2:
          matrix(i, j) = unite(up, left);
          break;
        }
      }
    }
  }

  // second pass
  std::vector<LabelIdType> new_labels(labels.size());
  for (auto j = 0; j < matrix.dim(1); j ++) {
    for (auto i = 0; i < matrix.dim(0); i ++) {
      if (matrix(i, j)) {
        LabelIdType x = find(matrix(i, j));
        if (new_labels[x] == 0) {
          new_labels[0] ++;
          new_labels[x] = new_labels[0];
        }
        matrix(i, j) = new_labels[x];
      }
    }
  }

  return new_labels[0];
}

}

#endif
