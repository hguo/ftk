#include <ftk/algorithms/uf.h>
#include <ftk/mesh/mesh.h>
#include <ftk/algorithms/uf.h>
#include <cmath>

namespace ftk {
struct Tree {
  void addNode(int i) {
    fprintf(stderr, "adding node %d\n", i);
  }
  void addArc(int i, int j) {
    fprintf(stderr, "adding arc %d<->%d\n", i, j);
  }
};
}

void joinTree(int nn, 
    const std::vector<size_t> &order,
    const std::vector<size_t> &inverseOrder,
    const std::function<double(size_t)> &value,
    const std::function<std::set<size_t>(size_t)> &neighbors)
{
  ftk::QuickUnion<size_t> uf(nn);
  ftk::Tree jt;

  for (int i=0; i<nn; i++) {
    // fprintf(stderr, "sweeping %d(%d)\n", i, order[i]); 
    // jt.addNode(i);
    fprintf(stderr, "====\nadding node %d(%d)\n", i, order[i]);

    for (auto j : neighbors(order[i])) {
      j = inverseOrder[j];
      // if (value(j) < value(i)) {
      if (j < i) {
        size_t ri = uf.root(i), rj = uf.root(j);
        fprintf(stderr, "--node %d(%d), root %d: neighbor %d(%d), root %d\n", i, order[i], ri, j, order[j], rj);
        if (ri != rj) {
          // jt.addArc(ri, rj);
          fprintf(stderr, "adding arc %d(%d)<-->%d(%d)\n", 
              ri, order[ri], rj, order[rj]);
          uf.unite(rj, ri);
        }
      }
    }
  }
}

double sinc(const double x)
{
  if (x==0) return 1;
  else return sin(x)/x;
}

int main(int argc, char **argv)
{
  const int W = 3, H = 3;
  std::vector<double> values(W*H);

#if 0
  for (int i=0; i<W*H; i++) {
    fprintf(stderr, "node %d\n", i);
    for (auto j : ftk::Get4Neighbors2DRegular<int>(W, H, i)) {
      fprintf(stderr, "--%d\n", j);
    }
  }
#endif

  for (int i=0; i<W; i++)
    for (int j=0; j<H; j++) {
      double x = i * 0.2, y = j * 0.2;
      // values[j*W+i] = sinc(sqrt(x*x + y*y));
      values[j*W+i] = 1. / sqrt(x*x + y*y + 1);
      // values[j*W+i] = (double)rand()/RAND_MAX; 
    }

  std::vector<size_t> order(W*H), inverseOrder(W*H);
  for (size_t i=0; i<W*H; i++) 
    order[i] = i;

  std::stable_sort(order.begin(), order.end(), 
      [&values](size_t p, size_t q) {
        return values[p] < values[q];
      });

  for (size_t i=0; i<W*H; i++)
    inverseOrder[order[i]] = i;

  for (size_t i=0; i<W*H; i++) {
    fprintf(stderr, "val=%f, order=%d, inverseOrder=%d\n", 
        values[i], order[i], inverseOrder[i]);
    // fprintf(stderr, "%f\n", values[order[i]]);
  }

  joinTree(W*H, order, inverseOrder,
      [&values](size_t i) {return values[i];}, 
      std::bind(ftk::Get4Neighbors2DRegular<size_t>, W, H, std::placeholders::_1));
}
