#include <ftk/mesh/mesh.h>
#include <ftk/mesh/access.h>
#include <ftk/basic/union_find.h>
#include <functional>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <map>

template <class IdType>
struct Tree {
  void addNode(IdType i) { nodes.insert(i); }

  bool hasNode(IdType i) const { return nodes.find(i) != nodes.end(); }

  void addArc(IdType lo, IdType hi) {
    upperLinks[lo].insert(hi);
    lowerLinks[hi].insert(lo);
  }
 
  bool isLeaf(IdType i) const {
    return upperDegree(i) + lowerDegree(i) == 1;
  }

  size_t upperDegree(IdType lo) const {
    auto it = upperLinks.find(lo);
    if (it != upperLinks.end()) return it->second.size();
    else return 0;
  }

  size_t lowerDegree(IdType hi) const {
    auto it = lowerLinks.find(hi);
    if (it != lowerLinks.end()) return it->second.size();
    else return 0;
  }

  IdType upperNode(IdType lo) const { // return the first upper node
    auto it = upperLinks.find(lo);
    if (it != upperLinks.end()) return *it->second.begin();
    else return IdType(-1);
  }

  IdType lowerNode(IdType hi) const { // return the first lower node
    auto it = lowerLinks.find(hi);
    if (it != lowerLinks.end()) return *it->second.begin();
    else return IdType(-1);
  }

  bool reduceNode(IdType i) { // remove a node whose up-degree and down-degree are both less or equal to 1
    size_t ud = upperDegree(i), ld = lowerDegree(i);
    // fprintf(stderr, "REDUCING NODE: ud=%zu, ld=%zu\n", ud, ld);
    if (ud > 1 || ld > 1) return false; // cannot do the reduction
    else if (ud == 1 && ld == 0) {
      lowerLinks[upperNode(i)].erase(i);
    } else if (ud == 0 && ld == 1) {
      upperLinks[lowerNode(i)].erase(i);
    } else { // ud == 1 && ld == 1
      auto hi = upperNode(i), lo = lowerNode(i);
      lowerLinks[hi].erase(i);
      upperLinks[lo].erase(i);
      lowerLinks[hi].insert(lo);
      upperLinks[lo].insert(hi);
    }
    nodes.erase(i);
    upperLinks.erase(i);
    lowerLinks.erase(i);
    return true;
  }

  std::set<IdType> getUpperLinks(IdType i) const {
    auto it = upperLinks.find(i);
    if (it != upperLinks.end()) return it.second;
    else return std::set<IdType>();
  }

  std::set<IdType> getLowerLinks(IdType i) const {
    auto it = lowerLinks.find(i);
    if (it != lowerLinks.end()) return it.second;
    else return std::set<IdType>();
  }

  void printNeighbor(IdType n, const std::function<double(size_t)> &value) const {
    std::stringstream ss;
    ss << value(n) << ": ";
    auto itu = upperLinks.find(n);
    if (itu != upperLinks.end()) {
      for (auto j : itu->second) {
        ss << value(j) << ", ";
      }
    }
    ss << ";";
    auto itl = lowerLinks.find(n);
    if (itl != lowerLinks.end()) {
      for (auto j : itl->second) {
        ss << value(j) << ", ";
      }
    }
    fprintf(stderr, "%s | ud=%zu, ld=%zu\n", ss.str().c_str(), upperDegree(n), lowerDegree(n));
  }

  void print(const std::function<double(size_t)> &value) const {
    for (auto n : nodes) {
      printNeighbor(n, value);
    }
  }

protected:
  std::set<IdType> nodes;
  std::map<IdType, std::set<IdType> > upperLinks, lowerLinks;
};

template <class IdType>
Tree<IdType> mergeTree(IdType nn, Tree<IdType>& jt, Tree<IdType>& st, const std::function<double(size_t)> &value) 
{
  Tree<IdType> ct;
  std::queue<IdType> Q;

  for (IdType i = 0; i < nn; i ++) {
    ct.addNode(i);
    if (jt.upperDegree(i) + st.lowerDegree(i) == 1) {
      // fprintf(stderr, "enqueuing %f, ud=%zu, ld=%zu\n", value(i), jt.upperDegree(i), st.lowerDegree(i));
      Q.push(i);
    }
  }

  while (!Q.empty()) {
    auto i = Q.front();
    Q.pop();

    // fprintf(stderr, "working on %f, ud=%zu, ld=%zu\n", value(i), jt.upperDegree(i), st.lowerDegree(i));
    size_t j;
    if (jt.upperDegree(i) == 0 && st.lowerDegree(i) == 1) { // upper leaf
      j = jt.lowerNode(i);
      ct.addArc(j, i);
      fprintf(stderr, "CTARC: %f <--> %f\n", value(j), value(i));
    } else if (st.lowerDegree(i) == 0 && jt.upperDegree(i) == 1) { // lower leaf
      j = st.upperNode(i);
      ct.addArc(i, j);
      fprintf(stderr, "CTARC: %f <--> %f\n", value(i), value(j));
    } else {
      // fprintf(stderr, "ERROR: %f, ud=%zu, ld=%zu\n", value(i), jt.upperDegree(i), st.lowerDegree(i));
      continue; // no longer a leaf?
      // assert(false);
    }
      
    jt.reduceNode(i);
    st.reduceNode(i);

    // fprintf(stderr, "incident j(%f): ud=%zu, ld=%zu\n", value(j), jt.upperDegree(j), st.lowerDegree(j));
    // jt.printNeighbor(j, value);
    // st.printNeighbor(j, value);
    if (jt.upperDegree(j) + st.lowerDegree(j) == 1) {
      // fprintf(stderr, "enqueuing %f\n", value(j));
      Q.push(j);
    }
  }

  return ct;
}

void buildTree(int nn,
    const std::vector<size_t> &order,
    const std::vector<size_t> &inverseOrder,
    const std::function<double(size_t)> &value,
    const std::function<std::set<size_t>(size_t)> &neighbors)
{
  ftk::quick_union<size_t> uf(nn);
  Tree<size_t> st, jt, ct;

  // split tree
  for (int i=0; i<nn; i++) {
    st.addNode(order[i]);
    // fprintf(stderr, "adding (%f)\n", value(order[i]));

    for (auto j : neighbors(order[i])) {
      j = inverseOrder[j];
      if (j < i) {
        size_t ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          st.addArc(order[rj], order[ri]);
          // fprintf(stderr, "adding arc (root=%f) <--> (root=%f)\n", 
          //     value(order[ri]), value(order[rj]));
          uf.unite(rj, ri);
        }
      }
    }
  }

  uf.reset();
  // join tree
  for (int i=nn-1; i>=0; i--) {
    jt.addNode(order[i]);
    for (auto j : neighbors(order[i])) {
      j = inverseOrder[j];
      if (j > i) {
        size_t ri = uf.root(i), rj = uf.root(j);
        if (ri != rj) {
          jt.addArc(order[ri], order[rj]);
          uf.unite(rj, ri);
        }
      }
    }
  }

  // fprintf(stderr, "join tree:\n");
  // jt.print(value);
  // fprintf(stderr, "split tree:\n");
  // st.print(value);

  // merge tree
  mergeTree<size_t>(nn, jt, st, value);
}

int main(int argc, char **argv)
{
#if 1
  const int W = 5, H = 3;
  std::vector<double> values = {
    0, 15, 4, 7, 2, 
    17, 20, 5, 10, 8, 
    13, 16, 1, 6, 3};
#endif
#if 0
  const int W = 3, H = 3;
  std::vector<double> values = {
    8, 7, 4, 
    6, 5, 2, 
    3, 1, 0};
#endif

  std::vector<size_t> order(W*H), inverseOrder(W*H);
  for (size_t i=0; i<W*H; i++) 
    order[i] = i;

  std::stable_sort(order.begin(), order.end(), 
      [&values](size_t p, size_t q) {
        return values[p] < values[q];
      });

  for (size_t i=0; i<W*H; i++)
    inverseOrder[order[i]] = i;

  buildTree(W*H, order, inverseOrder,
      [&values](size_t i) {return values[i];}, 
      std::bind(ftk::Get4Neighbors2DRegular<size_t>, W, H, std::placeholders::_1));
  
  return 0;
}
