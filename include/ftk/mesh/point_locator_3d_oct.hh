#ifndef _FTK_POINT_LOCATOR_3D_OCT_HH
#define _FTK_POINT_LOCATOR_3D_OCT_HH

#include <ftk/mesh/point_locator_3d.hh>
#include <ftk/mesh/bvh3d.hh>
#include <stack>

namespace ftk {

template <typename I=int, typename F=double>
struct point_locator_3d_oct : public point_locator_3d<I, F> {
  point_locator_3d_oct(const simplicial_unstructured_3d_mesh<I, F> &m) 
    : point_locator_3d<I, F>(m) { initialize(); }
  virtual ~point_locator_3d_oct();

  void initialize();
  // I locate(const F x[], F mu[]) const { return locate_point_recursive(x, root, mu); }
  I locate(const F x[], F mu[]) const { return locate_point_nonrecursive(x, mu); }

  std::vector<bvh3d_node_t<I, F>> to_bvh() const;

protected:
  struct oct_node {
    oct_node *parent = NULL;
    oct_node *children[8] = {NULL};
    AABB<3, F> aabb;
    std::vector<AABB<3, F>> elements; // valid only for leaf nodes

    ~oct_node();
    bool is_leaf() const { return elements.size() > 0; }
    void update_bounds();
    void subdivide();
    void print() const;
  } *root = NULL;

  // void subdivide_oct_node(oct_node*);

  static bool inside_tet(const F p[], const F p1[], const F p2[], const F p3[], const F p4[], F mu[]);
  I locate_point_recursive(const F x[], const oct_node *o, F mu[]) const;
  I locate_point_nonrecursive(const F x[], F mu[]) const;
};

/////  
template <typename I, typename F>
point_locator_3d_oct<I, F>::~point_locator_3d_oct()
{
  delete root;
}

template <typename I, typename F>
I point_locator_3d_oct<I, F>::locate_point_nonrecursive(const F x[], F mu[]) const
{
  // typedef std::chrono::high_resolution_clock clock;
  // auto t0 = clock::now();
 
  const ndarray<F> &coords = this->m3->get_coords();
  const auto &conn = this->m3->get_tets();
  
  // auto t1 = clock::now();

  std::stack<oct_node*> S;
  S.push(root);

  while (!S.empty()) {
    oct_node *o = S.top();
    S.pop();

    // fprintf(stderr, "checking %p\n", o);
    
    if (o->is_leaf()) {
      const int id = o->elements[0].id;
      const auto c = conn[id];
      const int i0 = std::get<0>(c), i1 = std::get<1>(c), i2 = std::get<2>(c), i3 = std::get<3>(c);
      int succ = inside_tet(x, &coords[i0*3], &coords[i1*3], &coords[i3*3], &coords[i4*3], mu);
      if (succ) {
        // auto t2 = clock::now();
        // float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
        // float tt1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        // fprintf(stderr, "tt0=%f, tt1=%f\n", tt0, tt1);
        return id;
      }
    } else if (o->aabb[x].IsDefined()) {
      for (int j=0; j<8; j++)
        if (o->children[j] != NULL)
          S.push(o->children[j]);
    }
  }
  
  return -1;
}

template <typename I, typename F>
I point_locator_3d_oct<I, F>::locate_point_recursive(const F x[], const oct_node *o, F mu[]) const
{
  // fprintf(stderr, "locating %f, %f in node %p\n", x[0], x[1], o);

  auto m3 = this->m3;
  const auto &coords = m3->get_coords();
  const auto &conn = m3->get_tets();
  
  if (o->aabb[x].IsDefined()) {
    if (o->is_leaf()) {
      const int id = o->elements[0].id;
      const auto c = conn[id];
      const int i0 = std::get<0>(c), i1 = std::get<1>(c), i2 = std::get<2>(c), i3 = std::get<3>(c);
      int succ = inside_tet(x, &coords[i0*3], &coords[i1*3], &coords[i2*3], &coords[i3*3], mu);
      if (succ) {
        // fprintf(stderr, "leaf node %d contains! tet check=%d\n", id, result);
        return id;
      }
    } else {
      for (int j=0; j<8; j++) {
        if (o->children[j] != NULL) {
          int result = locate_point_recursive(x, o->children[j], mu);
          if (result >= 0) return result;
        }
      }
    }
  }
  return -1;
}

template <typename I, typename F>
point_locator_3d_oct<I, F>::oct_node::~oct_node() 
{
  for (int j = 0; j < 8; j ++) 
    if (children[j])
      delete children[j];
}

template <typename I, typename F>
void point_locator_3d_oct<I, F>::oct_node::print() const
{
  fprintf(stderr, "parent=%p, A={%f, %f, %f}, B={%f, %f, %f}, \
      centroid={%f, %f, %f}, children={%p, %p, %p, %p, %p, %p, %p, %p}, element=%d\n", 
      parent,
      aabb.A[0], aabb.A[1], aabb.A[2],
      aabb.B[0], aabb.B[1], aabb.B[2], 
      aabb.C[0], aabb.C[1], aabb.C[2],
      children[0], children[1], children[2], children[3], 
      children[4], children[5], children[6], children[7], 
      elements.empty() ? -1 : elements[0].id);
}

template <typename I, typename F>
void point_locator_3d_oct<I, F>::oct_node::update_bounds() 
{
  for (int i=0; i<elements.size(); i++) {
    auto &aabb1 = elements[i];
    aabb.A[0] = std::min(aabb.A[0], aabb1.A[0]);
    aabb.A[1] = std::min(aabb.A[1], aabb1.A[1]);
    aabb.A[2] = std::min(aabb.A[2], aabb1.A[2]);
    aabb.B[0] = std::max(aabb.B[0], aabb1.B[0]);
    aabb.B[1] = std::max(aabb.B[1], aabb1.B[1]);
    aabb.B[2] = std::max(aabb.B[2], aabb1.B[2]);
  }
  aabb.update_centroid();
}

template <typename I, typename F>
void point_locator_3d_oct<I, F>::oct_node::subdivide() 
{
  if (elements.size() <= 1) {
    update_bounds();
    return;
  }

  // fprintf(stderr, "subdividing %p, parent=%p, #elements=%zu\n", o, o->parent, o->elements.size());
  for (int j=0; j<8; j++) {
    children[j] = new oct_node;
    children[j]->parent = this;
  }

  // left-bottom-B
  children[0]->aabb.A[0] = aabb.A[0];
  children[0]->aabb.A[1] = aabb.A[1];
  children[0]->aabb.A[2] = aabb.A[2];
  children[0]->aabb.B[0] = aabb.C[0];
  children[0]->aabb.B[1] = aabb.C[1];
  children[0]->aabb.B[2] = aabb.C[2];
  children[0]->aabb.update_centroid();
 
  // right-bottom-B
  children[1]->aabb.A[0] = aabb.C[0];
  children[1]->aabb.A[1] = aabb.A[1];
  children[1]->aabb.A[2] = aabb.A[2];
  children[1]->aabb.B[0] = aabb.B[0];
  children[1]->aabb.B[1] = aabb.C[1];
  children[1]->aabb.B[2] = aabb.C[1];
  children[1]->aabb.update_centroid();
  
  // right-top-B
  children[2]->aabb.A[0] = aabb.C[0];
  children[2]->aabb.A[1] = aabb.C[1];
  children[2]->aabb.A[2] = aabb.A[2];
  children[2]->aabb.B[0] = aabb.B[0];
  children[2]->aabb.B[1] = aabb.B[1];
  children[2]->aabb.B[2] = aabb.C[1];
  children[2]->aabb.update_centroid();
  
  // left-top-B
  children[3]->aabb.A[0] = aabb.A[0];
  children[3]->aabb.A[1] = aabb.C[1];
  children[3]->aabb.A[2] = aabb.A[2];
  children[3]->aabb.B[0] = aabb.C[0];
  children[3]->aabb.B[1] = aabb.B[1];
  children[3]->aabb.B[2] = aabb.C[1];
  children[3]->aabb.update_centroid();
  
  // left-bottom-U
  children[0]->aabb.A[0] = aabb.A[0];
  children[0]->aabb.A[1] = aabb.A[1];
  children[0]->aabb.A[2] = aabb.C[2];
  children[0]->aabb.B[0] = aabb.C[0];
  children[0]->aabb.B[1] = aabb.C[1];
  children[0]->aabb.B[2] = aabb.B[2];
  children[0]->aabb.update_centroid();
 
  // right-bottom-U
  children[1]->aabb.A[0] = aabb.C[0];
  children[1]->aabb.A[1] = aabb.A[1];
  children[1]->aabb.A[2] = aabb.C[2];
  children[1]->aabb.B[0] = aabb.B[0];
  children[1]->aabb.B[1] = aabb.C[1];
  children[1]->aabb.B[2] = aabb.B[1];
  children[1]->aabb.update_centroid();
  
  // right-top-U
  children[2]->aabb.A[0] = aabb.C[0];
  children[2]->aabb.A[1] = aabb.C[1];
  children[2]->aabb.A[2] = aabb.C[2];
  children[2]->aabb.B[0] = aabb.B[0];
  children[2]->aabb.B[1] = aabb.B[1];
  children[2]->aabb.B[2] = aabb.B[1];
  children[2]->aabb.update_centroid();
  
  // left-top-U
  children[3]->aabb.A[0] = aabb.A[0];
  children[3]->aabb.A[1] = aabb.C[1];
  children[3]->aabb.A[2] = aabb.C[2];
  children[3]->aabb.B[0] = aabb.C[0];
  children[3]->aabb.B[1] = aabb.B[1];
  children[3]->aabb.B[2] = aabb.B[1];
  children[3]->aabb.update_centroid();

  for (int i = 0; i < elements.size(); i++) {
    for (int j = 0; j < 8; j ++) {
      if (children[j]->aabb[elements[i].C].IsDefined()) {
        children[j]->elements.push_back(elements[i]);
        break;
      }
    }
  }

  if (parent != NULL) 
    update_bounds();
  elements.clear();

  for (int j=0; j<8; j++) {
    if (children[j]->elements.empty()) {
      delete children[j];
      children[j] = NULL;
    } else {
      children[j]->subdivide();
    }
  }
}

template <typename I, typename F>
bool point_locator_3d_oct<I, F>::inside_tet(const F p[], const F p1[], const F p2[], const F p3[], F mu[])
{
  // TODO FIXME
  mu[0] = ((p2[1] - p3[1])*(p[0] - p3[0]) + (p3[0] - p2[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  mu[1] = ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  mu[2] = ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  mu[3] = 1.0 - mu[0] - mu[1] - mu[2]; 
  // // fprintf(stderr, "barycentric: %f, %f, %f\n", alpha, beta, gamma);
  return mu[0] >= 0 && mu[1] >= 0 && mu[2] >= 0 && mu[3] >= 0;
}

template <typename I, typename F>
void point_locator_3d_oct<I, F>::initialize()
{
  auto m2 = this->m2;
  const auto &coords = m2.get_coords();
  const auto &conn = m2.get_tets();
  root = new oct_node;

  // global bounds
  AABB<3, F> &aabb = root->aabb;
  for (int i = 0; i < m2.n(0); i ++) {
    aabb.A[0] = std::min(aabb.A[0], coords[3*i]);
    aabb.A[1] = std::min(aabb.A[1], coords[3*i+1]);
    aabb.A[2] = std::min(aabb.A[2], coords[3*i+2]);
    aabb.B[0] = std::max(aabb.B[0], coords[3*i]);
    aabb.B[1] = std::max(aabb.B[1], coords[3*i+1]);
    aabb.B[2] = std::max(aabb.B[2], coords[3*i+2]);
  }
  
  aabb.update_centroid();
  // aabb.print();

  std::vector<AABB<3, F>> tets(m2.n(2));
  for (int i=0; i<m2.n(2); i++) {
    const auto c = conn[id];
    const int i0 = std::get<0>(c), i1 = std::get<1>(c), i2 = std::get<2>(c), i3 = std::get<3>(c);
    double x0 = coords[i0*3], x1 = coords[i1*3], x2 = coords[i2*3], 
           y0 = coords[i0*3+1], y1 = coords[i1*3+1], y2 = coords[i2*3+1],
           z0 = coords[i0*3+2], z1 = coords[i1*3+2], z2 = coords[i2*3+2];

    tet[i].C[0] = (x0 + x1 + x2) / 3;
    tet[i].C[1] = (y0 + y1 + y2) / 3;
    tet[i].C[2] = (z0 + z1 + z2) / 3;

    tet[i].A[0] = min3(x0, x1, x2); 
    tet[i].A[1] = min3(y0, y1, y2);
    tet[i].A[2] = min3(z0, z1, z2);

    tet[i].B[0] = max3(x0, x1, x2); 
    tet[i].B[1] = max3(y0, y1, y2);
    tet[i].B[2] = max3(z0, z1, z2);

    tet[i].id = i;

    root->elements.push_back(tet[i]);
  }

  root->subdivide();
  // subdivide_oct_node(root);
  // traverseoct_node(root);
  
  // std::vector<bvh3d_node_t<I, F>> rd = to_bvh_nodes();
  
#if 0
  const double X[2] = {2.3, -0.4};
  int r1 = locate(X);
  fprintf(stderr, "r1=%d\n", r1);
  
  std::vector<bvh3d_node_t<I, F>> rd = to_bvh_nodes(root, conn, coords);

#if 1
  fprintf(stderr, "BVH built.\n");
  typedef std::chrono::high_resolution_clock clock;

  const double X[2] = {2.3, -0.4};
  auto t0 = clock::now();
  int r0 = locatePointBruteForce(X, root, nNodes, nTriangles, coords, conn);
  auto t1 = clock::now();
  int r1 = locatePointRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t2 = clock::now();
  int r2 = locatePointNonRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t3 = clock::now();
#if 0 
  float alpha, beta, gamma;
  int r3 = bvh3d_node_t<I, F>_locatePoint(rd, X[0], X[1], alpha, beta, gamma);
  auto t4 = clock::now();
  int r4 = bvh3d_node_t<I, F>_locatePoint_recursive(rd, rd, X[0], X[1], alpha, beta, gamma);
  auto t5 = clock::now();
  fprintf(stderr, "r0=%d, r1=%d, r2=%d, r3=%d, r4=%d\n", r0, r1, r2, r3, r4);
#endif
  fprintf(stderr, "x={%f, %f}, r0=%d, r1=%d, r2=%d\n", X[0], X[1], r0, r1, r2);

  float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
  float tt1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
  float tt2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
  // float tt3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count();
  // float tt4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count();
  // fprintf(stderr, "tt0=%f, tt1=%f, tt2=%f, tt3=%f, tt4=%f\n", tt0, tt1, tt2, tt3, tt4);
  fprintf(stderr, "tt0=%f, tt1=%f, tt2=%f\n", tt0, tt1, tt2);
#endif
  deleteBVH(root);

  return rd;
#endif
}

template <typename I, typename F>
std::vector<bvh3d_node_t<I, F>> point_locator_3d_oct<I, F>::to_bvh() const {
  // oct_node* r, const std::vector<int> &conn, const std::vector<double> &coords) {
  oct_node *r = root;
  
  std::map<oct_node*, int> node_map;
  std::map<int, oct_node*> node_reverse_map;
  std::stack<oct_node*> S;
  S.push(r);

  int oct_node_count = 0;
  int max_stack_size = 0;
  while (!S.empty()) {
    max_stack_size = std::max(max_stack_size, static_cast<int>(S.size()));

    oct_node *o = S.top();
    S.pop();

    int nodeId = oct_node_count ++;
    node_map[o] = nodeId;
    node_reverse_map[nodeId] = o;

    for (int j=0; j<8; j++) 
      if (o->children[j] != NULL)
        S.push(o->children[j]);
  }

  std::vector<bvh3d_node_t<I, F>> rd(oct_node_count);

  for (int i=0; i<oct_node_count; i++) {
    oct_node *o = node_reverse_map[i];
    bvh3d_node_t<I, F> &d = rd[i];

    // parent
    if (o->parent == NULL) d.parentId = -1; // root
    else d.parentId = node_map[o->parent];

    // children
    for (int j=0; j<8; j++)
      if (o->children[j] == NULL) d.childrenIds[j] = -1;
      else d.childrenIds[j] = node_map[o->children[j]];

    // bounds
    d.Ax = o->aabb.A[0];
    d.Ay = o->aabb.A[1];
    d.Az = o->aabb.A[2];
    d.Bx = o->aabb.B[0];
    d.By = o->aabb.B[1];
    d.Bz = o->aabb.B[2];
    // fprintf(stderr, "%f, %f, %f, %f\n", d.Ax, d.Ay, d.Bx, d.By);

    // tet
    if (o->is_leaf()) {
      const int id = o->elements[0].id;
      d.tetId = id;

      I tet[4];
      this->m2.get_tet(id, tet);

      F p0[3], p1[3], p2[3], p3[3];
      this->m2.get_coords(tet[0], p0);
      this->m2.get_coords(tet[1], p1);
      this->m2.get_coords(tet[2], p2);
      this->m2.get_coords(tet[3], p3);

      d.i0 = tet[0];
      d.i1 = tet[1];
      d.i2 = tet[2];
      d.i3 = tet[3];
      d.x0 = p0[0];
      d.y0 = p0[1];
      d.z0 = p0[2];
      d.x1 = p1[0];
      d.y1 = p1[1];
      d.z1 = p1[2];
      d.x2 = p2[0];
      d.y2 = p2[1];
      d.z2 = p2[2];
      d.x3 = p3[0];
      d.y3 = p3[1];
      d.z3 = p3[2];
    } else 
      d.tetId = -1;
  }

  fprintf(stderr, "oct_node_count=%d, max_stack_size=%d\n", oct_node_count, max_stack_size);
  return rd; 
  // return oct_node_count;
}

}

#endif
