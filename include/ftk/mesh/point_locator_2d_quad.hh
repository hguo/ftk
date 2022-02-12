#ifndef _FTK_POINT_LOCATOR_2D_QUAD_HH
#define _FTK_POINT_LOCATOR_2D_QUAD_HH

#include <ftk/mesh/point_locator_2d.hh>
#include <ftk/mesh/bvh2d.hh>
#include <stack>

namespace ftk {

template <typename I=int, typename F=double>
struct point_locator_2d_quad : public point_locator_2d<I, F> {
  point_locator_2d_quad(const simplicial_unstructured_2d_mesh<I, F> &m) 
    : point_locator_2d<I, F>(m) { initialize(); }
  virtual ~point_locator_2d_quad();

  void initialize();
  // I locate(const F x[], F mu[]) const { return locate_point_recursive(x, root, mu); }
  I locate(const F x[], F mu[]) const { return locate_point_nonrecursive(x, mu); }

  std::vector<bvh2d_node_t<I, F>> to_bvh() const;

protected:
  struct quad_node {
    quad_node *parent = NULL;
    quad_node *children[4] = {NULL};
    AABB<F> aabb;
    std::vector<AABB<F>> elements; // valid only for leaf nodes

    ~quad_node();
    bool is_leaf() const { return elements.size() > 0; }
    void update_bounds();
    void subdivide();
    void print() const;
  } *root = NULL;

  // void subdivide_quad_node(quad_node*);

  static bool inside_triangle(const F p[], const F p1[], const F p2[], const F p3[], F mu[]);
  I locate_point_recursive(const F x[], const quad_node *q, F mu[]) const;
  I locate_point_nonrecursive(const F x[], F mu[]) const;
};

/////  
template <typename I, typename F>
point_locator_2d_quad<I, F>::~point_locator_2d_quad()
{
  delete root;
}

template <typename I, typename F>
I point_locator_2d_quad<I, F>::locate_point_nonrecursive(const F x[], F mu[]) const
{
  // typedef std::chrono::high_resolution_clock clock;
  // auto t0 = clock::now();
 
  const ndarray<F> &coords = this->m2.get_coords();
  const ndarray<I> &conn = this->m2.get_triangles();
  
  // auto t1 = clock::now();

  std::stack<quad_node*> S;
  S.push(root);

  while (!S.empty()) {
    quad_node *q = S.top();
    S.pop();

    // fprintf(stderr, "checking %p\n", q);
    
    if (q->is_leaf()) {
      const int id = q->elements[0].id;
      const int i0 = conn[id*3], i1 = conn[id*3+1], i2 = conn[id*3+2];
      int succ = inside_triangle(x, &coords[i0*2], &coords[i1*2], &coords[i2*2], mu);
      if (succ) {
        // auto t2 = clock::now();
        // float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
        // float tt1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        // fprintf(stderr, "tt0=%f, tt1=%f\n", tt0, tt1);
        return id;
      }
    } else if (q->aabb.contains(x)) {
      for (int j=0; j<4; j++)
        if (q->children[j] != NULL)
          S.push(q->children[j]);
    }
  }
  
  return -1;
}

template <typename I, typename F>
I point_locator_2d_quad<I, F>::locate_point_recursive(const F x[], const quad_node *q, F mu[]) const
{
  // fprintf(stderr, "locating %f, %f in node %p\n", x[0], x[1], q);

  auto m2 = this->m2;
  const auto &coords = m2.get_coords();
  const auto &conn = m2.get_triangles();
  
  if (q->aabb.contains(x)) {
    if (q->is_leaf()) {
      const int id = q->elements[0].id;
      const int i0 = conn[id*3], i1 = conn[id*3+1], i2 = conn[id*3+2];
      int succ = inside_triangle(x, &coords[i0*2], &coords[i1*2], &coords[i2*2], mu);
      if (succ) {
        // fprintf(stderr, "leaf node %d contains! triangle check=%d\n", id, result);
        return id;
      }
    } else {
      for (int j=0; j<4; j++) {
        if (q->children[j] != NULL) {
          int result = locate_point_recursive(x, q->children[j], mu);
          if (result >= 0) return result;
        }
      }
    }
  }
  return -1;
}

template <typename I, typename F>
point_locator_2d_quad<I, F>::quad_node::~quad_node() 
{
  for (int j = 0; j < 4; j ++) 
    if (children[j])
      delete children[j];
}

template <typename I, typename F>
void point_locator_2d_quad<I, F>::quad_node::print() const
{
  fprintf(stderr, "parent=%p, A={%f, %f}, B={%f, %f}, centroid={%f, %f}, children={%p, %p, %p, %p}, element=%d\n", 
      parent, aabb.A[0], aabb.A[1], aabb.B[0], aabb.B[1], aabb.C[0], aabb.C[1], 
      children[0], children[1], children[2], children[3], 
      elements.empty() ? -1 : elements[0].id);
}

template <typename I, typename F>
void point_locator_2d_quad<I, F>::quad_node::update_bounds() 
{
  for (int i=0; i<elements.size(); i++) {
    auto &aabb1 = elements[i];
    aabb.A[0] = std::min(aabb.A[0], aabb1.A[0]);
    aabb.A[1] = std::min(aabb.A[1], aabb1.A[1]);
    aabb.B[0] = std::max(aabb.B[0], aabb1.B[0]);
    aabb.B[1] = std::max(aabb.B[1], aabb1.B[1]);
  }
  aabb.update_centroid();
}

template <typename I, typename F>
void point_locator_2d_quad<I, F>::quad_node::subdivide() 
{
  if (elements.size() <= 1) {
    update_bounds();
    return;
  }

  // fprintf(stderr, "subdividing %p, parent=%p, #elements=%zu\n", q, q->parent, q->elements.size());
  for (int j=0; j<4; j++) {
    children[j] = new quad_node;
    children[j]->parent = this;
  }

  // left-bottom
  children[0]->aabb.A[0] = aabb.A[0];
  children[0]->aabb.A[1] = aabb.A[1];
  children[0]->aabb.B[0] = aabb.C[0];
  children[0]->aabb.B[1] = aabb.C[1];
  children[0]->aabb.update_centroid();
 
  // right-bottom
  children[1]->aabb.A[0] = aabb.C[0];
  children[1]->aabb.A[1] = aabb.A[1];
  children[1]->aabb.B[0] = aabb.B[0];
  children[1]->aabb.B[1] = aabb.C[1];
  children[1]->aabb.update_centroid();
  
  // right-top
  children[2]->aabb.A[0] = aabb.C[0];
  children[2]->aabb.A[1] = aabb.C[1];
  children[2]->aabb.B[0] = aabb.B[0];
  children[2]->aabb.B[1] = aabb.B[1];
  children[2]->aabb.update_centroid();
  
  // left-top
  children[3]->aabb.A[0] = aabb.A[0];
  children[3]->aabb.A[1] = aabb.C[1];
  children[3]->aabb.B[0] = aabb.C[0];
  children[3]->aabb.B[1] = aabb.B[1];
  children[3]->aabb.update_centroid();

  for (int i = 0; i < elements.size(); i++) {
    for (int j = 0; j < 4; j ++) {
      if (children[j]->aabb.contains(elements[i].C)) {
        children[j]->elements.push_back(elements[i]);
        break;
      }
    }
  }

  if (parent != NULL) 
    update_bounds();
  elements.clear();

  for (int j=0; j<4; j++) {
    if (children[j]->elements.empty()) {
      delete children[j];
      children[j] = NULL;
    } else {
      children[j]->subdivide();
    }
  }
}

template <typename I, typename F>
bool point_locator_2d_quad<I, F>::inside_triangle(const F p[], const F p1[], const F p2[], const F p3[], F mu[])
{
  mu[0] = ((p2[1] - p3[1])*(p[0] - p3[0]) + (p3[0] - p2[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  mu[1] = ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  mu[2] = 1.0 - mu[0] - mu[1]; 
  // fprintf(stderr, "barycentric: %f, %f, %f\n", alpha, beta, gamma);
  return mu[0] >= 0 && mu[1] >= 0 && mu[2] >= 0;
}

template <typename I, typename F>
void point_locator_2d_quad<I, F>::initialize()
{
  auto m2 = this->m2;
  const auto &coords = m2.get_coords();
  const auto &conn = m2.get_triangles();
  root = new quad_node;

  // global bounds
  AABB<F> &aabb = root->aabb;
  for (int i = 0; i < m2.n(0); i ++) {
    aabb.A[0] = std::min(aabb.A[0], coords[2*i]);
    aabb.A[1] = std::min(aabb.A[1], coords[2*i+1]);
    aabb.B[0] = std::max(aabb.B[0], coords[2*i]);
    aabb.B[1] = std::max(aabb.B[1], coords[2*i+1]);
  }
  aabb.update_centroid();
  // aabb.print();

  std::vector<AABB<F>> triangles(m2.n(2));
  for (int i=0; i<m2.n(2); i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    double x0 = coords[i0*2], x1 = coords[i1*2], x2 = coords[i2*2], 
           y0 = coords[i0*2+1], y1 = coords[i1*2+1], y2 = coords[i2*2+1];

    triangles[i].C[0] = (x0 + x1 + x2) / 3;
    triangles[i].C[1] = (y0 + y1 + y2) / 3;
    triangles[i].A[0] = min3(x0, x1, x2); 
    triangles[i].A[1] = min3(y0, y1, y2); 
    triangles[i].B[0] = max3(x0, x1, x2); 
    triangles[i].B[1] = max3(y0, y1, y2);
    triangles[i].id = i;

    root->elements.push_back(triangles[i]);
  }

  root->subdivide();
  // subdivide_quad_node(root);
  // traversequad_node(root);
  
  // std::vector<bvh2d_node_t<I, F>> rd = to_bvh_nodes();
  
#if 0
  const double X[2] = {2.3, -0.4};
  int r1 = locate(X);
  fprintf(stderr, "r1=%d\n", r1);
  
  std::vector<bvh2d_node_t<I, F>> rd = to_bvh_nodes(root, conn, coords);

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
  int r3 = bvh2d_node_t<I, F>_locatePoint(rd, X[0], X[1], alpha, beta, gamma);
  auto t4 = clock::now();
  int r4 = bvh2d_node_t<I, F>_locatePoint_recursive(rd, rd, X[0], X[1], alpha, beta, gamma);
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
std::vector<bvh2d_node_t<I, F>> point_locator_2d_quad<I, F>::to_bvh() const {
  // quad_node* r, const std::vector<int> &conn, const std::vector<double> &coords) {
  quad_node *r = root;
  
  std::map<quad_node*, int> node_map;
  std::map<int, quad_node*> node_reverse_map;
  std::stack<quad_node*> S;
  S.push(r);

  int quad_node_count = 0;
  int max_stack_size = 0;
  while (!S.empty()) {
    max_stack_size = std::max(max_stack_size, static_cast<int>(S.size()));

    quad_node *q = S.top();
    S.pop();

    int nodeId = quad_node_count ++;
    node_map[q] = nodeId;
    node_reverse_map[nodeId] = q;

    for (int j=0; j<4; j++) 
      if (q->children[j] != NULL)
        S.push(q->children[j]);
  }

  std::vector<bvh2d_node_t<I, F>> rd(quad_node_count);

  for (int i=0; i<quad_node_count; i++) {
    quad_node *q = node_reverse_map[i];
    bvh2d_node_t<I, F> &d = rd[i];

    // parent
    if (q->parent == NULL) d.parentId = -1; // root
    else d.parentId = node_map[q->parent];

    // children
    for (int j=0; j<4; j++)
      if (q->children[j] == NULL) d.childrenIds[j] = -1;
      else d.childrenIds[j] = node_map[q->children[j]];

    // bounds
    d.Ax = q->aabb.A[0];
    d.Ay = q->aabb.A[1];
    d.Bx = q->aabb.B[0];
    d.By = q->aabb.B[1];
    // fprintf(stderr, "%f, %f, %f, %f\n", d.Ax, d.Ay, d.Bx, d.By);

    // triangle
    if (q->is_leaf()) {
      const int id = q->elements[0].id;
      d.triangleId = id;

      I tri[3];
      this->m2.get_triangle(id, tri);

      F p0[3], p1[3], p2[3];
      this->m2.get_coords(tri[0], p0);
      this->m2.get_coords(tri[1], p1);
      this->m2.get_coords(tri[2], p2);

      d.i0 = tri[0];
      d.i1 = tri[1];
      d.i2 = tri[2];
      d.x0 = p0[0];
      d.y0 = p0[1];
      d.x1 = p1[0];
      d.y1 = p1[1];
      d.x2 = p2[0];
      d.y2 = p2[1];
    } else 
      d.triangleId = -1;
  }

  fprintf(stderr, "quad_node_count=%d, max_stack_size=%d\n", quad_node_count, max_stack_size);
  return rd; 
  // return quad_node_count;
}


}

#endif
