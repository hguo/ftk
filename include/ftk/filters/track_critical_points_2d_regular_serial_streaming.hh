#ifndef _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_STREAMING_HH
#define _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_STREAMING_HH

#include <ftk/filters/track_critical_points_2d_regular_serial.hh>
#include <deque>

namespace ftk {

struct track_critical_points_2d_regular_serial_streaming : public track_critical_points_2d_regular_serial {
  track_critical_points_2d_regular_serial_streaming() {}

  void push_input_scalar_field(const ndarray<double>& scalar0) {scalar.push_front(scalar0); scalar_updated = true;}
  void push_input_vector_field(const ndarray<double>& V0) {V.push_front(V0); vector_updated = true;}
  void push_input_jacobian_field(const ndarray<double>& gradV0) {gradV.push_front(gradV0); jacobian_updated = true;}

  void update();
  void advance_timestep();

protected:
  typedef regular_simplex_mesh_element element_t;

  std::deque<ndarray<double>> scalar, V, gradV;
  int start_timestep = 0, current_timestep = 0;

  bool scalar_updated = false, vector_updated = false, jacobian_updated = false;

protected:
  bool check_simplex_streaming(const element_t& s, critical_point_2dt_t& cp);
};


////////////////////
void track_critical_points_2d_regular_serial_streaming::advance_timestep()
{
  fprintf(stderr, "advancing timestep!\n");

  const int nt = 2;

  if (scalar.size() >= nt) scalar.pop_front();
  if (V.size() >= nt) V.pop_front();
  if (gradV.size() >= nt) gradV.pop_front();

  scalar_updated = false;
  vector_updated = false;
  jacobian_updated = false;
  
  current_timestep ++;
}

void track_critical_points_2d_regular_serial_streaming::update()
{
  // initializing vector fields
  if (scalar_updated) { 
    if (!vector_updated) push_input_vector_field(gradient2D(scalar[0])); // 0 is the current timestep; 1 is the last timestep
    if (!jacobian_updated) push_input_jacobian_field(jacobian2D(V[0]));
    has_jacobian = true;
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {
          static_cast<int>(V[0].dim(1)-3), 
          static_cast<int>(V[0].dim(2)-3), 
          std::numeric_limits<int>::max()});
    else
      m.set_lb_ub({0, 0, 0}, {
          static_cast<int>(V[0].dim(1)-1), 
          static_cast<int>(V[0].dim(2)-1), 
          std::numeric_limits<int>::max()});
  }

  // scan 2-simplices
  fprintf(stderr, "tracking 2D critical points...\n");
  auto func2 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (check_simplex_streaming(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
      }
    };

  if (V.size() >= 2)
    m.element_for_interval(2, current_timestep-1, current_timestep, func2);

  m.element_for_ordinal(2, current_timestep, func2);

#if 0
  // scan 3-simplices to get connected components
  fprintf(stderr, "trace intersections...\n");
  union_find<element_t> uf;
  for (const auto &kv : discrete_critical_points) 
    uf.add(kv.first);

  m.element_for(3, [&](const ftk::regular_simplex_mesh_element& f) {
    const auto sides = f.sides();
    std::set<element_t> intersected_sides;

    for (const auto& side : sides)
      if (uf.has(side)) 
        intersected_sides.insert(side);

    if (intersected_sides.size() > 1) {// the size of intersected_size should be 0 or 2
      for (auto it = std::next(intersected_sides.begin(), 1); it != intersected_sides.end(); it ++) {
        std::lock_guard<std::mutex> guard(mutex);
        uf.unite(*intersected_sides.begin(), *it);
      }
    }
  });
  uf.get_sets(connected_components);

  // convert connected components to traced critical points
  fprintf(stderr, "tracing critical points...\n");
  trace_connected_components();
#endif
}

bool track_critical_points_2d_regular_serial_streaming::check_simplex_streaming(
    const element_t& e,
    critical_point_2dt_t& cp)
{
  if (!e.valid()) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(); // obtain the vertices of the simplex
  
  double v[3][2]; // obtain vector values
  for (int i = 0; i < 3; i ++) {
    const int iv = vertices[i][2] == 0 ? 0 : 1;
    for (int j = 0; j < 2; j ++)
      v[i][j] = V[iv](j, vertices[i][0], vertices[i][1]);
  }
 
  double mu[3]; // check intersection
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // lerp position
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
  lerp_s2v3(X, mu, cp.x);

  fprintf(stderr, "%f, %f, %f\n", cp.x[0], cp.x[1], cp.x[2]);

  if (!scalar.empty()) {
    double values[3];
    for (int i = 0; i < 3; i ++) {
      const int iv = vertices[i][2] == 0 ? 0 : 1;
      values[i] = scalar[iv](vertices[i][0], vertices[i][1]);
    }
    cp.scalar = lerp_s2(values, mu);
  }

  if (has_jacobian) {
    // derive jacobian
    double J[2][2] = {0};
    if (gradV.empty()) {
      // TODO: jacobian is not given
    } else { // lerp jacobian
      double Js[3][2][2];
      for (int i = 0; i < 3; i ++) {
        const int iv = vertices[i][2] == 0 ? 0 : 1;
        for (int j = 0; j < 2; j ++) {
          for (int k = 0; k < 2; k ++) {
            Js[i][j][k] = gradV[iv](k, j, vertices[i][0], vertices[i][1]);
          }
        }
      }
      lerp_s2m2x2(Js, mu, J);
    }
    cp.type = critical_point_type_2d(J, symmetric_jacobian);
    if (cp.type & type_filter) return true;
    else return false;
  } else return true;
}

}

#endif
