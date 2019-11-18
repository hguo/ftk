#ifndef _FTK_TRACK_CRITICAL_POINTS_3D_REGULAR_SERIAL_STREAMING_HH
#define _FTK_TRACK_CRITICAL_POINTS_3D_REGULAR_SERIAL_STREAMING_HH

#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <deque>

namespace ftk {

struct critical_point_tracker_3d_regular_streaming : public critical_point_tracker_3d_regular {
  critical_point_tracker_3d_regular_streaming() {}

  void push_input_scalar_field(const ndarray<double>& scalar0) {scalar.push_front(scalar0); has_scalar_field = true;}
  void push_input_vector_field(const ndarray<double>& V0) {V.push_front(V0); has_vector_field = true;}
  void push_input_jacobian_field(const ndarray<double>& gradV0) {gradV.push_front(gradV0); has_jacobian_field = true;}

  void advance_timestep();
  virtual void update();

protected:
  virtual void update_timestep();

protected:
  typedef regular_simplex_mesh_element element_t;

  std::deque<ndarray<double>> scalar, V, gradV;
  int start_timestep = 0, current_timestep = 0;

protected:
  void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[4][3]) const;
  void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[4]) const;
  void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[4][3][3]) const;
};


////////////////////
void critical_point_tracker_3d_regular_streaming::advance_timestep()
{
  update_timestep();

  // fprintf(stderr, "advancing timestep!\n");
  const int nt = 2;
  if (scalar.size() > nt) scalar.pop_back();
  if (V.size() > nt) V.pop_back();
  if (gradV.size() > nt) gradV.pop_back();

  // fprintf(stderr, "scalar.size=%lu\n", scalar.size());

  has_scalar_field = false;
  has_vector_field = false;
  has_jacobian_field = false;
  
  current_timestep ++;
}

void critical_point_tracker_3d_regular_streaming::update()
{
#if 0 // doesn't make difference
  std::vector<int> lb, ub;
  m.get_lb_ub(lb, ub);
  ub[3] = current_timestep - 1;
  m.set_lb_ub(lb, ub);
#endif
  // fprintf(stderr, "trace1\n");
  // trace_intersections();
  fprintf(stderr, "trace2\n");
  trace_connected_components();
  fprintf(stderr, "done.\n");
}

void critical_point_tracker_3d_regular_streaming::update_timestep()
{
  // initializing vector fields
  if (has_scalar_field) { 
    if (!has_vector_field) push_input_vector_field(gradient3D(scalar[0])); // 0 is the current timestep; 1 is the last timestep
    if (!has_jacobian_field) push_input_jacobian_field(jacobian3D(V[0]));
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {
          static_cast<int>(V[0].dim(1)-3), 
          static_cast<int>(V[0].dim(2)-3), 
          static_cast<int>(V[0].dim(3)-3), 
          std::numeric_limits<int>::max()});
    else
      m.set_lb_ub({0, 0, 0}, {
          static_cast<int>(V[0].dim(1)-1), 
          static_cast<int>(V[0].dim(2)-1), 
          static_cast<int>(V[0].dim(3)-1), 
          std::numeric_limits<int>::max()});
  }

  // scan 3-simplices
  // fprintf(stderr, "tracking 3D critical points...\n");
  auto func3 = [=](element_t e) {
      critical_point_3dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
        fprintf(stderr, "%f, %f, %f, %f, type=%d\n", cp[0], cp[1], cp[2], cp[3], cp.type);
      }
    };

  if (V.size() >= 2)
    m.element_for_interval(3, current_timestep-1, current_timestep, func3);

  m.element_for_ordinal(3, current_timestep, func3);
}

void critical_point_tracker_3d_regular_streaming::simplex_vectors(
    const std::vector<std::vector<int>>& vertices, double v[4][3]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    for (int j = 0; j < 3; j ++)
      v[i][j] = V[iv](j, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
}

void critical_point_tracker_3d_regular_streaming::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[4]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    values[i] = scalar[iv](vertices[i][0], vertices[i][1], vertices[i][2]);
  }
}

void critical_point_tracker_3d_regular_streaming::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[4][3][3]) const
{
  for (int i = 0; i < 4; i ++) {
    const int iv = vertices[i][3] == 0 ? 0 : 1;
    for (int j = 0; j < 3; j ++) {
      for (int k = 0; k < 3; k ++) {
        Js[i][j][k] = gradV[iv](k, j, vertices[i][0], vertices[i][1], vertices[i][2]);
      }
    }
  }
}

}

#endif
