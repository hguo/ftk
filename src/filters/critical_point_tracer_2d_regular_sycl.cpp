#include <CL/sycl.hpp>
#include <vector>
#include <ftk/config.hh>
#include <ftk/mesh/lattice.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/critical_point_test.hh>
#include "common.cuh"

static inline bool check_simplex_cp2t_ordinal()
{
  return true;
}

static std::vector<cp_t> extract_cp2dt(
    int scope,
    int current_timestep,
    const lattice3_t& domain,
    const lattice3_t& core, 
    const lattice2_t& ext, 
    const double *V, // 4D array: 2*W*H*2
    bool use_explicit_coords,
    const double *coords)
{
  const unsigned int maxcps = 32768;

  std::vector<cp_t> results;
  cp_t *buf = (cp_t*)malloc(sizeof(cp_t) * maxcps); // TODO: there is a limit of #cps
  memset(buf, 0, sizeof(cp_t) * maxcps);

  unsigned int counter = 0;
  {
    cl::sycl::buffer<double, 1> dV(V, 2 * 2 * ext.n());
    auto dscope = cl::sycl::buffer<int>( &scope, 1 );
    auto dcounter = cl::sycl::buffer<unsigned int>( &counter, 1 );
    auto dbuf = cl::sycl::buffer<cp_t>( buf, maxcps );

    sycl::gpu_selector selector;
    // sycl::host_selector selector;
    cl::sycl::queue q(selector);

    std::cerr << "Running on "
              << q.get_device().get_info<sycl::info::device::name>()
              << std::endl;

    { // if (scope == scope_ordinal) {
      const size_t ntasks = core.n() * (scope == scope_ordinal ? ntypes_3_2<scope_ordinal>() : ntypes_3_2<scope_interval>());
      cl::sycl::range<1> work_items{ntasks};

      fprintf(stderr, "ntasks=%zu\n", ntasks);

      q.submit([&](cl::sycl::handler &cgh) {
        // auto V = dV.get_access<cl::sycl::access::mode::read>(cgh);
        // sycl::accessor V{dV, cgh, sycl::read_only, sycl::noinit};
        // auto scope = dscope.get_access<cl::sycl::access::mode::read>(cgh);
        auto V = dV.get_access<cl::sycl::access::mode::read>(cgh);
        auto counter = dcounter.get_access<cl::sycl::access::mode::write>(cgh);
        auto buf = dbuf.get_access<cl::sycl::access::mode::write>(cgh);
        sycl::stream out(1024, 256, cgh); // debug stream

        cgh.parallel_for<class cp2dtk>(work_items, [=](cl::sycl::id<1> tid) {
          const element32_t e = scope == 1 ? 
            element32_from_index<scope_ordinal>(core, tid) : 
            element32_from_index<scope_interval>(core, tid);

          int vertices[3][3], indices[3];
          size_t local_indices[3];
          for (int i = 0; i < 3; i ++) {
            for (int j = 0; j < 3; j ++) {
              vertices[i][j] = e.corner[j] 
                + (scope == 1 ? // scope_ordinal ? 
                  unit_simplex_offset_3_2<scope_ordinal>(e.type, i, j) :
                  unit_simplex_offset_3_2<scope_interval>(e.type, i, j));

              if (vertices[i][j] < domain.st[j] || 
                   vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
              // if (vertices[i][j] < core[0].st[j] || 
              //     vertices[i][j] > core[0].st[j] + core[0].sz[j] - 1)
                return; // false;
            }
            indices[i] = domain.to_index(vertices[i]);
            local_indices[i] = ext.to_index(vertices[i]);
          }
         
          printf("local_indices=%d, %d, %d\n", local_indices[0], local_indices[1], local_indices[2]);
          // out << "local_indices=" << local_indices[0] << "," << local_indices[1] << "," << local_indices[2] << cl::sycl::endl;
          // out << "vf[0]=" << vf[0][0] << "," << vf[0][1] << cl::sycl::endl;
          // return;

          double v[3][2];
          long long vf[3][2] = {0};
          for (int i = 0; i < 3; i ++) {
            // size_t k = ext.to_index(vertices[i]);
            const size_t k = local_indices[i];
            for (int j = 0; j < 2; j ++) {
              const int offset = (scope == 1) ? 
                unit_simplex_offset_3_2<scope_ordinal>(e.type, i, 2/*time dimension id*/) : 
                unit_simplex_offset_3_2<scope_interval>(e.type, i, 2/*time dimension id*/);

              v[i][j] = V[offset*2*ext.n() + k*2+j]; // V[unit_simplex_offset_3_2<scope_ordinal>(e.type, i, 2/*time dimension id*/)][k*2+j];
              // out << k*2+j << "," << v[i][j] << cl::sycl::endl;
              vf[i][j] = v[i][j] * 33554432; // 32768L; // WIP
            }
          }

          bool succ = ftk::robust_critical_point_in_simplex2(vf, indices);

          if (succ) {
            cl::sycl::atomic<unsigned> atomic_counter(cl::sycl::global_ptr<unsigned int>{&counter[0]});
            auto k = atomic_counter.fetch_add(1);
            cp_t &cp = buf[ k ];

            double mu[3], cond;
            bool succ2 = ftk::inverse_lerp_s2v2(v, mu, &cond);
            if (!succ2) ftk::clamp_barycentric<3>(mu);

            double X[3][3], x[3];
            for (int i = 0; i < 3; i ++)
              for (int j = 0; j < 3; j ++)
                X[i][j] = vertices[i][j];
            ftk::lerp_s2v3(X, mu, x);
            cp.x[0] = x[0];
            cp.x[1] = x[1];
            cp.t = x[2];
            cp.cond = cond;
            cp.tag = tid;
          }
        });
      });
      q.wait();
    } 
  }

  results.resize(counter);
  memcpy((void*)results.data(), buf, counter * sizeof(cp_t));

  free(buf);
  fprintf(stderr, "#results=%zu, counter=%d\n", results.size(), counter);

  return results;
}

std::vector<cp_t> extract_cp2dt_sycl(
    int scope, 
    int current_timestep,
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *Vc, 
    const double *Vn, 
    const double *Jc, 
    const double *Jn, 
    const double *Sc,
    const double *Sn, 
    bool use_explicit_coords,
    const double *coords)
{
  lattice3_t D(domain);
  lattice3_t C(core);
  lattice2_t E(ext);

  double *V = (double*)malloc(2*2*sizeof(double)*ext.n());
  if (Vc) memcpy(V, Vc, 2*sizeof(double)*ext.n());
  if (Vn) memcpy(V+2*ext.n(), Vn, 2*sizeof(double)*ext.n());

  // std::cerr << "domain=" << domain 
  //   << ", core=" << core << ", current_timestep=" 
  //   << current_timestep << std::endl;

  auto results = extract_cp2dt(scope, current_timestep, 
      D, C, E, 
      V, // Vc, Vn, Jc, Jn, Sc, Sn, 
      use_explicit_coords, coords);

  // free(V);

  return results;
}
