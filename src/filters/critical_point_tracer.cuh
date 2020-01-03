
struct cpt_cfg_t {
  double threshold;
  unsigned int type_filter;
};

template <int N> // N is 2 or 3
struct cpt_ctx_t {
  int current_timestep;

  double *hV[2], *hGradV[2], *hScalar[2];
  double *dV[2], *dGradV[2], *dScalar[2];
};
