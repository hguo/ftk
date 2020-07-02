
struct cpt_cfg_t {
  double threshold;
  unsigned int type_filter;
  bool use_explicit_coordinates;
  double *V[2], *gradV[2], *scalar[2];
  double *coords;
};

template <int N> // N is 2 or 3
struct cpt_ctx_t {
  int current_timestep;

  double *V[2], *gradV[2], *scalar[2];
};
