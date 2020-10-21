#ifndef _GLHEADER_H
#define _GLHEADER_H

enum {
  DTYPE_BDAT, 
  DTYPE_CA02
};

typedef struct {
  int ndims; 
  int dims[3];
  bool pbc[3];
  float zaniso;
  float lengths[3], origins[3], cell_lengths[3];
  float time;
  float B[3];
  float Jxext, Kex, Kex_dot, V;
  float fluctuation_amp;
  int dtype;
} GLHeader;

#endif
