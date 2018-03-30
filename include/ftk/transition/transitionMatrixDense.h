#ifndef _FTK_TRANSITION_MATRIX_DENSE_H
#define _FTK_TRANSITION_MATRIX_DENSE_H

#include "ftk/transition/transitionMatrix.h"

class ftkTransitionMatrixDense : public ftkTransitionMatrix 
{
public:
  explicit ftkTransitionMatrixDense(int t0, int t1, int n0, int n1) : ftkTransitionMatrix(t0, t1, n0, n1) {
    _matrix.resize(n0*n1);
  }
  
  // int& operator()(int i, int j) {return _matrix[i*n1() + j];}
  // int operator()(int i, int j) const {return _matrix[i*n1() + j];}
  int& at(int i, int j) {return _matrix[i*n1() + j];}
  int at(int i, int j) const {return _matrix[i*n1() + j];}

private:
  std::vector<int> _matrix; // dense matrix
};

#endif
