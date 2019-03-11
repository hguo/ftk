#ifndef _FTK_RAND_HH
#define _FTK_RAND_HH

#include <numeric>
#include <random>

namespace ftk {

template <typename T, int n>
void rand(T V[n])
{
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  for (int i = 0; i < n; i ++)
    V[i] = d(gen);
}

template <typename T, int m, int n>
void rand(T M[m][n])
{
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  for (int i = 0; i < m; i ++)
    for (int j = 0; j < n; j ++)
      M[i][j] = d(gen);
}

template <typename T, int m, int n>
void rand_symmetric(T M[m][n])
{
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  for (int i = 0; i < m; i ++)
    for (int j = i; j < n; j ++)
      M[i][j] = M[j][i] = d(gen);
}

template <typename T>
void rand2(T V[2])
{
  rand<T, 2>(V);
}

template <typename T>
void rand3(T V[3])
{
  rand<T, 3>(V);
}

template <typename T>
void rand2x2(T M[2][2])
{
  rand<T, 2, 2>(M);
}

template <typename T>
void rand_symmetric2x2(T M[2][2])
{
  rand_symmetric<T, 2, 2>(M);
}

template <typename T>
void rand2x3(T M[2][3])
{
  rand<T, 2, 3>(M);
}

template <typename T>
void rand3x3(T M[3][3])
{
  rand<T, 3, 3>(M);
}

template <typename T>
void rand3x2(T M[3][2])
{
  rand<T, 3, 2>(M);
}

template <typename T>
void rand4x3(T M[4][3])
{
  rand<T, 4, 3>(M);
}

template <typename T>
void rand_symmetric3x3(T M[3][3])
{
  rand_symmetric<T, 3, 3>(M);
}

}

#endif
