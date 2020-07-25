#ifndef _HYPERMESH_SYNTHETIC_DATA_H
#define _HYPERMESH_SYNTHETIC_DATA_H

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice.hh>

namespace ftk {

// the synthetic woven function
template <typename T> 
T woven_function_2Dt(T x, T y, T t) 
{
  return cos(x*cos(t)-y*sin(t))*sin(x*sin(t)+y*cos(t));
}

// generate 2D woven data
template <typename T>
ndarray<T> synthetic_woven_2D(int DW, int DH, T t = T(1e-4), T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH);

  // const T scaling_factor = 15; // the factor that controls the shape of the synthesize data
  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      const T x = ((T(i) / (DW-1)) - 0.5) * scaling_factor,
              y = ((T(j) / (DH-1)) - 0.5) * scaling_factor;
      scalar(i, j) = woven_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

// generate 2D woven data in a core
template <typename T>
ndarray<T> synthetic_woven_2D_part(const lattice& ext, const lattice& core, T t = T(1e-4), T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(core.sizes());

  const int DW = ext.size(0), DH = ext.size(1);
  for (int j = 0; j < core.size(1) ; j ++) {
    for (int i = 0; i < core.size(0); i ++) {
      const T x = ((T(i + core.start(0)) / (DW-1)) - 0.5) * scaling_factor,
              y = ((T(j + core.start(1)) / (DH-1)) - 0.5) * scaling_factor;
      scalar(i, j) = woven_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

// generate 2D time-varying woven data
template <typename T>
ndarray<T> synthetic_woven_2Dt(int DW, int DH, int DT, T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH, DT);

  // const T scaling_factor = 15; // the factor that controls the shape of the synthesize data
  for (int k = 0; k < DT; k ++) {
    for (int j = 0; j < DH; j ++) {
      for (int i = 0; i < DW; i ++) {
        const T x = ((T(i) / (DW-1)) - 0.5) * scaling_factor,
                y = ((T(j) / (DH-1)) - 0.5) * scaling_factor, 
                t = (T(k) / (DT-1)) + 1e-4;
        scalar(i, j, k) = woven_function_2Dt(x, y, t);
      }
    }
  }

  return scalar;
}

// double gyre 2D flow
template <typename T>
ndarray<T> synthetic_double_gyre(int DW, int DH, const T time, bool zchannel=false,
    const T A = 0.1, 
    const T omega = M_PI * 0.2, 
    const T epsilon = 0.25)
{
  const auto a = [&](T t) { return epsilon * sin(omega * t); };
  const auto b = [&](T t) { return 1 - 2 * epsilon * sin(omega * t); };
  const auto f = [&](T x, T t) { return a(t) * x * x + b(t) * x; };
  const auto dfdx = [&](T x, T t) {
    return 2 * a(t) * x * x + b(t);
  };
  const auto u = [&](T x, T y, T t) {
    return -M_PI * A * sin(M_PI * f(x, t)) * cos(M_PI * y);
  };
  const auto v = [&](T x, T y, T t) {
    return  M_PI * A * cos(M_PI * f(x, t)) * sin(M_PI * y) * dfdx(x, t);
  };

  ndarray<T> Vf;
  Vf.reshape(2 + zchannel, DW, DH);

  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      // the domain is [0, 2]x[0, 1]
      const T x = (T(i) / (DW-1)) * 2,
              y = (T(j) / (DH-1));

      Vf(0, i, j) = u(x, y, time);
      Vf(1, i, j) = v(x, y, time);
      // if (zchannel) Vf(2, i, j) = T(0);
    }
  }

  return Vf;
}

template <typename T>
ndarray<T> synthetic_time_varying_double_gyre(
    const int DW, const int DH, const int DT, 
    const T time_start, const T time_step, 
    const T A = 0.1, const T omega = M_PI * 0.2, const T epsilon = 0.25)
{
  ndarray<T> Vft;

  for (int t = 0; t < DT; t ++) {
    auto Vf = synthetic_double_gyre(DW, DH, 
        time_start + t * time_step, 
        A, omega, epsilon);
    Vft.p.insert(Vft.p.end(), Vf.p.begin(), Vf.p.end());
  }

  return Vft;
}

// ABC flow
template <typename T>
ndarray<T> synthetic_abc_flow(int DW, int DH, int DD, 
    T A=std::sqrt(T(3)), T B=std::sqrt(T(2)), T C=T(1))
{
  ndarray<T> Vf;
  Vf.reshape(3, DW, DH, DD);

  for (int k = 0; k < DD; k ++)
    for (int j = 0; j < DH; j ++)
      for (int i = 0; i < DW; i ++) {
        const T x = ((T(i) / (DW-1))) * 2 * M_PI,
                y = ((T(j) / (DH-1))) * 2 * M_PI,
                z = ((T(k) / (DD-1))) * 2 * M_PI;

        Vf(0, i, j, k) = A * sin(z) + C * cos(y);
        Vf(1, i, j, k) = B * sin(x) + A * cos(z);
        Vf(2, i, j, k) = C * sin(y) + B * cos(x);
      }

  return Vf;
}

// 2D merger
template <typename T>
T merger_function_2Dt(T x, T y, T t)
{
  auto f = [](T cx, T cy, T x, T y) {return exp(-((x-cx)*(x-cx) + (y-cy)*(y-cy)));};

  // add rotation
  T xp = x * cos(t) - y * sin(t), 
    yp = x * sin(t) + y * cos(t);
  x = xp;
  y = yp;

  T cx0 = sin(t - M_PI_2), // + 1e-4, 
    cx1 = sin(t + M_PI_2), // + 1e-4,
    cy0 = 1e-4,
    cy1 = 1e-4;
        
  return std::max(f(cx0, cy0, x, y), f(cx1, cy1, x, y));
}

template <typename T>
ndarray<T> synthetic_merger_2D(int DW, int DH, T t)
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH);
  
  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      // the domain is [-2, 2]x[-2, 2]
      const T x = ((T(i) / (DW-1)) - 0.5) * 4,
              y = ((T(j) / (DH-1)) - 0.5) * 4;
      scalar(i, j) = merger_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

template <typename T, int N>
ndarray<T> synthetic_moving_extremum(const std::vector<size_t>& shape, const T x0[N], const T dir[N], T t)
{
  ndarray<T> scalar(shape);
  const auto lattice = scalar.get_lattice();

  T xc[N]; // center
  for (int j = 0; j < N; j ++)
    xc[j] = x0[j] + dir[j] * t;

  for (auto i = 0; i < scalar.nelem(); i ++) {
    std::vector<int> xi = lattice.from_integer(i);
    T x[N];
    T d = 0;
    for (int j = 0; j < N; j ++) {
      d += pow(xi[j] - xc[j], T(2.0));
    }
    scalar[i] = d;
  }

  return scalar;
}

/// modified from https://web.cse.ohio-state.edu/~crawfis.3/Data/Tornado/tornadoSrc.c
// void gen_tornado( int xs, int ys, int zs, int time, float *tornado )
/*
 *  Gen_Tornado creates a vector field of dimension [xs,ys,zs,3] from
 *  a proceedural function. By passing in different time arguements,
 *  a slightly different and rotating field is created.
 *
 *  The magnitude of the vector field is highest at some funnel shape
 *  and values range from 0.0 to around 0.4 (I think).
 *
 *  I just wrote these comments, 8 years after I wrote the function.
 *  
 * Developed by Roger A. Crawfis, The Ohio State University
 *
 */
template <typename T>
ndarray<T> synthetic_tornado(int xs, int ys, int zs, int time)
{
  ndarray<T> array;
  array.reshape({3, static_cast<size_t>(xs), static_cast<size_t>(ys), static_cast<size_t>(zs)});
  T *tornado = &array[0];

  T x, y, z;
  int ix, iy, iz;
  T r, xc, yc, scale, temp, z0;
  T r2 = 8;
  T SMALL = 0.00000000001;
  T xdelta = 1.0 / (xs-1.0);
  T ydelta = 1.0 / (ys-1.0);
  T zdelta = 1.0 / (zs-1.0);

  for( iz = 0; iz < zs; iz++ )
  {
     z = iz * zdelta;                        // map z to 0->1
     xc = 0.5 + 0.1*sin(0.04*time+10.0*z);   // For each z-slice, determine the spiral circle.
     yc = 0.5 + 0.1*cos(0.03*time+3.0*z);    //    (xc,yc) determine the center of the circle.
     r = 0.1 + 0.4 * z*z + 0.1 * z * sin(8.0*z); //  The radius also changes at each z-slice.
     r2 = 0.2 + 0.1*z;                           //    r is the center radius, r2 is for damping
     for( iy = 0; iy < ys; iy++ )
     {
		y = iy * ydelta;
		for( ix = 0; ix < xs; ix++ )
		{
			x = ix * xdelta;
			temp = sqrt( (y-yc)*(y-yc) + (x-xc)*(x-xc) );
			scale = fabs( r - temp );
/*
 *  I do not like this next line. It produces a discontinuity 
 *  in the magnitude. Fix it later.
 *
 */
           if ( scale > r2 )
              scale = 0.8 - scale;
           else
              scale = 1.0;
			z0 = 0.1 * (0.1 - temp*z );
		   if ( z0 < 0.0 )  z0 = 0.0;
		   temp = sqrt( temp*temp + z0*z0 );
			scale = (r + r2 - temp) * scale / (temp + SMALL);
			scale = scale / (1+z);
           *tornado++ = scale * (y-yc) + 0.1*(x-xc);
           *tornado++ = scale * -(x-xc) + 0.1*(y-yc);
           *tornado++ = scale * z0;
		}
     }
  }

  return array;
}

}

#endif
