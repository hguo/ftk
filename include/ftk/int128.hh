#ifndef _FTK_INT128_HH
#define _FTK_INT128_HH

#include <ftk/config.hh>
#include <iostream>

namespace ftk {

#ifdef HAVE_int128_t_HAVE__int128_t
typedef __int128_t int128_t;
#elif defined(HAVE_int128_t_HAVEint128_t)
// We're okay
#elif defined(HAVE_int128_t_HAVE__int128)
typedef __int128 int128_t;
#elif defined(HAVE_int128_t_HAVEint128)
typedef int128 int128_t;
#else
#include <boost/multiprecision/cpp_int.hpp>
  typedef boost::multiprecision::int128_t int128_t;
#endif

#ifdef HAVE_uint128_t_HAVE__uint128_t
typedef __uint128_t uint128_t;
#elif defined(HAVE_uint128_t_HAVEuint128_t)
// We're okay
#elif defined(HAVE_uint128_t_HAVE__uint128)
typedef __uint128 uint128_t;
#elif defined(HAVE_uint128_t_HAVEuint128)
typedef uint128 uint128_t;
#elif defined(HAVE_uint128_t_HAVEunsigned__int128_t)
typedef unsigned __int128_t uint128_t;
#elif defined(HAVE_uint128_t_HAVEunsignedint128_t)
typedef unsigned int128_t uint128_t;
#elif defined(HAVE_uint128_t_HAVEunsigned__int128)
typedef unsigned __int128 uint128_t;
#elif defined(HAVE_uint128_t_HAVEunsignedint128)
typedef unsigned int128 uint128_t;
#else
#include <boost/multiprecision/cpp_int.hpp>
  typedef boost::multiprecision::uint128_t uint128_t;
#endif 


// https://stackoverflow.com/questions/25114597/how-to-print-int128-in-g
inline std::ostream& operator<<( std::ostream& dest, int128_t value )
{
  std::ostream::sentry s( dest );
  if ( s ) {
    uint128_t tmp = value < 0 ? -value : value;
    char buffer[ 128 ];
    char* d = std::end( buffer );
    do
    {
      -- d;
      *d = "0123456789"[ tmp % 10 ];
      tmp /= 10;
    } while ( tmp != 0 );
    if ( value < 0 ) {
      -- d;
      *d = '-';
    }
    int len = std::end( buffer ) - d;
    if ( dest.rdbuf()->sputn( d, len ) != len ) {
      dest.setstate( std::ios_base::badbit );
    }
  }
  return dest;
}

}

#endif
