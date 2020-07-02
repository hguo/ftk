#ifndef _FTK_FIXED_INT_HH
#define _FTK_FIXED_INT_HH

namespace ftk {

template <int N> // stiching N integers together
struct fixed_integer {
  fixed_integer() {for (int i = 0; i < N; i ++) num[i] = 0;}
  fixed_integer(const std::string& str);
  fixed_integer(const fixed_integer& x) {
    for (int i = 0; i < N; i ++) 
      num[i] = x.num[i];
  }
  fixed_integer(int x) {I[0] = x;}

  fixed_integer& operator=(const fixed_integer& x) {
    for (int i = 0; i < N; i ++)
      num[i] = x.num[i];
    return *this;
  }
  fixed_integer& operator=(int x) {
    for (int i = 1; i < N; i ++)
      num[i] = 0;
    num[0] = x;
    return *this;
  }

  fixed_integer& operator+=(const fixed_integer& x) {
    for (int i = 0; i < N; i ++) {

    }
  }

private:
  int num[N];
};

} // namespace ftk

#endif
