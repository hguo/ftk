#ifndef _PUNCTURE_H
#define _PUNCTURE_H

#include <cereal/cereal.hpp>

struct punctured_face_t {
  float x[3];
  float cond;

  template<class Archive>
  void serialize(Archive & ar) {
    ar(
        cereal::make_nvp("x", x[0]), 
        cereal::make_nvp("y", x[1]), 
        cereal::make_nvp("z", x[2]), 
        cereal::make_nvp("cond", cond)
    );
  }
};

#endif
