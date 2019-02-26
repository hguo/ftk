#ifndef _FTK_EVENTS_H
#define _FTK_EVENTS_H

#include <vector>
#include <set>
#include <ftk/external/json.hh>

namespace ftk {

enum {
  FTK_EVENT_NONE = 0, // no event
  FTK_EVENT_BIRTH = 1,
  FTK_EVENT_DEATH = 2,
  FTK_EVENT_MERGE = 3,
  FTK_EVENT_SPLIT = 4,
  FTK_EVENT_RECOMBINATION = 5, 
  FTK_EVENT_COMPOUND = 6
};

template <class IdType, class LabelType>
struct Event {
  std::pair<IdType, IdType> interval;
  std::set<LabelType> lhs, rhs; // local ids on left and right hand sides

  int type() const {
    if (lhs.size() == 1 && rhs.size() == 1) { // no event
      return FTK_EVENT_NONE;
    } else if (lhs.size() == 0 && rhs.size() == 1) {
      return FTK_EVENT_BIRTH;
    } else if (lhs.size() == 1 && rhs.size() == 0) {
      return FTK_EVENT_DEATH;
    } else if (lhs.size() == 1 && rhs.size() == 2) {
      return FTK_EVENT_SPLIT;
    } else if (lhs.size() == 2 && rhs.size() == 1) { 
      return FTK_EVENT_MERGE;
    } else if (lhs.size() == 2 && rhs.size() == 2) { 
      return FTK_EVENT_RECOMBINATION;
    } else {
      return FTK_EVENT_COMPOUND;
    }
  }

  static std::string eventTypeToString(int e) {
    static const char* strs[7] = {
      "none", "birth", "death", "merge", "split", "recombination", "compound"};
    return strs[e];
  }
};

}

// json
namespace nlohmann {
  template <class IdType, class LabelType>
  struct adl_serializer<ftk::Event<IdType, LabelType> > {
    static void to_json(json& j, const ftk::Event<IdType, LabelType> &e) {
      j["interval"] = e.interval;
      j["type"] = ftk::Event<IdType, LabelType>::eventTypeToString(e.type());
      j["lhs"] = e.lhs;
      j["rhs"] = e.rhs;
    }

    static void from_json(const json& j, ftk::Event<IdType, LabelType> &e) {
      // TODO
    }
  };
}

#endif
