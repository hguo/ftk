#ifndef _FTK_OBJECT
#define _FTK_OBJECT

#include <string>
#include "ftk/external/json.hh"

class ftkObject {
public:
  using json = nlohmann::json;

  virtual json toJson() const = 0; 
  virtual void fromJson(json j) = 0;

  std::string toJsonString(int indent=-1) const {
    return toJson().dump(indent);
  }

  void fromJsonString(const std::string &str) {
    auto j = json::parse(str);
    fromJson(j);
  }
};

#endif
